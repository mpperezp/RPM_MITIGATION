#*********************************************************************************************#
#                            PRECIPITATION MITIGATION FACTOR FOR PAVED RAODS                  #
#              Written by: Maria Paula Perez P.                                               #
#                          Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs the mitigation caused by rain events on paved  #
#                          road emissions                                                     #
#              References: EPA - AP42 Chapter 13 Miscellaneous sources                        #
#*********************************************************************************************#

import os
import datetime
import numpy as np
from netCDF4 import Dataset

max_hours_since_precip = 12

def creditsf(ISRAIN, HSP):
    """
    Requires
        ISRAIN = vector of booleans rain is higher than threshold
        HSP = vector of hours since precipitation
    
    Returns
        CREDITS = vector of migitation credits for application to future hours
    """    
    out = np.zeros(ISRAIN.shape, dtype = 'i')
    if (ISRAIN == False).all():
        return out
    for i, israin in enumerate(ISRAIN):
        if i == 0 or HSP[i] > max_hours_since_precip:# or (ISRAIN[i - 1] > 0 and israin < 0):
            out[i] = israin.astype('i')
        else:
            out[i] = israin + out[i - 1]
    return out * 0.2


def nmf(rain, threshold = 0.0254, verbose = False):
    """
    Requires
        rain = vector of hourly precipitation amounts (units assumed to be cm)
    Optional
        threshold = minimum rain for dust suppression in units consistent with rain
        verbose = set to true to print intermediate results

    Returns
        nmf = natural mitigation factor to be applied to dust emissions
    """
    israin = (rain > threshold)
    if (israin == False).all():
        return np.ones(israin.shape, dtype = 'f')
    nat_mit_factor = 1 - 1.2 * israin

    hours_since_rain = -np.ma.masked_greater((np.where(israin)[0][:, None] - np.indices(israin.shape)[0]), 0).max(0)
    debit_eligible = israin == False

    credits = np.round(creditsf(israin, hours_since_rain), 1)
    debits = np.round((hours_since_rain * 0.2).filled(1), 1)

    nat_mit_factor[debit_eligible & (credits >= debits)] -= 0.2
    nat_mit_factor = np.ma.masked_less(nat_mit_factor, 0).filled(0)
    if verbose:
        np.set_printoptions(precision = 2, formatter = dict(float = lambda x: '%0.2f' % x))
        print 'Rain'
        print rain
        print 'Credits'
        print credits
        print 'Debits'
        print debits
        print 'NMF'
        print nat_mit_factor

    return nat_mit_factor


# Create a file for storage
def getoutfile(metpath, outpath, ntimes):
    """
    Requires
        metpath - path to METCRO2D file
        outpath - path for output file
    Returns
        open NetCDF output file
    """
    os.system('ncks -O -dVAR,0 -vRN,TFLAG  %s mettemp.nc' % metpath)
    os.system('ncdump -vTFLAG mettemp.nc > mettemp.cdl')
    os.system('sed -E -e \'s/VAR-LIST = .*/VAR-LIST = "NMF             ";/\' -e "s/RN/NMF/g" -e \'s/units = "CM      /units = "fraction/g\' -e \'s/nonconvec. pcpn per met TSTEP/natural mitigation factor/g\' mettemp.cdl > nmf.cdl')
    os.system('ncgen -o %s  nmf.cdl' % outpath)
    outf = Dataset(outpath, 'a')
    outf.NVARS = 1
    tflag = outf.variables['TFLAG']
    sdate = outf.SDATE
    stime = outf.STIME
    datestart = datetime.datetime.strptime('%d %06d' % (sdate, stime), '%Y%j %H%M%S')
    for ti in range(ntimes):
        #datenow = datetime.time(hours = ti)
        datenow = datestart + datetime.timedelta(hours = ti)
        tflag[ti, 0, :] = int(datenow.strftime('%Y%j')), int(datenow.strftime('%H%M%S'))
        #tflag[ti, 0, :] = int(datenow.strfdelta('%Y%j')), int(datenow.strfdelta('%H%M%S'))
    outf.sync()    
    return outf 

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog='python %s' % sys.argv[0])
    parser.add_argument('METCRO2D', help='path to METCRO2D file; used for file structure (25 hours)')
    parser.add_argument('RAININPUT', help='path to RAIN input file; used for file data (nday * 24 + 1 hours')
    args = parser.parse_args()
    
    metpath = args.METCRO2D
    rainpath = args.RAININPUT
    rainf = Dataset(rainpath)
    rain = rainf.variables['RN'][:] # + rainf.variables['RC'][:]

    # Use the nmf function to calculate NMF for each time vector within rain
    nmf4d = np.apply_along_axis(nmf, axis = 0, arr = rain)
    np.set_printoptions(precision = 2, formatter = dict(float = lambda x: '%0.2f' % x))
#print nmf4d

    # Check results for a random cell
    print 'Rain (mm)'
    print rain[:, 32, 32] * 10
    print 'NMF'
    print nmf4d[:, 32, 32]

    # Open the empty file and put results in it.
    nmff = getoutfile(metpath, 'NMF.nc', nmf4d.shape[0])
    nmff.variables['NMF'][:, 0] = nmf4d
    nmff.sync()
