#*********************************************************************************************#
#                            PRECIPITATION MITIGATION FACTOR FOR UNPAVED ROADS                #
#              Written by: Maria Paula Perez P.                                               #
#                          Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs the mitigation caused by rain events on unpaved#
#                          road emissions                                                     #
#              References: EPA - AP42 Chapter 13 Miscellaneous sources                        #
#*********************************************************************************************#


import netCDF4
from netCDF4 import Dataset
import os
import numpy as np
from numpy import array
import datetime


def nmf_precipitation(rain, threshold = 0.0254):
    """
    Requires
        rain = vector of hourly precipitation amounts (units assumed to be cm)
        threshold = minimum rain for dust suppression in units consistent with rain
    Returns
        nmf_precipitation = natural mitigation factor due to precipitation occurance per hour
    """
    #for RN in rain.variables['RN'][:]:
    #    if RN.any() > threshold:
    #        nmfp = 0.2
    #    else:
    #        nmfp = 1
    nmfp = np.where(rain.variables['RN'][:] > threshold, 0.2, 1.0)
    return nmfp

# Create a netCDF file to store the NMF per day
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
        datenow = datestart + datetime.timedelta(hours = ti)
        tflag[ti, 0, :] = int(datenow.strftime('%Y%j')), int(datenow.strftime('%H%M%S'))
    outf.sync()
    return outf

# Opening file with hourly rain RN
if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog='python %s' % sys.argv[0])
    parser.add_argument('METCRO2D', help='path to file meteorology format')
    args = parser.parse_args()
    
    metpath = args.METCRO2D
    met_inn = Dataset(metpath)

    # Call the function
    nmfpi = nmf_precipitation(met_inn)
   
    # Open the empty file and put results in it.
    nmff = getoutfile(metpath, 'UNPAVED_NMF', nmfpi.shape[0])
    nmff.variables['NMF'][:, 0] = nmfpi
    nmff.sync() 
    
    # Check the outputs for a cell
    #np.set_printoptions(precision = 4)
    #print met_inn.variables['RN'][:][7, :, 32, 32]
    #print nmfpi[7, :, 32, 32]
