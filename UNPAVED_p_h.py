#*********************************************************************************************#
#                           TOTAL MITIGATION FACTOR FOR UNPAVED ROADS (p(x,t))&(h(x,t))       #
#              Written by: Maria Paula Perez P.                                               #
#                          Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs the mitigation caused by rain r(x,t)  events   #
#                          and soil moisture m(x,t) on unpaved roads                          #
#              References: EPA - AP42 Chapter 13 Miscellaneous sources                        #
#*********************************************************************************************#


import netCDF4
from netCDF4 import Dataset
import os
import numpy as np
from numpy import array
import datetime

#*****************Total Natural Mitigation Factor for UNPAVED ROADS*****************

def nat_mit_fact(nmf_precipitation, nmf_moisture):
    """

    """
#    import pdb; pdb.set_trace()
    NMFf = nmf_precipitation.variables['NMF'][:] * nmf_moisture.variables['NMFM'][:]
    #print 'Nat mit fac (prec*moist)', NMFf[1, :, 32, 32]    # Check functions performance
    return NMFf

# Create a netCDF file to store the NMF per day
def getoutfile(path, outpath, ntimes):
    """
    Requires
        path - path to nmfp file to give structure
        outpath - path for output file
    Returns
        open NetCDF output file
    """
    os.system('ncks -O -dVAR,0 -vNMFM,TFLAG  %s mettemp.nc' % path)
    os.system('ncdump -vTFLAG mettemp.nc > mettemp.cdl')
    os.system('sed -E -e \'s/VAR-LIST = .*/VAR-LIST = "NMFM             ";/\' -e "s/NMFM/NMF/g" -e \'s/units = "fraction      /units = "fraction/g\' -e \'s/natural mitigation factor due to soil moisture fraction/total natural mitigation factor due to precipitation and moisture/g\' mettemp.cdl > nmf.cdl')
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

# Opening files with mitigation effect due to precipitation and mitigation effect due to surface moisture content
if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog='python %s' % sys.argv[0])
    parser.add_argument('PREC_NMFIN', help = 'file with the nat. mit. fact. due to precipitation')
    parser.add_argument('SOIM_NMFIN', help = 'file with the nat. mit. fact. due to surface moisture')
    args = parser.parse_args()

    nmfp_path = args.PREC_NMFIN
    nmfm_path = args.SOIM_NMFIN
    nmfpv = Dataset(nmfp_path)
    nmfmv = Dataset(nmfm_path)

    # Call the functions
    #soilmoist_correc(nmfmv)    
    NMFfin = nat_mit_fact(nmfpv, nmfmv)
    
    # Open the empty file and put results in it.
    nmff = getoutfile(nmfm_path, 'UNPAVED_NMF.nc', NMFfin.shape[0])
    nmff.variables['NMF'][:, 0] = NMFfin
    nmff.sync()

    # Check functions performance
    np.set_printoptions(precision = 7)
    #print nmff.variables['NMF'][:][1, :, 32, 32]
    print 'Nat mit fact due to surface moisture h(x,t) at cell 32,32', nmfmv.variables['NMFM'][:][1, :, 32, 32]
    print 'Nat mit fact due to precipitation p(x,t) at cell 32,32', nmfpv.variables['NMF'][:][1, :, 32, 32]
