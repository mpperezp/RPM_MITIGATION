#*********************************************************************************************#
#                      SURFACE HUMIDITY FACTOR h(x,t) 			                      #
#              Written by: Maria Paula Perez P.                                               #
#                          Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs the mitigation caused by surface humidity on   #
#                          unpaved road emissions                                             #
#              References: EPA - AP42 Chapter 13 Miscellaneous sources                        #
#*********************************************************************************************#

import netCDF4
from netCDF4 import Dataset
import os
import numpy as np
from numpy import array
import datetime

def soilmoist_correc(soil_mass_frac):
    """
    Requires
        soil_perc = netCDF file with SOIM1 as a mass percentage
    Returns
        soilmois = soil moisture correction term
    """
    nmfm = ((soil_mass_frac.variables['SFMOIS'][:] / 0.5) ** -0.2)
    #print 'Soil moisture correction term M', nmfm[:, :, 32, 32]    # Check functions performance
    return nmfm

# Create a netCDF file to store the NMF per day
def getoutfile(sfpath, outpath, ntimes):
    """
    Requires
        sfpath - path to SFMOIS file 
        outpath - path for output file
    Returns
        open NetCDF output file
    """
    os.system('ncks -O -dVAR,0 -vSFMOIS,TFLAG  %s mettemp.nc' % sfpath)
    os.system('ncdump -vTFLAG mettemp.nc > mettemp.cdl')
    os.system('sed -E -e \'s/VAR-LIST = .*/VAR-LIST = "NMFM             ";/\' -e "s/SFMOIS/NMFM/g" -e \'s/units = "percentage      /units = "fraction/g\' -e \'s/surface moisture/natural mitigation factor due to soil moisture/g\' mettemp.cdl > nmf.cdl')
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
    parser.add_argument('SOIM_NMFIN', help = 'file with the soil moisture mass percentage to calculate the nat. mit. fact. due to soil moisture')
    args = parser.parse_args()

    nmfm_path = args.SOIM_NMFIN
    nmfmv = Dataset(nmfm_path)

    # Call the functions
    #soilmoist_correc(nmfmv)    
    NMFM = soilmoist_correc(nmfmv)

    # Open the empty file and put results in it.
    nmff = getoutfile(nmfm_path, 'NMF_MOIS.nc', NMFM.shape[0])
    nmff.variables['NMFM'][:, 0] = NMFM
    nmff.sync()

