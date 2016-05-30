
#*********************************************************************************************#
#                            TRANSPORTABLE FRACTION ESTIMATION SCRIPT                         #
#              Written by: Maria Paula Perez P.						      #
#			   Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs the transportable fraction considered as the   #
#                          amount of fugitive dust that has not been affected by the capture  # 
#                          fraction - removal factor from local surfaces (e.g. vegetation,    # 
#                          buildings) of re-suspended particulate matter -                    #
#              References: Pouliot, el.al.,(2009) available at                                #
#                          http://www.epa.gov/ttnchie1/conference/ei19/session9/pouliot.pdf   #
#*********************************************************************************************#

# Imports
import netCDF4
from netCDF4 import Dataset
import numpy as np
import os
import datetime

# Values of transportable fraction per land use available in Pouliot, et.al.,(2009) Table 1: BELD3 categories, Capture Fraction Class, and Transportable Fraction
TFF = 0.05   #Trans Frac for Forest USGS
TFU = 0.5    #Trans Frac for Urban USGS
TFGS = 0.75  #Trans Frac for Grassland
TFA = 0.75   #Trans Frac for Agricultural USGS
TFBW = 1.0   #Trans Frac for Barren/Water USGS

# Function to estimate transportable fraction
def trans_frac(landpath):
    """
    Requires
        land use - file of GRIDCRO2D with fraction of land use per grid cell USGS land use type
    Returns 
        TF - Total Transportable Fraction for the modeling domain
    """
    landpathf = Dataset(landpath)
    # USGS land use types selected from GRIDCRO2D, add classes if neccesary
    urban = landpathf.variables['LUFRAC_01'][:]
    drycroppst = landpathf.variables['LUFRAC_02'][:]
    cropgrassmos = landpathf.variables['LUFRAC_05'][:]
    grassland = landpathf.variables['LUFRAC_07'][:]
    shrubland = landpathf.variables['LUFRAC_08'][:]
    savanna = landpathf.variables['LUFRAC_10'][:]
    evgrnbrdfor = landpathf.variables['LUFRAC_13'][:]
    mixfor = landpathf.variables['LUFRAC_15'][:]
    woodtndr = landpathf.variables['LUFRAC_21'][:]
    
    URTF = urban * TFU
    GRTF = grassland * TFGS + cropgrassmos * TFGS + drycroppst * TFGS + savanna * TFGS  + woodtndr * TFGS
    WTBR = shrubland * TFBW
    FRST = evgrnbrdfor * TFF + mixfor * TFF
    TF = URTF + GRTF + WTBR + FRST 
    return TF

# Create a netCDF file to store the Transportable Fraction
def getoutfile(landpath, outpath, ntimes):
    """
    Requires
	landpath - path to GRIDCRO2D file
        outpath - path for output file
    Returns
        open NetCDF output file
    """
    os.system('ncks -O -dVAR,0 -vLUFRAC_01,TFLAG  %s trtemp.nc' % landpath)
    os.system('ncdump -vTFLAG trtemp.nc > trtemp.cdl')
    os.system('sed -E -e \'s/VAR-LIST = .*/VAR-LIST = "TF             ";/\' -e "s/LUFRAC_01/TF/g" -e \'s/units = "FRACTION      /units = "FRACTION/g\' -e \'s/USGS24: Urban Land/transportable fraction/g\' trtemp.cdl > tr_fra.cdl')
    os.system('ncgen -o %s  tr_fra.cdl' % outpath)
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

# Opening GRIDCRO2D file with land use type fraction to calculate TF
if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog='python %s' % sys.argv[0])
    parser.add_argument('GRIDCRO2D', help = 'format file with land use fraction')
    args = parser.parse_args()

    landpath = args.GRIDCRO2D
    

    # Call the function
    transp = trans_frac(landpath)

    # Check the result for a random grid cell
    print transp[:, :, 32, 32]

    # Open the empty file and place results in it.
    transpfraction = getoutfile(landpath, 'TRANSFRAC.nc', transp.shape[0])
    transpfraction.variables['TF'][:, 0] = transp
    transpfraction.sync()
