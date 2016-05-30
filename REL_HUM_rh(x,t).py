#*********************************************************************************************#
#                          RELATIVE HUMIDITY CALCULATION		                      #
#              Written by: Maria Paula Perez P.                                               #
#                          Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs relative humidity neccesary to use AP42_2.py   #
#                          Works on-line with the code                                        #
#*********************************************************************************************#

import netCDF4
from netCDF4 import Dataset
import os
import numpy as np
from numpy import array
import datetime

#**********RELATIVE HUMIDITY CALCULATION TO ESTIMATE SURFACE MOISTURE FRACTION************

# Constants to calculate Relative Humidity from Water vapor mixing Ratio (Q2) and surface Pressure (PRSFC) and Temperature at 2m from NCL formula 

A2 = 17.2693882
A3 = 273.16
A4 = 35.86
PQ0 = 379.90516

def relhum(metpath):
    """
    Requires
        metpath - File of METCRO2D with all meteorology variables
    Returns
        R - Relative Humidity %
    """
    R = metpath.variables['Q2'][:] / ((PQ0 / metpath.variables['PRSFC'][:]) * np.exp(A2 * (metpath.variables['TEMP2'][:] - A3) / (metpath.variables['TEMP2'][:] - A4)))
    return R

# Create a netCDF file to store the Relative Humidity per day and hour
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
    os.system('sed -E -e \'s/VAR-LIST = .*/VAR-LIST = "HREL             ";/\' -e "s/RN/HREL/g" -e \'s/units = "CM      /units = "percentage/g\' -e \'s/nonconvec. pcpn per met TSTEP/relative humidity/g\' mettemp.cdl > rel_hum.cdl')
    os.system('ncgen -o %s  rel_hum.cdl' % outpath)
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
        
# Opening Meteorology file with variables required to calculate relative humidity
if __name__ == '__main__':
    import argparse 
    import sys 
    parser = argparse.ArgumentParser(prog='python %s' % sys.argv[0])
    parser.add_argument('METCRO2D', help = 'format file with structure to the Relative Humidity and with values to calculate it')
    args = parser.parse_args()
    
    metpath = args.METCRO2D
    metpathf = Dataset(metpath)    

    # Call the function
    relf = relhum(metpathf)
    #print relf[:, :, 32, 32]
    
    # Open the empty file and place results in it.
    rel_h = getoutfile(metpath, 'HREL.nc', relf.shape[0])
    rel_h.variables['HREL'][:, 0] = relf
    rel_h.sync()
