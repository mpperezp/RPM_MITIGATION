#*********************************************************************************************#
#                         CALCULATION OF SURFACE MOISTURE CONTENT IN MASS %                   #
#              Written by: Maria Paula Perez P.                                               #
#                          Barron Henderson, PhD.                                             #
#                   email: mpperezpena@hotmail.com                                            #
#             Description: This script outputs surface moisture m(x,t) in % to use in         #
#                          the calculation of h(x,t)                                          #
#              References: EPA - AP42 Chapter 13 Miscellaneous sources, excel sheet           #
#*********************************************************************************************#

import netCDF4
from netCDF4 import Dataset
import numpy as np
from numpy import array
import os
import datetime
from REL_HUM import relhum

#*************************Surface Moisture for Unpaved Roads***********************

#*********Calculation (Based on AP42 - spreadsheet Chapter 13.2.2 - http://www.epa.gov/ttn/chief/ap42/ch13/related/c13s02-2.html) 

#CHANGE VARIABLES WHEN NECCESARY
#**** Variables needed to calculate surface moisture fraction
TDF = 1.0        # Threshold level for dew point to form
HRNtoRN = 0.0254 # Relationship between threshold value and Relative Humidity units = cm
ro = 1440.0      # Road surface density units = kg/m**3 (Appendix A - AP42 Densities of Selected Substances)
T = 0.25         # Road surface thickness
MAXMp = 20       # Maximum moisture percentage
MAXM = 0.1786    #1.125    # Maximum moisture units = cm
MINMp = 0.2      # Minimum moisture percentage
AHE = 12./8760   # Average Hourly Evaporation (month = 12 / 8760)
EC = 0.75 * AHE  # Evaporation Constant
EM = 7.60        # Evaporation per month units = cm Average for Bogota for both February and October (IDEAM, 2005)
VH = 19.18   # Average Vehicle per Hour per lane for 16 unpaved measured stations
HOURLY_EVAP = ((EC * EM * VH) / MAXM) # hourly water loss in fraction of maximum water in cm
Mp_LOWER_BOUND = MINMp / MAXMp        # minimum percent of maximum moisture
#PHM Previous hour surface moisture

#********Function of surface moisture
def surf_moist(rain25, relhum25, soilmois25, spin = None):
    """
    Requires
        rain - path to file with RN in cm
        re_humidity - path to file with Relative Humidity
    Returns
        surface moisture - hourly value of surface moisture as a fraction  
    """
    surfmois25 = np.zeros(rain25.shape)
    for hour_idx in range(25):
        rain = rain25[hour_idx]
        relhum = relhum25[hour_idx]
        surfmois = surfmois25[hour_idx]
        if hour_idx == 0:
            if spin is None:
                PHM = soilmois25[hour_idx] / 1.44 * 100 / MAXMp #Calculate the surface moisture as a mass percentage from the SOIM1 [M**3/M**3] volumetric soil moisture using the soil density
                PHM[PHM > 1] = 1
                PHM[PHM < Mp_LOWER_BOUND] = Mp_LOWER_BOUND
            else:
                PHM = spin[:][-2] / MAXMp
                assert(round(PHM.min(), 2) >= Mp_LOWER_BOUND)
                assert(round(PHM.max(), 2) <= 1)
        else:
            PHM = surfmois25[hour_idx - 1]
        #print hour_idx, PHM.max()
        not_raining = (rain[:] == 0.0) 
        rh_condensing = (relhum[:] >= TDF)
        rh_condensing_saturates = (PHM * MAXM + HRNtoRN) > MAXM
        rain_saturates = ((rain[:] + (PHM * MAXM)) > MAXM)
        evap_desicates = (PHM - HOURLY_EVAP) < Mp_LOWER_BOUND #Is the Pevious Hours Moisture minus the Evaporation greater than the minimum moisture level for the surface material?
    #********************************If the RAIN IS equal to 0.0***************************
        #print 'Start', 'mask'
        rhsat = not_raining&rh_condensing&rh_condensing_saturates
        rhunsat = not_raining&rh_condensing&(~rh_condensing_saturates)
        evaponly = not_raining&(~rh_condensing)&(~evap_desicates)
        evapmin = not_raining&(~rh_condensing)&(evap_desicates)
        rainsat = (~not_raining)&rain_saturates
        rainunsat = (~not_raining)&(~rain_saturates)
        surfmois[rhsat] = 1.0 #If the rain is equal to 0 and the RELHUM IS greater than the threshold for DEW TO FORM and the addition would saturate surface material, THEN material is saturated 100%
        #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
        #if len(check) > 0: print hour_idx, 'RH saturation', np.percentile(check, [0, 10, 25, 50, 75, 90, 100]), check.mean()
        surfmois[rhunsat] = (HRNtoRN / MAXM) + PHM[rhunsat] #If rain is equal to 0 and the RELHUM IS greater than the threshold for DEW TO FORM and the adittion would not saturate surface material, THEN add the water to the surface material
        #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
        #if len(check) > 0: print hour_idx, 'RH unsaturation', np.percentile(check, [0, 10, 25, 50, 75, 90, 100]), check.mean()
        surfmois[evaponly] = PHM[evaponly] - HOURLY_EVAP #If the rain is equal to 0 and RELHUM is NOT greater than threshold, and the previous hours moisture minus the evaporation is GREATER than the minimum moisture level for the surface material, THEN reduce the moisture level
        #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
        #if len(check) > 0: print hour_idx, 'EVAP only', np.percentile(check, [0, 10, 25, 50, 75, 90, 100]), check.mean()
        surfmois[evapmin] = Mp_LOWER_BOUND  #If the rain is equal to 0 and RELHUM is NOT greater than threshold, and the previous hours moisture minus the evaporation is LESS than the minimum moisture level for the surface material, THEN set moisture level at minimum
        #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
        #if len(check) > 0: print hour_idx, 'setmin', np.percentile(check, [0, 10, 25, 50, 75, 90, 100]), check.mean()

    #******************************If the RAIN is NOT equal to 0.0***************************
        surfmois[rainsat] = 1.0 #If it rained and the addition of the amount of rain would exceed the maximum water that the surface material could hold, THEN set surface moisture at maximum 100%
        #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
        #if len(check) > 0: print hour_idx, 'rain sat', np.percentile(check, [0, 10, 25, 50, 75, 90, 100]), check.mean()
        rain_split = rain[:]/MAXM
        surfmois[rainunsat] = PHM[rainunsat] + rain_split[rainunsat] #If it rained but the addition of the amount of rain would not exceed the maximum water that the surface material could hold, THEN add the rain inpunt to the previous surface moisture    
        #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
        #if len(check) > 0: print hour_idx, 'rain unsat', np.percentile(check, [0, 10, 25, 50, 75, 90, 100]), check.mean()
        conditions = rhsat.astype('i') + rhunsat + evaponly + evapmin + rainsat + rainunsat
        assert(conditions.sum() == np.prod(conditions.shape))
        #print hour_idx, conditions.min(), conditions.max()
    assert(round(surfmois25.min(), 2) >= Mp_LOWER_BOUND)
    assert(round(surfmois25.max(), 2) <= 1)
    #check = np.ma.compressed(np.ma.masked_values(surfmois, 0))
    #if len(check) > 0: print 'Day', 'rain unsat', np.percentile(check, [0, 10, 25, 50, 75, 90, 100])
    
    retval = surfmois25 * MAXMp

#    meantssm = retval.reshape(-1, np.prod(retval.shape[1:])).mean(1)
#    nmfm = (retval/.5)**(-.2)   # Moisture term NMFm = (M/0.5)^-0.2
#    meannmfm = nmfm.reshape(-1, np.prod(nmfm.shape[1:])).mean(1)
#    print 'SM', meantssm
#    print 'NMFm', meannmfm

    return retval

# Create a netCDF file to store the Surface moisture fraction per day and hour
def getoutfile(metpath, outpath, ntimes):
    """
    Requires
        metpath - path to METCRO2D file
        outpath - path for output file
    Returns
        open NetCDF output file
    """
    os.system('ncks -O -dVAR,0 -vRN,TFLAG  %s sfmtemp.nc' % metpath)
    os.system('ncdump -vTFLAG sfmtemp.nc > sfmtemp.cdl')
    os.system('sed -E -e \'s/VAR-LIST = .*/VAR-LIST = "SFMOIS             ";/\' -e "s/RN/SFMOIS/g" -e \'s/units = "CM      /units = "percentage/g\' -e \'s/nonconvec. pcpn per met TSTEP/surface moisture/g\' sfmtemp.cdl > srfc_mois.cdl')
    os.system('ncgen -o %s  srfc_mois.cdl' % outpath)
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

# Open files with information to calculate surface moisture

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog='python %s' % sys.argv[0])
    parser.add_argument('METCRO2D', help = 'format file with Meteorological parameters to calculate the surface moisture and to give structure to the Soil Moisture containig file')
    #parser.add_argument('RELHUM_IN', help = 'file containing reltative humidity calculated from water vapor mixing ratio')
    parser.add_argument('SPIN', nargs='?', help = 'optional file to spin-up the surface moisture')
    args = parser.parse_args()    
    
    metpath = args.METCRO2D
    #relhumpath = args.RELHUM_IN
    #relhump = Dataset(relhumpath)
    #relhum = relhump.variables['HREL'][:]
    metfile = Dataset(metpath)
    relhum = relhum(metfile)
    rain = metfile.variables['RN'][:]
    soilm1 = metfile.variables['SOIM1'][:]
    if not args.SPIN is None:
        spinpath = args.SPIN
        spinfile = Dataset(spinpath)
        spin = spinfile.variables['SFMOIS'][:]
    else:
        spin = None

    # Call the function to calculate surface moisture
    srfcmois = surf_moist(rain, relhum, soilm1, spin)  

    # Open the empty file and place results in it.
    surf_m = getoutfile(metpath, 'SFMOIS.nc', srfcmois.shape[0])
    surf_m.variables['SFMOIS'][:, 0] = srfcmois
    surf_m.sync()

   # Check values for cells
    #print 'rain', rain[:]
    #print 'relhum', relhum[:]
    #print 'surfmois', srfcmois[:]
    
