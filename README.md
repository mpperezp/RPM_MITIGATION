# RPM_MITIGATION_FACTORS
Python scripts to develop mitigation fraction for RPM emissions from paved and unpaved roads using MCIP outputs
SCRIPT FOLDER                                     
       Authors: Maria Paula Perez P.                                                   
                Barron Henderson, PhD                                                  
       email : mpperezpena@hotmail.com                                                 
       Period: February 2016                                                           
       Contents: This folder contains the scripts developed to estimate Natural        
       Mitigation Factors for paved and unpaved roads, codes require MCIP ouputs       
       (METCRO2D and GRIDCRO2D), they work under any TSTEP, LAY, ROW, COL              
       within the AQM model set-up                                                     
       Last update: April 2016                                                         
___________
|--FOLDERS
___________

          ****THE SCRIPTS ARE LISTED IN THE SEQUENTIAL ORDER AS SHOULD BE USED
        *RQ = input requirements
        *OP = output

-----------------
|-- Paved_Roads
-----------------
        1. PAVED_NMF.py - Python script to calculate the NMF for paved roads considering moisture credits following cessation of precipitation to a maximum up to 12 hours after a rain event
                1. RQ:   METCRO2D - File to give structure to the NMF file
                         RAIN - File with the precipitation in cm for the entire modeling period (Precipitation variables from MCIP are RN and RC, for this case only RN is used because WRF is expected to solve RC within RN on the scale used in this work (d04) and passed it to MCIP). This file is obtained by extracting RN variable from METCRO2D files available for the modeling period and concatenating them using ($cdo copy $(for m in METCRO2D_onlyRN*.nc; do echo - seltimestep,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 $m; done; echo -seltimestep,25 $m) OUTFILE.nc)
                1. OP:   NMF - netCDF file with natural mitigation factor due to precipitation

-----------------
|-- Unpaved_Roads
-----------------

        1. PREC_NMF.py - Python script to determine the Natural Mitigation Factor per hour due to precipitation (80% reduction if a precipitation event >0.254mm occurs) NMFp
                1. RQ:   METCRO2D - file with meteorological info. (precipitation information (RN in cm))
                1. OP:   PRECNMF - netCDF file with the natural mitigation factor due to precipitation

        2. REL_HUM.py - Python script to calculate the ambient relative humidity per hour neccesary in the surface moisture calculation - The function in this script is called by AP42_2.py
                2. RQ:   METCRO2D - file with meteorological info. (Temperature, pressure, water vapor mixing ratio)
                2. OP:   RELHUM - file with relative humidity

        3. AP42_2.py - Python script to calculate the surface moisture per hour as suggested by USEPA's spreadsheet using REL_HUM.py script, evaporation, vehicle flux and the moisture that is added from any precipitation event
                3. RQ:   METCRO2D - file with meteorological info.
                         RELHUM - *Optional input - default is calculated online with function from REL_HUM.py; if provided, read from file.
                         SPIN - *Optional input- file to spin-up the results
                3. OP:   SFMOIS - netCDF file with moisture values to calculate NMFm

        4. MOIS.py - Python script to calculate the Natural Mitigation Factor per hour due to surface moisture ((M/0.5)^-0.2)=NMFm
                4. RQ:   SFMOIS - file with the surface moisture (output by AP42_2.py)
                4. OP:   NMF_MOIS - netCDF file with the natural mitigation factor due to moisture

        5. UNPAVED_NMF.py - Python script to calculate the total Natural Mitigation Factor of both NMFp and NMFm NMF=(NMFp*NMFm)
                5. RQ:   PRECNMF - netCDF file with NMFp (output by REC_NMF.py)
                         NMF_MOIS - netCDF file with NMFm (output by MOIS.py)
                5. OP:   NMF - netCDF file with the Natural Mitigation Factor due to precipitation and moisture

---------------
|-- Trans_Fract
---------------

        1. trans_fract.py - Python script to estimate the transportable fraction
                1. RQ:   GRIDCRO2D - File with land use fraction in USGS (otput from MCIP)
                1. OP:   TRANS.nc - netCDF file with the transportable fraction
