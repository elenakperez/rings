"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to process CMEMS geostrophic velocity data
 
    1) create_masks          : creates and saves masks for the Northwest Atlantic, pre-processes daily vels data to monthly
    2) geo_vels_anom         : creates and saves DataArray of geostrophic velocity anomalies, EKE, KE, and speed 
    3) calc_area_regions     : calculates time-varying regional area of Slope, Gulf Stream, Sargasso and saves as dataArrays 

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# import functions for data analysis 
import requests
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import xesmf as xe
import numpy as np
import pandas as pd
import netCDF4 as nc
from scipy.io import loadmat
from scipy.interpolate import griddata
import pickle
import datetime
from datetime import date
from scipy.stats.stats import pearsonr   

# import the util functions
from utils.ring_data_utils import * # for getting Gulf Stream paths

# turn off warnings
import warnings
warnings.filterwarnings("ignore")


#-------------------------------------------------------------------------------------------------------------------------------
# 1) 
def create_masks():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        This function create and saves masks for the Northwest Atlantic based on lat/lon boundaries
        (e.g. Zone 1: 75-70W, Zone 2: 70-65W, Zone 3: 65-60W, Zone 4: 60-55W). As well as,
        regional masks based on the time-varying monthly mean Gulf Stream boundary. 
        Additionally, this function pre-processes the geostrophic velocity file so that daily data is
        resampled to monthly.

        Input:
            * None, it loads necessary files when called 

        Output:
            * None, it saves all masks in respective folders (folder for xarray DataArrays, np Arrays, etc.)
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    ####################################### OPEN DATA ##########################################
    # geostrophic velocity file from CMEMS, daily data
    geo_vels = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/CMEMS/cmems_geos.nc').sel(time = slice(None, "2017"))
    # resample geostrophic velocity file from daily data to monthly data
    geo_vels = geo_vels.resample(time = "1M").mean()
    
    bathy = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/nwa_bathy.nc')
    # regrid bathy file into grid like the geostrophic velocity file
    bathy = bathy.z.interp(lat=geo_vels.latitude).interp(lon=geo_vels.longitude)
    ############################################################################################    
    
    
    ################################## BATHY MASK #############################################
    # create a mask for bathy that only includes values deeper than the 300-m isobath
    mask_bathy = (xr.where(bathy<=-300, 1, np.nan))

    # create a meshgrid from the regridded bathy file 
    LON, LAT = np.meshgrid(bathy.lon, bathy.lat)

    
    ################################### LON MASKS #############################################
    # create a mask that includes all values east of 75W
    mask_75 = 1.0 *(LON>=-75)
    mask_75[np.where(mask_75==0)] = np.nan
    
    # create a mask that includes all values east of 75W
    mask_75_west = 1.0 *(LON>=-75)
    mask_75_west[np.where(mask_75_west==0)] = np.nan
    
    # create masks for values east/west of 70W
    mask_70_east = 1.0 * (LON<=-70)
    mask_70_east[np.where(mask_70_east==0)] = np.nan
    mask_70_west = 1.0 * (LON>=-70)
    mask_70_west[np.where(mask_70_west==0)] = np.nan
    
    # create masks for values east/west of 65W
    mask_65_east = 1.0 * (LON<=-65)
    mask_65_east[np.where(mask_65_east==0)] = np.nan
    mask_65_west = 1.0 * (LON>=-65)
    mask_65_west[np.where(mask_65_west==0)] = np.nan

    # create masks for values east/west of 60W
    mask_60_east = 1.0 * (LON<=-60)
    mask_60_east[np.where(mask_60_east==0)] = np.nan
    mask_60_west = 1.0 * (LON>=-60)
    mask_60_west[np.where(mask_60_west==0)] = np.nan

    # create a mask that includes all values west of 55W
    mask_55_east = 1.0 * (LON<=-55.25)
    mask_55_east[np.where(mask_55_east==0)] = np.nan
    
    
    ################################## ZONE MASKS #############################################
    # create masks for the Zones
    mask_zone1 = mask_75_west*mask_70_east
    mask_zone2 = mask_70_west*mask_65_east
    mask_zone3 = mask_65_west*mask_60_east
    mask_zone4 = mask_60_west*mask_55_east

    ################################## NWA MASK ###############################################
    # create a mask for the NWA region (30-45N, )
    mask_NWA = mask_75_west * mask_55_east * mask_bathy
    
    
    ################################ SAVE LAT/LON MASKS ########################################
    mask_bathy.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/masks/mask_bathy.nc')
    
    np.save('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_75_west.npy', mask_75_west) 
    
    np.save('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_55_east.npy', mask_55_east)
    
    np.save('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone1.npy', mask_zone1)
    np.save('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone2.npy', mask_zone2)
    np.save('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone3.npy', mask_zone3)
    np.save('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone4.npy', mask_zone4)
    
    mask_NWA.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/masks/mask_NWA.nc') 
    ############################################################################################
    
    
    ################################## REGIONAL MASKS ##########################################
    ## create mask for SLOPE
    # create empty grid for Slope mask and store in geostrophic velocity file
    geo_vels["mask_slope"] = geo_vels.ugos * np.nan
    geo_vels["mask_sag"] = geo_vels.ugos * np.nan

    # for loop to create mask for all months, years 1993–2017
    for i in range(len(geo_vels.time)):
            # pull month, year for current time step
            month = geo_vels.time.dt.month[i]
            year = geo_vels.time.dt.year[i]

            # pull Gulf Stream monthly mean position for current time step
            LON_GS, LAT_GS = get_gs_month(year,month)[0], get_gs_month(year,month)[1]


            ################################## SLOPE SEA #######################################
            # create temporary Slope Sea mask dataset
            ds_locs = xr.Dataset()
            ds_locs['lon'] = xr.DataArray(data=LON_GS, dims=('location'))
            ds_locs['lat'] = xr.DataArray(data=LAT_GS, dims=('location'))
            ds_locs['mask_slope'] = xr.DataArray(data= np.ones(len(LON_GS)), dims=('location'))

            # regrid temporary Slope mask to Gulf Stream resolution
            regrid_2_ds = xe.Regridder(ds_locs, bathy, "nearest_d2s", locstream_in = True)
            slope_mask = xr.where(regrid_2_ds(ds_locs)['mask_slope'] != 0.0, np.nan, 1.0)

            # create second temporary mask that will be used to subtract everything below the GS path from the mask
            slope_mask2 = 1 * slope_mask
            slope_mask2.values = np.cumsum(slope_mask2.values[::-1, :], 0)[::-1, :]
            slope_mask.values = np.cumsum(slope_mask.values, 0)

            # set empty values in the mask to NaNs
            slope_mask = xr.where(np.isnan(slope_mask), 1.0, np.nan)
            slope_mask2 = xr.where(np.isnan(slope_mask2), np.nan, 1.0)

            # merge the masks to create Slope Sea mask (75-55W, 300-m isobath–Gulf Stream monthly north wall path)
            slope_mask = slope_mask * slope_mask2 * mask_bathy * mask_55_east * mask_75_west

            # store the Slope Sea mask in the geostrophic velocity file 
            geo_vels["mask_slope"].values[i, :, :] =  slope_mask.values


            ################################ SARGASSO SEA ######################################
            # create temporary Sargasso Sea mask dataset
            ds_locs = xr.Dataset()
            ds_locs['lon'] = xr.DataArray(data=LON_GS, dims=('location'))
            ds_locs['lat'] = xr.DataArray(data=LAT_GS, dims=('location'))
            ds_locs['mask_sag'] = xr.DataArray(data= np.ones(len(LON_GS)), dims=('location'))

            # regrid temporary Sargasso Sea mask to Gulf Stream resolution
            regrid_2_ds = xe.Regridder(ds_locs, bathy, "nearest_d2s", locstream_in = True, periodic = True)
            sag_mask = xr.where(regrid_2_ds(ds_locs)['mask_sag'] != 0.0, np.nan, 1.0)
            sag_mask = xr.where(np.isnan(sag_mask), np.nan, 1.0)

            # select offset for Gulf Stream southern boundary
            offset = 2.0; 
            lat_subs = LAT_GS[(LON_GS > -75) * (LON_GS < -55.25)]; 
            lons_subs = LON_GS[(LON_GS > -75) * (LON_GS < -55.25)]; 

            # create a secondary temporary Sargasso Sea mask dataset
            ds_locs = xr.Dataset()
            ds_locs['lon'] = xr.DataArray(data=LON_GS, dims=('location'))
            ds_locs['lat'] = xr.DataArray(data=LAT_GS-2, dims=('location'))
            ds_locs['mask_sag'] = xr.DataArray(data= np.ones(len(LON_GS)), dims=('location'))
            regrid_2_ds = xe.Regridder(ds_locs, bathy, "nearest_d2s", locstream_in = True, periodic = True)
            sag_mask2 = xr.where(regrid_2_ds(ds_locs)['mask_sag'] != 0.0, np.nan, 1.0)

            # create the northern Gulf Stream boundary mask
            sag_mask = mask_55_east * mask_75_west * sag_mask2
            sag_mask = xr.where(np.isnan(sag_mask), np.nan, 1.0)

            # create the southern Gulf Stream boundary mask
            sag_mask2.values = np.cumsum(sag_mask2.values[::-1, :], 0)[::-1, :]
            sag_mask.values = np.cumsum(sag_mask.values, 0)

            sag_mask = xr.where(np.isnan(sag_mask), np.nan, 1.0)

            geo_vels["mask_sag"].values[i, :, :] =  sag_mask.values

            
    ################################### GULF STREAM ############################################
    # combine the northern and southern boundary masks to create the Gulf Stream mask
    gs_south_mask =  xr.where(np.isnan(geo_vels["mask_sag"]), 1.0, np.nan)
    gs_north_mask =  xr.where(np.isnan(geo_vels["mask_slope"]), 1.0, np.nan)
    gs_mask = gs_south_mask + gs_north_mask

    # create the final Gulf Stream mask using the Slope Sea mask and Sargasso Sea mask as the north/south boundaries
    geo_vels["mask_gs"] = gs_mask * mask_55_east * mask_75_west * mask_bathy

    # fill empty values with NaNs and mask values with 1s
    # store the Gulf Stream mask in the geostrophic velocity file
    geo_vels["mask_gs"] =  xr.where(np.isnan(geo_vels["mask_gs"]), np.nan, 1.0)
    
    
    
    ###################################### SAVE MASKs ##########################################
    geo_vels.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/geo_vels_proc.nc')
    # save the new geo_vels as a DataArray
    ############################################################################################

#-------------------------------------------------------------------------------------------------------------------------------
# 2)
def geo_vels_anom():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        This function creates and saves an xarray DataArray of geostrophic velocity anomalies and 
        corresponding 

        Input:
            * None, it loads necessary files when called 

        Output:
            * None, it saves all masks in respective folders (folder for xarray DataArrays, np Arrays, etc.)
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""    
    ####################################### OPEN DATA ##########################################
    # geostrophic velocity pre-processed
    geo_vels = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/geo_vels_proc.nc')

    bathy = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/nwa_bathy.nc')
    # regrid bathy file into grid like the geostrophic velocity file
    bathy = bathy.z.interp(lat=geo_vels.latitude).interp(lon=geo_vels.longitude)

    # open lat/lon masks
    mask_bathy = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/masks/mask_bathy.nc')

    mask_75_west = np.load('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_75_west.npy')
    mask_55_east = np.load('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_55_east.npy')

    mask_zone1 = np.load('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone1.npy')
    mask_zone2 = np.load('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone2.npy')
    mask_zone3 = np.load('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone3.npy')
    mask_zone4 = np.load('/Users/elenaperez/Desktop/rings/data/np_arrays/masks/mask_zone4.npy')

    mask_NWA = xr.open_dataarray('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/masks/mask_NWA.nc') 
    ############################################################################################


    ####################### COMPUTE GEOSTROPHIC VELOCITY ANOMALIES #############################
    # create a new xarray Dataset for computation of EKE, KE, and speed
    # geostrophic velocity anomalies = monthly geostrophic velocity - climatology of geostrophic velocity (1993–2017)
    ga = geo_vels.groupby("time.month") - geo_vels.groupby("time.month").mean("time")
    ga = ga[["ugos", "vgos"]]

    ####################################### COMPUTE EKE ########################################
    # compute EKE of the whole region (85–55W, 30-45N)
    ga["EKE"] = ((ga.ugos**2) + (ga.vgos**2))
    ga["EKE"] = 1025 * ga["EKE"] / 2 

    # mask regions 
    ga['EKE_NWA'] = ga["EKE"] * mask_NWA
    ga['EKE_slope'] = geo_vels["mask_slope"] * ga["EKE"]
    ga['EKE_gs'] = geo_vels["mask_gs"] * ga["EKE"]
    ga['EKE_sag'] = geo_vels["mask_sag"] * ga["EKE"]
    ga['EKE_zone1'] = ga["EKE"] * mask_zone1
    ga['EKE_zone2'] = ga["EKE"] * mask_zone2
    ga['EKE_zone3'] = ga["EKE"] * mask_zone3
    ga['EKE_zone4'] = ga["EKE"] * mask_zone4



    ####################################### COMPUTE KE #########################################
    # compute KE of the whole region 
    ga['KE'] = (((geo_vels.ugos)**2) + ((geo_vels.vgos)**2))
    ga['KE'] = ((1025 * ga['KE'])/2)

    # mask regions 
    ga['KE_NWA'] = ga["KE"] * mask_NWA
    ga['KE_slope'] = geo_vels["mask_slope"] * ga["KE"]
    ga['KE_gs'] = geo_vels["mask_gs"] * ga["KE"]
    ga['KE_sag'] = geo_vels["mask_sag"] * ga["KE"]
    ga['KE_zone1'] = ga["KE"] * mask_zone1
    ga['KE_zone2'] = ga["KE"] * mask_zone2
    ga['KE_zone3'] = ga["KE"] * mask_zone3
    ga['KE_zone4'] = ga["KE"] * mask_zone4


    ################################### COMPUTE SPEED #########################################
    # compute speed of the whole region
    ga['speed'] = np.sqrt(((geo_vels.ugos)**2) + ((geo_vels.vgos)**2))

    # mask regions 
    ga['speed_NWA'] = ga["speed"] * mask_NWA
    ga['speed_slope'] = geo_vels["mask_slope"] * ga["speed"]
    ga['speed_gs'] = geo_vels["mask_gs"] * ga["speed"]
    ga['speed_sag'] = geo_vels["mask_sag"] * ga["speed"]
    ga['speed_zone1'] = ga["speed"] * mask_zone1
    ga['speed_zone2'] = ga["speed"] * mask_zone2
    ga['speed_zone3'] = ga["speed"] * mask_zone3
    ga['speed_zone4'] = ga["speed"] * mask_zone4


    ###################################### SAVE MASKs ##########################################
    ga.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/geo_vels_anom.nc')
    # save the new geostrophic velocity anomaly DataArray 
    ############################################################################################
    
    
#-------------------------------------------------------------------------------------------------------------------------------
# 3)
def calc_area_regions():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    * None, it loads all necessary files to compute area of the regions                                "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * None, saves the regional, time-varying areas as DataArrays in the xr_dataarrays folder           "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    ####################################### OPEN DATA ##########################################
    geo_vels = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/geo_vels_proc.nc')
    
    
    ######################## CREATE GRID OF SPHERICAL EARTH ############################
    boxlo,boxla=np.array(np.meshgrid(geo_vels.longitude,geo_vels.latitude))
    grid=np.cos(np.radians(abs(boxla)))*(111.1*111.1*0.25*0.25)
    
    
    ################################## SLOPE SEA #######################################
    regionMask = 'mask_slope'
    slopeArea_monthly = ((grid*geo_vels[regionMask]).sum(("longitude", "latitude"))  / 1000)
    slopeArea_annual = ((grid*geo_vels[regionMask]).sum(("longitude", "latitude"))  / 1000).resample(time='1Y')
    
    
    ################################ GULF STREAM ########################################
    regionMask = 'mask_gs'
    gsArea_monthly = ((grid*geo_vels[regionMask]).sum(("longitude", "latitude"))  / 1000)
    gsArea_annual = ((grid*geo_vels[regionMask]).sum(("longitude", "latitude"))  / 1000).resample(time='1Y')

    
    ################################ SARGASSO SEA ######################################
    regionMask = 'mask_sag'
    sagArea_monthly = ((grid*geo_vels[regionMask]).sum(("longitude", "latitude"))  / 1000)
    sagArea_annual = ((grid*geo_vels[regionMask]).sum(("longitude", "latitude"))  / 1000).resample(time='1Y')
    
    
    ############################## SAVE DATAARRAYS ####################################
    slopeArea_monthly.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/slopeArea_monthly.nc')
    gsArea_monthly.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/gsArea_monthly.nc')
    sagArea_monthly.to_netcdf('/Users/elenaperez/Desktop/rings/data/xr_dataarrays/sagArea_monthly.nc')


#-------------------------------------------------------------------------------------------------------------------------------
# 4) example function
# def func_name(variable1, variable2, ...):
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
#     This function __[describe what the function does]___________________
    
#     Input:
#         variable1 (type)           : 
#         variable2 (type)           : 

#     Output:
#         variable_out (type)        : returns 

#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#-------------------------------------------------------------------------------------------------------------------------------
