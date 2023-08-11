"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to manipulate Chelton Eddy Tracks ()
 
    ) eddy_df_old_version            : !OLD! create a pandas dataframe for eddies in the nwa
    ) make_eddy_df_old_version       : !OLD! create zonal eddy dataframe from nwa tracks
    ) create_nwa_tracks              : saves nwa_tracks as an np array of track ids for nwa eddies
    ) get_gs                         : returns overall mean position of the Gulf Stream
    ) get_gs_year                    : takes a given year and returns the average GS path for that year 
    ) get_gs_month                   : takes a given year, month and returns the average GS path for that month
    ) get_gs_day                     : takes a given year, month, day and returns the GS path for that day
    ) get_eddy_formation_loc         : returns lon, lat for a single eddy
    ) get_eddy_formation_time        : returns year, month for a single eddy
    ) get_eddy_lifespan              : returns lifespan of the eddy as integer
    ) is_wcr                         : returns true if given eddy is warm core ring, false otherwise
    ) is_ccr                         : returns true if given eddy is cold core ring, false otherwise
    ) eddy_moves_west                : returns true is eddy moves westward
    ) closest_gs_lat                 : returns closest lon, lat of GS position for a given eddy
    ) is_geq_500m_isobath            : returns true if eddy depth is greater than or equal to 100-m isobath
    ) make_eddy_region_df            : takes Chelton tracks and saves a processed eddy DataFrame as pickled file
    ) make_eddy_zone_df              : cuts down eddy DataFrame to smaller zone + saves new DataFrame as pickle
    ) eddy_df_to_ring_df             : takes an eddy DataFrame and saves a new DataFrame of just WCRs or CCRs
    ) count_monthly_ring_formations  : saves DataFrame of number of monthly eddy formations for an eddy DataFrame
    ) count_annual_ring_formations   : saves DataFrame of number of annual eddy formations for an eddy DataFrame
    ) count_all_ring_formations      : saves DataFrame of number of monthly, yearly eddy formations for an eddy DataFrame
    ) merge_ring_monthly_counts      : saves merged DataFrame of monthly ring formations by zone
    ) merge_ring_annual_counts       : saves merged DataFrame of annual ring formations by zone
    ) merge_ring_allf_counts          : saves merged DataFrame of monthly, yearly ring formations by zone
    ) count_monthly_eddy_formations  : saves DataFrame of number of monthly eddy formations for an eddy DataFrame
    ) count_annual_eddy_formations   : saves DataFrame of number of annual eddy formations for an eddy DataFrame
    ) count_all_eddy_formations      : saves DataFrame of number of monthly, yearly eddy formations for an eddy DataFrame
    ) merge_eddy_monthly_counts      : saves merged DataFrame of monthly ring formations by zone
    ) merge_eddy_annual_counts       : saves merged DataFrame of annual ring formations by zone
    ) merge_eddy_all_counts          : saves merged DataFrame of monthly, yearly ring formations by zone

    
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# import necessary packages and functions
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
from scipy.io import loadmat
from scipy.interpolate import griddata
import pickle
import datetime
from datetime import date

# turn off warnings
import warnings
warnings.filterwarnings("ignore")


#-------------------------------------------------------------------------------------------------------------------------------

