"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to handle META eddy trajectory data
 
    1) is_wcrlike                       : returns true if given eddy is warm core ring-like, false otherwise
    2) is_ccrlike                       : returns true if given eddy is cold core ring-like, false otherwise
    3) meta_eddy_to_nwa_ringlike_eddies : saves pandas dataframe of eddy trajectories for wcr-like & ccr-like eddies
    4) get_gs                           : returns overall mean position of the Gulf Stream
    5) get_gs_year                      : takes a given year and returns the average GS path for that year 
    6) get_gs_month                     : takes a given year, month and returns the average GS path for that month
    7) get_gs_day                       : takes a given year, month, day and returns the GS path for that day
    8) get_eddy_formation_loc           : returns lon, lat for a single eddy
    9) get_eddy_formation_time          : returns year, month for a single eddy
    10) get_eddy_lifespan                : returns lifespan of the eddy as integer
    11) eddy_moves_west                  : returns true is eddy moves westward
    12) closest_gs_lat                   : returns closest lon, lat of GS position for a given eddy
    13) is_geq_500m_isobath              : returns true if eddy depth is greater than or equal to 100-m isobath
    14) count_annual_ring_formations     : returns DataFrame of number of annual eddy formations for an eddy DataFrame
    15) count_all_ring_formations        : returns DataFrame of number of monthly, yearly eddy formations for an eddy DataFrame
    16) eddy_df_to_formation_counts_df   : saves DataFrames of merged formation counts for ring-like eddies
    
    
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
# 1) 
def is_wcrlike(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function uses mutliple criterion to filter out warm core rings from mesoscale eddies. The criteria is:
        1) is the eddy anti-cyclonic
        2) is the eddy formed north of -0.25 deg of the GS daily path
        3) is the eddy formed south of +3 deg of the GS daily path
        4) does the eddy propagate westwards
    Essentially, we expect to find WCRs between the GS path and the continental shelf.
    Since the resolution of the datasets is coarse, we use -0.25 deg margin since our first criteria of
    anti-cylonicity would filter out cold core rings from being counted as warm core rings.
    
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)            : returns true if eddy is north of Gulf Stream for given date, else returns false

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    eddy_year = get_eddy_formation_time(eddy)[0] # year of eddy formation
    eddy_month = get_eddy_formation_time(eddy)[1] # month of eddy formation
    eddy_day = get_eddy_formation_time(eddy)[2] # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon

    return ((eddy_lat >= gs_lat_np.T[min_index]-0.25) & (is_geq_100m_isobath(eddy) & (eddy['cyclonic_type']==1).all()) & (eddy_lat <= gs_lat_np.T[min_index]+3) & eddy_moves_west(eddy)) # if true then eddy formation is north of gs path


#-------------------------------------------------------------------------------------------------------------------------------
# 2)
def is_ccrlike(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function uses mutliple criterion to filter out cold core rings from mesoscale eddies. The criteria is:
        1) is the eddy cyclonic
        2) is the eddy formed south of +0.25 deg of the Gulf path
        3) is the eddy formed north of 34N
        4) does the eddy propagate westward
    Essentially, we expect to find CCRs between the GS path and the Sargasso Sea.
    Since the resolution of the datasets is coarse, we use +0.25 deg margin since our first criteria of
    cylonicity would filter out warm core rings from being counted as cold core rings.
    
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)            : returns true if eddy is south of Gulf Stream for given date, else returns false

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    eddy_year = get_eddy_formation_time(eddy)[0] # year of eddy formation
    eddy_month = get_eddy_formation_time(eddy)[1] # month of eddy formation
    eddy_day = get_eddy_formation_time(eddy)[2] # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon
    
    westward_propagation = eddy_moves_west(eddy) # True if eddy moves west
    
    return ((eddy_lat <= gs_lat_np.T[min_index]+0.25) & (is_geq_100m_isobath(eddy)) & (eddy_lat >= gs_lat_np.T[min_index]-3) & ((eddy['cyclonic_type']==-1).all()) & (eddy_lat >= 34) & westward_propagation)


#-------------------------------------------------------------------------------------------------------------------------------
# 3)
def meta_eddy_to_nwa_ringlike_eddies(path):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    " Function Description:                                                                                   "
    "    The meta_eddy_to_nwa_ringlike_eddies function reads in META eddy datasets and converts it to         " 
    "    a pandas dataframe, filters out eddies that don't qualify as WCR-like or CCR-like, and saves         "
    "    the WCR-like and CCR-like eddy dataframes in the data/dataframes folder.                             "
    "                                                                                                         "
    " Input:                                                                                                  " 
    "    path (String)         : path to where the META eddy dataset is stored                                "
    "                                                                                                         "
    " Output:                                                                                                 "
    "    (None)                : saves pandas dataframe of eddy trajectories for wcr-like and ccr-like eddies "
    "                                                                                                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                 read in META eddy trajectory datasets                                   "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # read in META2.0 .nc file as xarray dataset
    meta2_ds = xr.open_dataset(path)

    # convert xr dataset to pandas dataframe
    meta2_df = meta2_ds.to_dataframe()
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                        restrict dataset to nwa box & change lon to -180 to 180                          "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # change longitude from 0 to 360 --> -180 to 180
    meta2_df['longitude'] = ((((meta2_df['longitude'] + 180) % 360) - 180))

    # bounding box for NWA based on Gangopadhyay et al., 2019
    zone_lat = [30,45] # N
    zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]] # W ... 0 –> all zones, 1 –> zone 1, etc.

    # cut eddy trajectory dataframe down to NWA region
    meta2_df = meta2_df[(meta2_df['longitude'] >= zone_lon[0][0]) & (meta2_df['longitude'] <= zone_lon[0][1]) & (meta2_df['latitude'] >= zone_lat[0]) & (meta2_df['latitude'] <= zone_lat[1])]

    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                 filter for wcr-like & ccr-like eddies                                   "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # separate anticyclonic and cyclonic eddies into their own dataframes
    anticyclonic_eddy_df = meta2_df[meta2_df['cyclonic_type']==1]
    cyclonic_eddy_df = meta2_df[meta2_df['cyclonic_type']==-1]


    """"""""""""""""""""""""""""""""""""""" wcr-like eddies """""""""""""""""""""""""""""""""""""""""""""""""""
    wcrlike_tracks = []
    # loop through each eddy to determine if it is wcr-like
    for i in np.array(anticyclonic_eddy_df['track'].unique()):
        eddy = anticyclonic_eddy_df[anticyclonic_eddy_df['track']==i]
        if ((eddy['time'].dt.year<2018).all()): # stops at 2017 because that is common cut-off year
            if is_wcrlike(eddy).all():
                wcrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_wcrlike_df = anticyclonic_eddy_df[anticyclonic_eddy_df['track'].isin(wcrlike_tracks)]


    # """"""""""""""""""""""""""""""""""""""" ccr-like eddies """""""""""""""""""""""""""""""""""""""""""""""""""
    ccrlike_tracks = []
    # loop through each eddy to determine if it is ccr-like
    for i in np.array(cyclonic_eddy_df['track'].unique()):
        eddy = cyclonic_eddy_df[cyclonic_eddy_df['track']==i]
        if ((eddy['time'].dt.year<2018).all()): # stops at 2017 because that is common cut-off year
            if is_ccrlike(eddy).all():
                ccrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_ccrlike_df = cyclonic_eddy_df[cyclonic_eddy_df['track'].isin(ccrlike_tracks)]
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                       save filtered dataframes for wcr- & ccr-like eddies                               "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # save wcr-like df as pickled file
    eddy_wcrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/nwa_wcrlike_eddies.pkl')
    
    # save ccr-like df as pickled file
    eddy_ccrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/nwa_ccrlike_eddies.pkl') 


#-------------------------------------------------------------------------------------------------------------------------------
# 4)




#-------------------------------------------------------------------------------------------------------------------------------
# 5)




#-------------------------------------------------------------------------------------------------------------------------------
# 6)




#-------------------------------------------------------------------------------------------------------------------------------
# 7)



#-------------------------------------------------------------------------------------------------------------------------------
# 8)




#-------------------------------------------------------------------------------------------------------------------------------
# 9)




#-------------------------------------------------------------------------------------------------------------------------------
# 10)




#-------------------------------------------------------------------------------------------------------------------------------
# 11)




#-------------------------------------------------------------------------------------------------------------------------------
# 12)




#-------------------------------------------------------------------------------------------------------------------------------
# 13)




#-------------------------------------------------------------------------------------------------------------------------------
# 14)




#-------------------------------------------------------------------------------------------------------------------------------
# 15)




#-------------------------------------------------------------------------------------------------------------------------------
# 16)




#-------------------------------------------------------------------------------------------------------------------------------
# 17)




#-------------------------------------------------------------------------------------------------------------------------------


