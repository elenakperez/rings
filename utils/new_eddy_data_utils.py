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

def meta_eddy_to_nwa_ringlike_eddies():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                 read in META eddy trajectory datasets                                   "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # read in META2.0 .nc file as xarray dataset
    meta2_ds = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/META/META2_0/eddy_trajectory_dt_2.0_19930101_20200307.nc')

    # convert xr dataset to pandas dataframe
    meta2_df = meta2_ds.to_dataframe()
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                        restrict dataset to nwa box & change lon to -180 to 180                          "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # bounding box for NWA based on Gangopadhyay et al., 2019
    zone_lat = [30,45] # N
    zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]] # W ... 0 –> all zones, 1 –> zone 1, etc.
    
    # cut eddy trajectory dataframe down to NWA region
    meta2_df = meta2_df[(meta2_df['longitude'] >= zone_lon[0][0]) & (meta2_df['longitude'] <= zone_lon[0][1]) & (meta2_df['latitude'] >= zone_lat[0]) & (meta2_df['latitude'] <= zone_lat[1])]

    # change longitude from 0 to 360 --> -180 to 180
    meta2_df['longitude'] = ((((meta2_df['longitude'] + 180) % 360) - 180))
     
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
        anticyclonic_eddy = anticyclonic_eddy_df[anticyclonic_eddy_df['track']==i]
        if ((anticyclonic_eddy['time'].dt.year<2018).all()): # stops at 2017 because that is common cut-off year
            if is_wcrlike(anticyclonic_eddy).all():
                wcrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_wcrlike_df = anticyclonic_eddy[anticyclonic_eddy['track'].isin(wcrlike_tracks)]
    
    # save wcr-like df as pickled file
    eddy_wcrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/nwa_wcrlike_eddies.pkl') 
    
    
    """"""""""""""""""""""""""""""""""""""" ccr-like eddies """""""""""""""""""""""""""""""""""""""""""""""""""
    ccrlike_tracks = []
    # loop through each eddy to determine if it is ccr-like
    for i in np.array(cyclonic_eddy_df['track'].unique()):
        cyclonic_eddy = anticyclonic_eddy_df[cyclonic_eddy_df['track']==i]
        if ((cyclonic_eddy['time'].dt.year<2018).all()): # stops at 2017 because that is common cut-off year
            if is_ccrlike(cyclonic_eddy).all():
                ccrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_ccrlike_df = cyclonic_eddy[cyclonic_eddy['track'].isin(ccrlike_tracks)]

    # save ccr-like df as pickled file
    eddy_ccrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/nwa_ccrlike_eddies.pkl') 
    

# meta_eddy_to_nwa_eddylike_rings()

#-------------------------------------------------------------------------------------------------------------------------------


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

    return ((eddy_lat >= gs_lat_np.T[min_index]-0.25) & (is_geq_500m_isobath(eddy) & (eddy['cyclonic_type']==1).all()) & (eddy_lat <= gs_lat_np.T[min_index]+3) & eddy_moves_west(eddy)) # if true then eddy formation is north of gs path

#-------------------------------------------------------------------------------------------------------------------------------
# ) returns true if given eddy is cold core ring, false otherwise

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
    
    return ((eddy_lat <= gs_lat_np.T[min_index]+0.25) & (is_geq_500m_isobath(eddy)) & (eddy_lat >= gs_lat_np.T[min_index]-3) & ((eddy['cyclonic_type']==-1).all()) & (eddy_lat >= 34) & westward_propagation)


#-------------------------------------------------------------------------------------------------------------------------------
# ) returns overall mean position of the Gulf Stream 1993 - 2018

def get_gs():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
    
    Output:
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/gulf_stream/GS.mat')
    
    # save Gulf Stream data as xarray
    return (gs['lon_cell'][0][0][0] - 360),(gs['lat_cell'][0][0][0])

#-------------------------------------------------------------------------------------------------------------------------------
# ) return the overall position of the GS for that given year

def get_gs_year(year):
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input: 
        year (Int)          : the year for which you want to extract GS position
    
    Output: 
        gs_lon_year (array) : 1D array of GS longitude 
        gs_lat_year (array) : 1D array of GS latitude 
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/gulf_stream/GS.mat')
    
    return (gs['lon_cell_yearly'][year-1993][0][0] - 360),(gs['lat_cell_yearly'][year-1993][0][0])

#-------------------------------------------------------------------------------------------------------------------------------
# ) returns overall position of GS for that month of the given year

def get_gs_month(year,month):
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        year (Int)           : the year for which you want to extract GS position
        month (Int)          : the month for which you want to extract GS position
    
    Output:
        gs_lon_month (array) : 1D array of GS longitude for that month
        gs_lat_month (array) : 1D array of GS latitude for that month
        
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/gulf_stream/GS.mat')
    
    return (gs['lon_cell_monthly'][(((year-1993)*12)+month)][0][0] - 360),(gs['lat_cell_monthly'][(((year-1993)*12)+month)][0][0])


#-------------------------------------------------------------------------------------------------------------------------------

# ) create get_gs_day to return gs position for given day 1-Jan-1993 to 5-Oct-2022'

def get_gs_day(year,month,day):
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        year (Int)           : the year for which you want to extract GS position (1993-2022)
        month (Int)          : the month for which you want to extract GS position (Jan 1993 to Oct 2022)
        day (Int)            : the day for which you want to extract GS position (Jan 1 1993 to Oct 5 2022)
    
    Output:
        gs_lon_month (array) : 1D array of GS longitude for that month
        gs_lat_month (array) : 1D array of GS latitude for that month
        
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/gulf_stream/GS_daily_CMEMS_047_50cm_contours_1993_to_nrt.mat')
    
    'convert time array to ordinal dates'
    for d in range(len(gs['time'][0])-1):
        gs['time'][0][d] = gs['time'][0][d]+date.toordinal(date(1950,1,1))
    
    'get index in lon/lat arrays from time array'
    index_from_date = np.where(gs['time'][0]==date.toordinal(date(year,month,day)))[0][0]
    
    return (gs['lon_cell'][index_from_date][0][0] - 360),(gs['lat_cell'][index_from_date][0][0])


#-------------------------------------------------------------------------------------------------------------------------------

# ) returns lon, lat for a single eddy

def get_eddy_formation_loc(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns lon/lat location of eddy formation as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    return eddy[eddy['observation_number']==0]['longitude'],eddy[eddy['observation_number']==0]['latitude']

#-------------------------------------------------------------------------------------------------------------------------------

# ) returns location (lon, lat) of demise for a single eddy

def get_eddy_demise_loc(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns lon/lat location of eddy demise as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    return eddy.iloc[-1].longitude,eddy.iloc[-1].latitude

#-------------------------------------------------------------------------------------------------------------------------------

# ) returns year, month for a single eddy

def get_eddy_formation_time(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns year, month, day of eddy formation as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # old eddy['time'].iloc[0].year, eddy['time'].iloc[0].month 
    return eddy[eddy['observation_number']==0]['time'].dt.year,eddy[eddy['observation_number']==0]['time'].dt.month,eddy[eddy['observation_number']==0]['time'].dt.day

#-------------------------------------------------------------------------------------------------------------------------------

# ) returns lifespan of eddy

def get_eddy_lifespan(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns integer equal to lifespan of eddy

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    return eddy['observation_number'].iloc[-1]

#-------------------------------------------------------------------------------------------------------------------------------