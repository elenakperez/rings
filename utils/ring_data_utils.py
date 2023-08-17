"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to handle META eddy trajectory data
 
    1) is_wcrlike                        : returns true if given eddy is warm core ring-like, false otherwise
    2) is_ccrlike                        : returns true if given eddy is cold core ring-like, false otherwise
    3) meta_eddy_to_nwa_ringlike_eddies  : saves pandas dataframe of eddy trajectories for wcr-like & ccr-like eddies
    4) get_gs                            : returns overall mean position of the Gulf Stream
    5) get_gs_year                       : takes a given year and returns the average GS path for that year 
    6) get_gs_month                      : takes a given year, month and returns the average GS path for that month
    7) get_gs_day                        : takes a given year, month, day and returns the GS path for that day
    8) get_eddy_formation_loc            : returns lon, lat for a single eddy
    9) get_eddy_formation_time           : returns year, month for a single eddy
   10) get_eddy_lifespan                 : returns lifespan of the eddy as integer
   11) eddy_moves_west                   : returns true is eddy moves westward
   12) eddy_moves_east                   : returns true is eddy moves eastward
   13) closest_gs_lat                    : returns closest lon, lat of GS position for a given eddy
   14) is_geq_500m_isobath               : returns true if eddy depth is greater than or equal to 100-m isobath
   15) count_annual_ring_formations      : returns DataFrame of number of annual eddy formations for an eddy DataFrame
   16) count_all_ring_formations         : returns DataFrame of number of monthly, yearly eddy formations for an eddy DataFrame
   17) eddy_df_to_formation_counts_df    : saves DataFrames of merged formation counts for ring-like eddies
    
    
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
def get_gs():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
    
    Output:
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/GS.mat')
    
    # save Gulf Stream data as xarray
    return (gs['lon_cell'][0][0][0] - 360),(gs['lat_cell'][0][0][0])


#-------------------------------------------------------------------------------------------------------------------------------
# 5) return the overall position of the GS for that given year
def get_gs_year(year):
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input: 
        year (Int)          : the year for which you want to extract GS position
    
    Output: 
        gs_lon_year (array) : 1D array of GS longitude 
        gs_lat_year (array) : 1D array of GS latitude 
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/GS.mat')
    
    return (gs['lon_cell_yearly'][year-1993][0][0] - 360),(gs['lat_cell_yearly'][year-1993][0][0])


#-------------------------------------------------------------------------------------------------------------------------------
# 6) returns overall position of GS for that month of the given year
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
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/GS.mat')
    
    return (gs['lon_cell_monthly'][(((year-1993)*12)+month)][0][0] - 360),(gs['lat_cell_monthly'][(((year-1993)*12)+month)][0][0])


#-------------------------------------------------------------------------------------------------------------------------------
# 7)
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
# 8)
def get_eddy_formation_time(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns year, month, day of eddy formation as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # v0: eddy['time'].iloc[0].year, eddy['time'].iloc[0].month 
    # v1: eddy[eddy['observation_number']==0]['time'].dt.year,eddy[eddy['observation_number']==0]['time'].dt.month,eddy[eddy['observation_number']==0]['time'].dt.day
    return np.array(eddy['time'].dt.year)[0], np.array(eddy['time'].dt.month)[0], np.array(eddy['time'].dt.day)[0]


#-------------------------------------------------------------------------------------------------------------------------------
# 9)
def get_eddy_formation_loc(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns lon/lat location of eddy formation as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # v0: eddy[eddy['observation_number']==0]['longitude'],eddy[eddy['observation_number']==0]['latitude']
    return np.array(eddy['longitude'])[0], np.array(eddy['latitude'])[0]


#-------------------------------------------------------------------------------------------------------------------------------
# 10)
def get_eddy_lifespan(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns integer equal to lifespan of eddy

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    return eddy['observation_number'].iloc[-1]

#-------------------------------------------------------------------------------------------------------------------------------
# 11)
def eddy_moves_west(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function determines if an eddy propagates net westward or eastward by comparing the starting 
    position to the ending position. If start is east of the end position (i.e. eddy cumulatively moved west), 
    then the function returns True.
    
    Input:
        eddy (DataFrame) : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)           : returns true if eddy propagates westwards, false otherwise

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    eddy_lon_start = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lon_end = get_eddy_demise_loc(eddy)[0] # lon of eddy demise
    
    return (eddy_lon_start > eddy_lon_end)


#-------------------------------------------------------------------------------------------------------------------------------
# 12)
def eddy_moves_east(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function determines if an eddy propagates net eastward by comparing the starting 
    position to the ending position. If start is west of the end position (i.e. eddy cumulatively moved east), 
    then the function returns True.
    
    Input:
        eddy (DataFrame) : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)           : returns true if eddy propagates eastward, false otherwise

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    eddy_lon_start = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lon_end = get_eddy_demise_loc(eddy)[0] # lon of eddy demise
    
    return (eddy_lon_start < eddy_lon_end)


#-------------------------------------------------------------------------------------------------------------------------------
# 13)
def is_geq_100m_isobath(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function determins if an eddy is *not* on the shelf. If its depth at bird is greater than or equal to
    the 100-m isobath, return True since it is *not* on the shelf. If its depth is less than 100-m,
    return False because it is on the shelf.
    
    Input:
        eddy (DataFrame) : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)           : returns true if eddy is deeper than 500m isobath (on shelf), returns false otherwise

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                         Bathymetry for Northwest Atlantic [24,53]N & [-82,-48]W                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    bathy = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/nwa_bathy.nc')
    
    target_bathy = 100 # desired isobath

    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    bathy_lon_np = bathy.lon.values
    bathy_lon_len = len(bathy_lon_np) # length of bathy lon array

    eddy_lon_np = np.full(bathy_lon_len, eddy_lon)

    min_lon_index = np.argmin(abs(bathy_lon_np-eddy_lon_np)) # index of closest lon in bathy file

    bathy_lat_np = bathy.z.sel(lon=bathy_lon_np[min_lon_index]).lat.values # array of bathy lats at the eddy lon
    bathy_lat_len = len(bathy_lat_np)

    eddy_lat_np = np.full(bathy_lat_len, eddy_lat) 

    min_lat_index = np.argmin(abs(bathy_lat_np-eddy_lat_np)) # index of closest lat in bathy file

    bathy_depth = bathy.z.sel(lat=bathy_lat_np[min_lat_index],lon=bathy_lon_np[min_lon_index]).values # depth of 

    return (target_bathy <= abs(bathy_depth))




#-------------------------------------------------------------------------------------------------------------------------------
# 14)
def is_geq_100m_isobath(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function determins if an eddy is *not* on the shelf. If its depth at bird is greater than or equal to
    the 100-m isobath, return True since it is *not* on the shelf. If its depth is less than 100-m,
    return False because it is on the shelf.
    
    Input:
        eddy (DataFrame) : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)           : returns true if eddy is deeper than 500m isobath (on shelf), returns false otherwise

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    bathy = xr.open_dataset('/Users/elenaperez/Desktop/rings/data/nwa_bathy.nc')
    
    target_bathy = 100 # desired isobath

    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    bathy_lon_np = bathy.lon.values
    bathy_lon_len = len(bathy_lon_np) # length of bathy lon array

    eddy_lon_np = np.full(bathy_lon_len, eddy_lon)

    min_lon_index = np.argmin(abs(bathy_lon_np-eddy_lon_np)) # index of closest lon in bathy file

    bathy_lat_np = bathy.z.sel(lon=bathy_lon_np[min_lon_index]).lat.values # array of bathy lats at the eddy lon
    bathy_lat_len = len(bathy_lat_np)

    eddy_lat_np = np.full(bathy_lat_len, eddy_lat) 

    min_lat_index = np.argmin(abs(bathy_lat_np-eddy_lat_np)) # index of closest lat in bathy file

    bathy_depth = bathy.z.sel(lat=bathy_lat_np[min_lat_index],lon=bathy_lon_np[min_lon_index]).values # depth of 

    return (target_bathy <= abs(bathy_depth))


#-------------------------------------------------------------------------------------------------------------------------------
# 15)
def count_annual_formations(ring_df, which_zone):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of annual formations for 
    that type of rings (e.g. wcr or ccr)
    
    Input:
        ring_df (DataFrame)              : pandas dataframe of rings
        
    Output:
        ring_annual_count_df (DataFrame) : pandas dataframe of annual ring formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    year_range = np.arange(min(ring_df['time'].dt.year), max(ring_df['time'].dt.year)+1)
    var = ['year', which_zone]

    df_structure = np.zeros((len(year_range), len(var)))
    ring_annual_count_df = pd.DataFrame(df_structure, columns = var)

    counter = 0
    for i in year_range:
        annual_formations = len((ring_df[(ring_df['time'].dt.year == i) & (ring_df['observation_number']==0)])['track'].unique())
        ring_annual_count_df.iloc[counter]=[i, annual_formations]
        counter += 1
        
    return ring_annual_count_df


#-------------------------------------------------------------------------------------------------------------------------------
# 16)
def count_all_formations(ring_df, which_zone):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of monthly formations 
    for that type of ring (e.g. wcr or ccr)
    
    Input:
        ring_df (DataFrame)              : pandas dataframe of rings
        ring_type (String)               : 'wcr' for warm core ring, 'ccr' for cold core ring
        
    Output:
        ring_month_count_df (DataFrame)  : pandas dataframe of monthly ring formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    year_range = np.arange(min(ring_df['time'].dt.year), max(ring_df['time'].dt.year)+1)
    month_range = np.arange(1,13)

    var = ['year','month', which_zone]

    df_structure = np.zeros((len(month_range)*len(year_range), len(var)))
    ring_all_count_df = pd.DataFrame(df_structure, columns = var)

    counter=0
    for year in year_range:
        for month in month_range:
            ring_formations = len((ring_df[(ring_df['time'].dt.month == month) & (ring_df['time'].dt.year == year) & (ring_df['observation_number']==0)])['track'].unique())
            ring_all_count_df.iloc[counter]=[year, month, ring_formations]
            counter += 1

    return ring_all_count_df


#-------------------------------------------------------------------------------------------------------------------------------
# 17)
def eddy_df_to_formation_counts_df():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                      define zone boundaries & create dataframe for each zone                            "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    eddy_wcr_df = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/nwa_wcrlike_eddies.pkl')
    eddy_ccr_dr = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/nwa_ccrlike_eddies.pkl') 
        
    # Gangopadhyay et al., 2019 bounds
    zone_lat = [30,45] # N
    zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]]

    # WCR-like eddies
    zone1_wcrs = eddy_wcr_df[(eddy_wcr_df['longitude']>=zone_lon[1][0]) & (eddy_wcr_df['longitude']<=zone_lon[1][1])]
    zone2_wcrs = eddy_wcr_df[(eddy_wcr_df['longitude']>=zone_lon[2][0]) & (eddy_wcr_df['longitude']<=zone_lon[2][1])]
    zone3_wcrs = eddy_wcr_df[(eddy_wcr_df['longitude']>=zone_lon[3][0]) & (eddy_wcr_df['longitude']<=zone_lon[3][1])]
    zone4_wcrs = eddy_wcr_df[(eddy_wcr_df['longitude']>=zone_lon[4][0]) & (eddy_wcr_df['longitude']<=zone_lon[4][1])]

    # CCR-like eddies
    zone1_ccrs = eddy_ccr_dr[(eddy_ccr_dr['longitude']>=zone_lon[1][0]) & (eddy_ccr_dr['longitude']<=zone_lon[1][1])]
    zone2_ccrs = eddy_ccr_dr[(eddy_ccr_dr['longitude']>=zone_lon[2][0]) & (eddy_ccr_dr['longitude']<=zone_lon[2][1])]
    zone3_ccrs = eddy_ccr_dr[(eddy_ccr_dr['longitude']>=zone_lon[3][0]) & (eddy_ccr_dr['longitude']<=zone_lon[3][1])]
    zone4_ccrs = eddy_ccr_dr[(eddy_ccr_dr['longitude']>=zone_lon[4][0]) & (eddy_ccr_dr['longitude']<=zone_lon[4][1])]
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                  create dictionary with zone name as key and dataframe as value                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  
    zones_wcrs = {'zone1':zone1_wcrs, 'zone2':zone2_wcrs, 'zone3':zone3_wcrs, 'zone4':zone4_wcrs}
    zones_ccrs = {'zone1':zone1_ccrs, 'zone2':zone2_ccrs, 'zone3':zone3_ccrs, 'zone4':zone4_ccrs}

    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                        ANNUAL – count eddy annual formations for each zone                              "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  
    # WCR-like eddies
    zone_wcr_yyyy_formations = count_annual_formations(eddy_wcr_df, 'all_zones')
    for zone in zones_wcrs:
        zone_wcr_yyyy_formations[zone] = count_annual_formations(zones_wcrs[zone], zone)[zone]

    # CCR-like eddies
    zone_ccr_yyyy_formations = count_annual_formations(eddy_ccr_dr, 'all_zones')
    for zone in zones_ccrs:
        zone_ccr_yyyy_formations[zone] = count_annual_formations(zones_ccrs[zone], zone)[zone]
        
    # save zone_wcr_yyyy_formations df as pickled file
    zone_wcr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/zone_wcrlike_yyyy_formations.pkl') 

    # save zone_ccr_yyyy_formations df as pickled file
    zone_ccr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/zone_ccrlike_yyyy_formations.pkl') 

    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                       ALL – count eddy annual & monthly formations for each zone                        "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # WCR-like eddies
    zone_wcr_yyyy_mm_formations = count_all_formations(eddy_wcr_df, 'all_zones')
    for zone in zones_wcrs:
        zone_wcr_yyyy_mm_formations[zone] = count_all_formations(zones_wcrs[zone], zone)[zone]

    # CCR-like eddies
    zone_ccr_yyyy_mm_formations = count_all_formations(eddy_ccr_dr, 'all_zones')
    for zone in zones_ccrs:
        zone_ccr_yyyy_mm_formations[zone] = count_all_formations(zones_ccrs[zone], zone)[zone]
    
    # save zone_wcr_yyyy_mm_formations df as pickled file
    zone_wcr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/zone_wcrlike_yyyy_mm_formations.pkl') 

    # save zone_ccr_yyyy_mm_formations df as pickled file
    zone_ccr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/zone_ccrlike_yyyy_mm_formations.pkl') 


#-------------------------------------------------------------------------------------------------------------------------------


