"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to handle META eddy trajectory data
 
    1) is_wcrlike                           : returns true if given eddy is warm core ring-like, false otherwise
    2) is_ccrlike                           : returns true if given eddy is cold core ring-like, false otherwise
    3) meta2_eddy_to_nwa_ringlike_eddies    : saves pandas df of META2.0 eddy trajectories for wcr-like & ccr-like eddies
    4) get_gs                               : returns overall mean position of the Gulf Stream
    5) get_gs_year                          : takes a given year and returns the average GS path for that year 
    6) get_gs_month                         : takes a given year, month and returns the average GS path for that month
    7) get_gs_day                           : takes a given year, month, day and returns the GS path for that day
    8) get_eddy_formation_loc               : returns lon, lat for a single eddy
    9) get_eddy_formation_time              : returns year, month for a single eddy
   10) get_eddy_lifespan                    : returns lifespan of the eddy as integer
   11) eddy_moves_west                      : returns true is eddy moves westward
   12) eddy_moves_east                      : returns true is eddy moves eastward
   13) is_geq_500m_isobath                  : returns true if eddy depth is greater than or equal to 100-m isobath
   14) get_eddy_demise_loc                  : returns lon, lat for demise location of eddy
   15) count_annual_ring_formations         : returns DataFrame of number of annual eddy formations for an eddy DataFrame
   16) count_all_ring_formations            : returns DataFrame of number of monthly, yearly eddy formations for an eddy DF
   17) eddy_df_to_formation_counts_df       : saves DataFrames of merged formation counts for ring-like eddies
   18) meta31_eddy_to_nwa_ringlike_eddies   : saves pandas df of META3.1exp eddy trajectories for wcr-like & ccr-like eddies
   19) formation_counts_df_to_excel         : saves xlxs files of annual formation counts for regime shift detection
   20) is_wcrlike_gled                      : same as is_wcrlike, but specifically for GLED eddies
   21) is_ccrlike_gled                      : same as is_ccrlike, but specifically for GLED eddies
   22) gled_count_all_formations            : same as count_annual_ring_formations, but specifically for GLED eddies
   23) gled_count_annual_formations         : same as count_all_ring_formations, but specifically for GLED eddies
   24) gled_eddy_df_to_nwa_ringlike_edddies : same as meta2_eddy_to_nwa_ringlike_eddies + count_xx_ring_formations, but for GLED
    
    
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
def is_wcrlike(eddy, gs):
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
        gs (mat)          : .mat file of daily Gulf Stream paths 

    Output:
        (bool)            : returns true if eddy is north of Gulf Stream for given date, else returns false

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    eddy_year = get_eddy_formation_time(eddy)[0] # year of eddy formation
    eddy_month = get_eddy_formation_time(eddy)[1] # month of eddy formation
    eddy_day = get_eddy_formation_time(eddy)[2] # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day, gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day, gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon

    return ((eddy_lat >= gs_lat_np.T[min_index]-0.25) & (is_geq_100m_isobath(eddy) & (eddy['cyclonic_type']==1).all()) & (eddy_lat <= gs_lat_np.T[min_index]+3) & eddy_moves_west(eddy)) # if true then eddy formation is north of gs path


#-------------------------------------------------------------------------------------------------------------------------------
# 2)
def is_ccrlike(eddy, gs):
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
        gs (mat)          : .mat file of daily Gulf Stream paths 

    Output:
        (bool)            : returns true if eddy is south of Gulf Stream for given date, else returns false

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    eddy_year = get_eddy_formation_time(eddy)[0] # year of eddy formation
    eddy_month = get_eddy_formation_time(eddy)[1] # month of eddy formation
    eddy_day = get_eddy_formation_time(eddy)[2] # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day, gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day, gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon
    
    westward_propagation = eddy_moves_west(eddy) # True if eddy moves west
    
    return ((eddy_lat <= gs_lat_np.T[min_index]+0.25) & (is_geq_100m_isobath(eddy)) & (eddy_lat >= gs_lat_np.T[min_index]-3) & ((eddy['cyclonic_type']==-1).all()) & (eddy_lat >= 34) & westward_propagation)


#-------------------------------------------------------------------------------------------------------------------------------
# 3)
def meta2_eddy_to_nwa_ringlike_eddies(path):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    " Function Description:                                                                                   "
    "    The meta2_eddy_to_nwa_ringlike_eddies function reads in META2.0 eddy datasets and converts it to         " 
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
    "                          read in META eddy trajectory datasets & Gulf Stream paths                      "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # read in META2.0 .nc file as xarray dataset
    meta2_ds = xr.open_dataset(path)

    # convert xr dataset to pandas dataframe
    meta2_df = meta2_ds.to_dataframe()

    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/gulf_stream/GS_daily_CMEMS_047_50cm_contours_1993_to_nrt.mat')

    # convert time array to ordinal dates
    for d in range(len(gs['time'][0])-1):
        gs['time'][0][d] = gs['time'][0][d]+date.toordinal(date(1950,1,1))
        
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                        restrict dataframe to nwa box & change lon to -180 to 180                        "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # change longitude from 0 to 360 --> -180 to 180
    meta2_df['longitude'] = ((((meta2_df['longitude'] + 180) % 360) - 180))

    # bounding box for NWA based on Gangopadhyay et al., 2019
    zone_lat = [30,45] # N
    zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]] # W ... 0 –> all zones, 1 –> zone 1, etc.

    # cut eddy trajectory dataframe down to NWA region and time period (1993-2017)
    meta2_df = meta2_df[(meta2_df['longitude'] >= zone_lon[0][0]) & (meta2_df['longitude'] <= zone_lon[0][1]) & (meta2_df['latitude'] >= zone_lat[0]) & (meta2_df['latitude'] <= zone_lat[1]) & (meta2_df['time'].dt.year<2018)]
    
    
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
        if is_wcrlike(eddy, gs).all():
            wcrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_wcrlike_df = anticyclonic_eddy_df[anticyclonic_eddy_df['track'].isin(wcrlike_tracks)]


    # """"""""""""""""""""""""""""""""""""""" ccr-like eddies """""""""""""""""""""""""""""""""""""""""""""""""""
    ccrlike_tracks = []
    # loop through each eddy to determine if it is ccr-like
    for i in np.array(cyclonic_eddy_df['track'].unique()):
        eddy = cyclonic_eddy_df[cyclonic_eddy_df['track']==i]
        if is_ccrlike(eddy, gs).all():
            ccrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_ccrlike_df = cyclonic_eddy_df[cyclonic_eddy_df['track'].isin(ccrlike_tracks)]
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                       save filtered dataframes for wcr- & ccr-like eddies                               "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # save wcr-like df as pickled file
    eddy_wcrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/meta2_nwa_wcrlike_eddies.pkl')
    
    # save ccr-like df as pickled file
    eddy_ccrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/meta2_nwa_ccrlike_eddies.pkl') 


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
def get_gs_day(year, month, day, gs):
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        year (Int)           : the year for which you want to extract GS position (1993-2022)
        month (Int)          : the month for which you want to extract GS position (Jan 1993 to Oct 2022)
        day (Int)            : the day for which you want to extract GS position (Jan 1 1993 to Oct 5 2022)
        gs (mat)             : .mat file of daily Gulf Stream paths 
    
    Output:
        gs_lon_month (array) : 1D array of GS longitude for that month
        gs_lat_month (array) : 1D array of GS latitude for that month
        
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
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
# 14) returns location (lon, lat) of demise for a single eddy

def get_eddy_demise_loc(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns lon/lat location of eddy demise as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    return eddy.iloc[-1].longitude,eddy.iloc[-1].latitude


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
def eddy_df_to_formation_counts_df(whichDataset):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                      define zone boundaries & create dataframe for each zone                            "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 
    eddy_wcr_df = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_nwa_wcrlike_eddies.pkl')
    eddy_ccr_dr = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_nwa_ccrlike_eddies.pkl') 
    
        
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
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                   save files for specified dataset                                      "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
# save zone_wcr/ccr_yyyy_formations df as pickled file
    zone_wcr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_zone_wcrlike_yyyy_formations.pkl')
    zone_ccr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_zone_ccrlike_yyyy_formations.pkl') 

    # save zone_wcr/ccr_yyyy_mm_formations df as pickled file
    zone_wcr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_zone_wcrlike_yyyy_mm_formations.pkl')
    zone_ccr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_zone_ccrlike_yyyy_mm_formations.pkl') 
    

#-------------------------------------------------------------------------------------------------------------------------------
# 18)
def meta31_eddy_to_nwa_ringlike_eddies(path_anticyclonic, path_cyclonic):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    " Function Description:                                                                                   "
    "   The meta31_eddy_to_nwa_ringlike_eddies function reads in META3.1exp eddy datasets and converts it to  " 
    "   a pandas dataframe, filters out eddies that don't qualify as WCR-like or CCR-like, and saves          "
    "   the WCR-like and CCR-like eddy dataframes in the data/dataframes folder.                              "
    "                                                                                                         "
    " Input:                                                                                                  " 
    "    path (String)         : path to where the META eddy dataset is stored                                "
    "                                                                                                         "
    " Output:                                                                                                 "
    "    (None)                : saves pandas dataframe of eddy trajectories for wcr-like and ccr-like eddies "
    "                                                                                                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                      read in META eddy trajectory datasets & Gulf Stream paths                          "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # read in META3.1exp .nc files as xarray dataset
    meta31_a_ds = xr.open_dataset(path_anticyclonic)
    meta31_c_ds = xr.open_dataset(path_cyclonic)

    # change longitude (0-360) --> (-180 to 180)
    meta31_a_ds = meta31_a_ds.assign_coords(longitude=((meta31_a_ds.longitude + 180) % 360) - 180, latitude=meta31_a_ds.latitude)
    meta31_c_ds = meta31_c_ds.assign_coords(longitude=((meta31_c_ds.longitude + 180) % 360) - 180, latitude=meta31_c_ds.latitude)

    # for convience call the Dataet x_ds
    a_ds = meta31_a_ds
    c_ds = meta31_c_ds

    # create cyclonic_type column to put in pandas DataFrame
    anticyclonic_arr = np.full(shape=len(meta31_a_ds.longitude), fill_value=1)
    cyclonic_arr = np.full(shape=len(meta31_c_ds.longitude), fill_value=-1)

    # convert xarray DataSet to pandas DataFrame
    anticyclonic_df = pd.DataFrame({'amplitude': np.array(a_ds.amplitude), 'cyclonic_type': anticyclonic_arr, 'latitude': np.array(a_ds.latitude), 'longitude': np.array(a_ds.longitude), 'observation_flag' :  np.array(a_ds.observation_flag), 'observation_number' : np.array(a_ds.observation_number), 'speed_average' : np.array(a_ds.speed_average), 'speed_radius' : np.array(a_ds.speed_radius), 'time' : np.array(a_ds.time), 'track' : np.array(a_ds.track)})
    cyclonic_df = pd.DataFrame({'amplitude': np.array(c_ds.amplitude), 'cyclonic_type': cyclonic_arr, 'latitude': np.array(c_ds.latitude), 'longitude': np.array(c_ds.longitude), 'observation_flag' :  np.array(c_ds.observation_flag), 'observation_number' : np.array(c_ds.observation_number), 'speed_average' : np.array(c_ds.speed_average), 'speed_radius' : np.array(c_ds.speed_radius), 'time' : np.array(c_ds.time), 'track' : np.array(c_ds.track)})
    
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/rings/data/gulf_stream/GS_daily_CMEMS_047_50cm_contours_1993_to_nrt.mat')

    # convert time array to ordinal dates
    for d in range(len(gs['time'][0])-1):
        gs['time'][0][d] = gs['time'][0][d]+date.toordinal(date(1950,1,1))
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                        restrict dataframe to nwa box & change lon to -180 to 180                        "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # bounding box for NWA based on Gangopadhyay et al., 2019
    zone_lat = [30,45] # N
    zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]] # W ... 0 –> all zones, 1 –> zone 1, etc.

    # cut eddy trajectory dataframe down to NWA region
    anticyclonic_eddy_df = anticyclonic_df[(anticyclonic_df['longitude'] >= zone_lon[0][0]) & (anticyclonic_df['longitude'] <= zone_lon[0][1]) & (anticyclonic_df['latitude'] >= zone_lat[0]) & (anticyclonic_df['latitude'] <= zone_lat[1]) & (anticyclonic_df['time'].dt.year<2018)]
    cyclonic_eddy_df = cyclonic_df[(cyclonic_df['longitude'] >= zone_lon[0][0]) & (cyclonic_df['longitude'] <= zone_lon[0][1]) & (cyclonic_df['latitude'] >= zone_lat[0]) & (cyclonic_df['latitude'] <= zone_lat[1]) & (cyclonic_df['time'].dt.year<2018)]


    """"""""""""""""""""""""""""""""""""""" wcr-like eddies """""""""""""""""""""""""""""""""""""""""""""""""""
    wcrlike_tracks = []
    # loop through each eddy to determine if it is wcr-like
    for i in np.array(anticyclonic_eddy_df['track'].unique()):
        eddy = anticyclonic_eddy_df[anticyclonic_eddy_df['track']==i]
        if is_wcrlike(eddy, gs).all():
            wcrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_wcrlike_df = anticyclonic_eddy_df[anticyclonic_eddy_df['track'].isin(wcrlike_tracks)]


    # """"""""""""""""""""""""""""""""""""""" ccr-like eddies """""""""""""""""""""""""""""""""""""""""""""""""""
    ccrlike_tracks = []
    # loop through each eddy to determine if it is ccr-like
    for i in np.array(cyclonic_eddy_df['track'].unique()):
        eddy = cyclonic_eddy_df[cyclonic_eddy_df['track']==i]
        if is_ccrlike(eddy, gs).all():
            ccrlike_tracks.append(i)

    # only incude eddies that are wcr-like
    eddy_ccrlike_df = cyclonic_eddy_df[cyclonic_eddy_df['track'].isin(ccrlike_tracks)]
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                       save filtered dataframes for wcr- & ccr-like eddies                               "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # save wcr-like df as pickled file
    eddy_wcrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/meta31_nwa_wcrlike_eddies.pkl')
    
    # save ccr-like df as pickled file
    eddy_ccrlike_df.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/meta31_nwa_ccrlike_eddies.pkl') 


#-------------------------------------------------------------------------------------------------------------------------------
# 19) export formation counts to excel to run the Rodinov Regime Shift Detection Software 
def formation_counts_df_to_excel(whichDataset):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    " Function Description:                                                                                   "
    "    The formation_counts_df_to_excel function combines annual formation counts of META dataframes        " 
    "    and saves the counts as excel .xlsx files to be used in the regime shift detection code.             "
    "                                                                                                         "
    " Input:                                                                                                  " 
    "    whichDataset (String) : which META eddy dataset to run the function for                              "
    "                                                                                                         "
    " Output:                                                                                                 "
    "    (None)                : saves .xlsx file of merged annual formation counts                           "
    "                                                                                                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    ccr_formations_df = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_zone_ccrlike_yyyy_formations.pkl') 
    wcr_formations_df = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichDataset+'_zone_wcrlike_yyyy_formations.pkl') 
    
    year_range = np.arange(min(ccr_formations_df['year']), max(ccr_formations_df['year'])+1)
    var = ['year', 'ccr-like', 'wcr-like']

    df_structure = np.zeros((len(year_range), len(var)))
    ring_annual_count_df = pd.DataFrame(df_structure, columns = var)

    ring_annual_count_df['year'] = ccr_formations_df['year']
    ring_annual_count_df['ccr-like'] = ccr_formations_df['all_zones']
    ring_annual_count_df['wcr-like'] = wcr_formations_df['all_zones']
    
    ring_annual_count_df.to_excel('/Users/elenaperez/Desktop/rings/data/excel/'+whichDataset+'_zone_ringlike_yyyy_formations.xlsx')
    


#-------------------------------------------------------------------------------------------------------------------------------

def is_wcrlike_gled(eddy, gs, bathy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Same as is_wcrlike, but altered for GLED eddies, which are a special case of their own.
    
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
        gs (mat)          : .mat file of daily Gulf Stream paths 

    Output:
        (bool)            : returns true if eddy is north of Gulf Stream for given date, else returns false

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    eddy_lon = ((((eddy['longitude'].iloc[0][0])+180)%360)-180) # formation longitude
    eddy_lat = eddy['latitude'].iloc[0][0] # formation latitude

    eddy_year = int((pd.to_datetime(eddy['time'])).dt.year) # year of eddy formation
    eddy_month = int((pd.to_datetime(eddy['time'])).dt.month) # month of eddy formation
    eddy_day = int((pd.to_datetime(eddy['time'])).dt.day) # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year, eddy_month, eddy_day, gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year, eddy_month, eddy_day, gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon

    # greater than or equal to 100m isobath
    target_bathy = 100 # desired isobath
    bathy_lon_np = bathy.lon.values
    bathy_lon_len = len(bathy_lon_np) # length of bathy lon array
    eddy_lon_np = np.full(bathy_lon_len, eddy_lon)
    min_lon_index = np.argmin(abs(bathy_lon_np-eddy_lon_np)) # index of closest lon in bathy file
    bathy_lat_np = bathy.z.sel(lon=bathy_lon_np[min_lon_index]).lat.values # array of bathy lats at the eddy lon
    bathy_lat_len = len(bathy_lat_np)
    eddy_lat_np = np.full(bathy_lat_len, eddy_lat) 
    min_lat_index = np.argmin(abs(bathy_lat_np-eddy_lat_np)) # index of closest lat in bathy file
    bathy_depth = bathy.z.sel(lat=bathy_lat_np[min_lat_index],lon=bathy_lon_np[min_lon_index]).values # depth of 
    geq100m = (target_bathy <= abs(bathy_depth))

    # already checked for westward propagation in intial trim thru of dataset
    # this checks eddy w/in northern bounds of GS path & doesn't go onto shelf
    return ((eddy_lat >= gs_lat_np.T[min_index]-0.25) & (geq100m) & (eddy_lat <= gs_lat_np.T[min_index]+3) & (eddy_lat >= 37))


#-------------------------------------------------------------------------------------------------------------------------------

def is_ccrlike_gled(eddy, gs, bathy):
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Same as is_ccrlike, but altered for GLED eddies, which are a special case of their own.

    
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
        gs (mat)          : .mat file of daily Gulf Stream paths 

    Output:
        (bool)            : returns true if eddy is south of Gulf Stream for given date, else returns false

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    eddy_lon = ((((eddy['longitude'].iloc[0][0])+180)%360)-180) # formation longitude
    eddy_lat = eddy['latitude'].iloc[0][0] # formation latitude

    eddy_year = int((pd.to_datetime(eddy['time'])).dt.year) # year of eddy formation
    eddy_month = int((pd.to_datetime(eddy['time'])).dt.month) # month of eddy formation
    eddy_day = int((pd.to_datetime(eddy['time'])).dt.day) # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year, eddy_month, eddy_day, gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year, eddy_month, eddy_day, gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon

    # greater than or equal to 100m isobath
    target_bathy = 100 # desired isobath
    bathy_lon_np = bathy.lon.values
    bathy_lon_len = len(bathy_lon_np) # length of bathy lon array
    eddy_lon_np = np.full(bathy_lon_len, eddy_lon)
    min_lon_index = np.argmin(abs(bathy_lon_np-eddy_lon_np)) # index of closest lon in bathy file
    bathy_lat_np = bathy.z.sel(lon=bathy_lon_np[min_lon_index]).lat.values # array of bathy lats at the eddy lon
    bathy_lat_len = len(bathy_lat_np)
    eddy_lat_np = np.full(bathy_lat_len, eddy_lat) 
    min_lat_index = np.argmin(abs(bathy_lat_np-eddy_lat_np)) # index of closest lat in bathy file
    bathy_depth = bathy.z.sel(lat=bathy_lat_np[min_lat_index],lon=bathy_lon_np[min_lon_index]).values # depth of 
    geq100m = (target_bathy <= abs(bathy_depth))

    # already checked for westward propagation in intial trim thru of dataset
    # this checks eddy w/in northern bounds of GS path & doesn't go onto shelf
    return ((eddy_lat <= gs_lat_np.T[min_index]+0.25) & (eddy_lat >= gs_lat_np.T[min_index]-3) & (eddy_lat >= 34) & (geq100m))


#-------------------------------------------------------------------------------------------------------------------------------

def gled_count_all_formations(ring_df, which_zone):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of monthly formations 
    for that type of ring (e.g. wcr or ccr)
    
    Input:
        ring_df (DataFrame)              : pandas dataframe of rings
        ring_type (String)               : 'wcr' for warm core ring, 'ccr' for cold core ring
        
    Output:
        ring_month_count_df (DataFrame)  : pandas dataframe of monthly ring formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    year_range = np.arange(1993, 2018)
    month_range = np.arange(1,13)

    var = ['year','month', which_zone]

    df_structure = np.zeros((len(month_range)*len(year_range), len(var)))
    ring_all_count_df = pd.DataFrame(df_structure, columns = var)

    counter=0
    for year in year_range:
        for month in month_range:
            if month<=9: 
                ring_formations = len(ring_df[(ring_df['id'].str.contains(str(year))) & (ring_df['id'].str.contains('-0'+str(month)+'-'))])
            else:
                ring_formations = len(ring_df[(ring_df['id'].str.contains(str(year))) & (ring_df['id'].str.contains('-'+str(month)+'-'))])
            ring_all_count_df.iloc[counter]=[year, month, ring_formations]
            counter += 1

    return ring_all_count_df

#-------------------------------------------------------------------------------------------------------------------------------

def gled_count_annual_formations(ring_df, which_zone):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of annual formations for 
    that type of rings (e.g. wcr or ccr)
    
    Input:
        ring_df (DataFrame)              : pandas dataframe of rings
        
    Output:
        ring_annual_count_df (DataFrame) : pandas dataframe of annual ring formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    year_range = np.arange(1993, 2018)
    var = ['year', which_zone]

    df_structure = np.zeros((len(year_range), len(var)))
    ring_annual_count_df = pd.DataFrame(df_structure, columns = var)

    counter = 0
    for i in year_range:
        annual_formations = len(ring_df[ring_df['id'].str.contains(str(i))])
        ring_annual_count_df.iloc[counter]=[i, annual_formations]
        counter += 1
        
    return ring_annual_count_df

#-------------------------------------------------------------------------------------------------------------------------------

def gled_eddy_df_to_nwa_ringlike_edddies(path30d, path90d, path180d):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    " Function Description:                                                                                   "
    "    The ??????? function reads in GLED eddy dataset and converts it to                                   " 
    "    a pandas dataframe, filters out eddies that don't qualify as WCR-like or CCR-like, and saves         "
    "    the WCR-like and CCR-like eddy dataframes in the data/dataframes folder.                             "
    "                                                                                                         "
    " Input:                                                                                                  " 
    "    path (String)         : path to where the GLED eddy dataset is stored                                "
    "                                                                                                         "
    " Output:                                                                                                 "
    "    (None)                : saves pandas dataframe of eddy trajectories for wcr-like and ccr-like eddies "
    "                                                                                                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                               read in GLED eddy trajectory datasets                                     "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # read in GLED files as pandas dataframe
    gled_df30 = pd.read_json(path30d)
    gled_df90 = pd.read_json(path90d)
    gled_df180 = pd.read_json(path180d)

    # define lat/lon box of NWA, Gangopadhyay et al., 2019 bounds
    zone_lon = [[285,305],[285,290],[290,295],[295,300],[300,305]] # W
    zone_lat = [30,45] # N

    # which dataset

    Files = {'gled30d': gled_df30, 'gled90d': gled_df90,'gled180d': gled_df180}

    for whichFile in Files:
        gled_df = Files[whichFile]

        # filter to lat eddy start is in lat/lon box and net movement is westward
        bbox = []
        # checks eddy start is in lat/lon box and net movement is westward
        for i in np.arange(len(gled_df['center_lon'])):
            if ((gled_df['center_lon'][i][0]>=zone_lon[0][0]) & (gled_df['center_lon'][i][0]<=zone_lon[0][1]) & (gled_df['center_lat'][i][0] >= zone_lat[0]) & (gled_df['center_lat'][i][0] <= zone_lat[1]) & (int(gled_df['date_start'][i][0:4])<2018) & (gled_df['center_lon'][i][0] > gled_df['center_lon'][i][-1])):
                bbox.append(i)

        # nwa eddies
        gled_df = gled_df.iloc[bbox]

        # anticyclonic
        a_gled_df = gled_df[gled_df['cyc']==1]
        # cyclonic
        c_gled_df = gled_df[gled_df['cyc']==-1]

        # create pandas dataframe similar to META eddy dataframe

        # create cyclonic_type column to put in pandas DataFrame
        anticyclonic_arr = np.full(shape=len(a_gled_df.center_lon), fill_value=1)
        cyclonic_arr = np.full(shape=len(c_gled_df.center_lon), fill_value=-1)

        # to make the following lines not as long
        a_df = a_gled_df
        c_df = c_gled_df

        # change into more META type format
        a_gled_df = pd.DataFrame({'id': np.array(a_df.id), 'time' : np.array(a_df.date_start), 'duration' : np.array(a_df.duration), 'radius' : np.array(a_df.radius),'cyclonic_type': np.array(a_df.cyc), 'longitude': np.array(a_df.center_lon), 'latitude': np.array(a_df.center_lat), 'dx' :  np.array(a_df.dx), 'speed_x' : np.array(a_df.speed_x), 'dy' : np.array(a_df.dy), 'speed_y' : np.array(a_df.speed_y), 'vort' : np.array(a_df.vort), 'lavd' : np.array(a_df.lavd)})
        c_gled_df = pd.DataFrame({'id': np.array(c_df.id), 'time' : np.array(c_df.date_start), 'duration' : np.array(c_df.duration), 'radius' : np.array(c_df.radius),'cyclonic_type': np.array(c_df.cyc), 'longitude': np.array(c_df.center_lon), 'latitude': np.array(c_df.center_lat), 'dx' :  np.array(c_df.dx), 'speed_x' : np.array(c_df.speed_x), 'dy' : np.array(c_df.dy), 'speed_y' : np.array(c_df.speed_y), 'vort' : np.array(c_df.vort), 'lavd' : np.array(c_df.lavd)})

        #### FILTER
        # filter WCR like eddies
        wcrlike_tracks = []
        ind = (np.array(a_gled_df['id']))
        for i in np.arange(len(a_gled_df['id'])):
            eddy = a_gled_df[a_gled_df['id']==ind[i]]
            if is_wcrlike_gled(eddy, gs, bathy).all():
                wcrlike_tracks.append(i)

        # filter CCR like eddies
        ccrlike_tracks = []
        ind = (np.array(c_gled_df['id']))
        for i in np.arange(len(c_gled_df['id'])):
            eddy = c_gled_df[c_gled_df['id']==ind[i]]
            if is_ccrlike_gled(eddy, gs, bathy).all():
                ccrlike_tracks.append(i)

        # BY ZONE

        # gled_zone_wcrlike_gled_eddy_df = a_gled_df30.iloc[wcrlike_tracks]
        zone_wcrlike_eddies = a_gled_df.iloc[wcrlike_tracks]       
        zone_ccrlike_eddies = c_gled_df.iloc[ccrlike_tracks]       

        # save eddy file
        zone_wcrlike_eddies.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_nwa_wcrlike_eddies.pkl') 
        zone_ccrlike_eddies.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_nwa_ccrlike_eddies.pkl') 


        eddy_df = zone_wcrlike_eddies
        zone1_wcrs_arr = []
        zone2_wcrs_arr = []
        zone3_wcrs_arr = []
        zone4_wcrs_arr = []
        # checks eddy start is in lat/lon box and net movement is westward
        for i in np.arange(len(eddy_df['longitude'])):
            # zone 1
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[1][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[1][1])):
                zone1_wcrs_arr.append(i)
            # zone 2
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[2][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[2][1])):
                zone2_wcrs_arr.append(i)
            # zone 3
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[3][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[3][1])):
                zone3_wcrs_arr.append(i)
            # zone 4
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[4][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[4][1])):
                zone4_wcrs_arr.append(i)


        eddy_df = zone_ccrlike_eddies
        zone1_ccrs_arr = []
        zone2_ccrs_arr = []
        zone3_ccrs_arr = []
        zone4_ccrs_arr = []
        # checks eddy start is in lat/lon box and net movement is westward
        for i in np.arange(len(eddy_df['longitude'])):
            # zone 1
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[1][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[1][1])):
                zone1_ccrs_arr.append(i)
            # zone 2
            elif ((eddy_df['longitude'].iloc[i][0]>=zone_lon[2][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[2][1])):
                zone2_ccrs_arr.append(i)
            # zone 3
            elif ((eddy_df['longitude'].iloc[i][0]>=zone_lon[3][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[3][1])):
                zone3_ccrs_arr.append(i)
            # zone 4
            elif ((eddy_df['longitude'].iloc[i][0]>=zone_lon[4][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[4][1])):
                zone4_ccrs_arr.append(i)

        zone1_wcrs = zone_wcrlike_eddies.iloc[zone1_wcrs_arr]
        zone2_wcrs = zone_wcrlike_eddies.iloc[zone2_wcrs_arr]
        zone3_wcrs = zone_wcrlike_eddies.iloc[zone3_wcrs_arr]
        zone4_wcrs = zone_wcrlike_eddies.iloc[zone4_wcrs_arr]

        # cyclones
        zone1_ccrs = zone_ccrlike_eddies.iloc[zone1_ccrs_arr]
        zone2_ccrs = zone_ccrlike_eddies.iloc[zone2_ccrs_arr]
        zone3_ccrs = zone_ccrlike_eddies.iloc[zone3_ccrs_arr]
        zone4_ccrs = zone_ccrlike_eddies.iloc[zone4_ccrs_arr]


        ## counts

        # WCR LIKE

        zones_wcrs = {'zone1':zone1_wcrs, 'zone2':zone2_wcrs, 'zone3':zone3_wcrs, 'zone4':zone4_wcrs}
        # # zones_ccrs = {'zone1':zone1_ccrs, 'zone2':zone2_ccrs, 'zone3':zone3_ccrs, 'zone4':zone4_ccrs}

        zone_wcr_yyyy_formations = gled_count_annual_formations(zone_wcrlike_eddies, 'all_zones')
        for zone in zones_wcrs:
            zone_wcr_yyyy_formations[zone] = gled_count_annual_formations(zones_wcrs[zone], zone)[zone]

        # save zone_wcr_yyyy_formations df as pickled file
        zone_wcr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_wcrlike_yyyy_formations.pkl') 

        # yearly, monthly formations wcrs
        zone_wcr_yyyy_mm_formations = gled_count_all_formations(zone_wcrlike_eddies, 'all_zones')
        for zone in zones_wcrs:
            zone_wcr_yyyy_mm_formations[zone] = gled_count_all_formations(zones_wcrs[zone], zone)[zone]
        zone_wcr_yyyy_mm_formations

        # save zone_wcr_yyyy_formations df as pickled file
        zone_wcr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_wcrlike_yyyy_mm_formations.pkl') 



        ### CCR LIKE
        zones_ccrs = {'zone1':zone1_ccrs, 'zone2':zone2_ccrs, 'zone3':zone3_ccrs, 'zone4':zone4_ccrs}

        zone_ccr_yyyy_formations = gled_count_annual_formations(zone_ccrlike_eddies, 'all_zones')
        for zone in zones_ccrs:
            zone_ccr_yyyy_formations[zone] = gled_count_annual_formations(zones_ccrs[zone], zone)[zone]

        # save zone_ccr_yyyy_formations df as pickled file
        zone_ccr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_ccrlike_yyyy_formations.pkl') 

        # yearly, monthly formations ccrs
        zone_ccr_yyyy_mm_formations = gled_count_all_formations(zone_ccrlike_eddies, 'all_zones')
        for zone in zones_ccrs:
            zone_ccr_yyyy_mm_formations[zone] = gled_count_all_formations(zones_ccrs[zone], zone)[zone]
        zone_ccr_yyyy_mm_formations

        # save zone_ccr_yyyy_formations df as pickled file
        zone_ccr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_ccrlike_yyyy_mm_formations.pkl') 


    # end of loop

    # load all the counts dataframes
    # GLED 30d
    # annual formation counts by zone for zone_wcrs/ccrs
    gled30d_zone_wcrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_wcrlike_yyyy_formations.pkl') 
    gled30d_zone_ccrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_ccrlike_yyyy_formations.pkl') 

    # monthly formation counts by zone for zone_wcrs/ccrs
    gled30d_zone_wcrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled30d_zone_ccrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_ccrlike_yyyy_mm_formations.pkl') 

    # GLED 90d
    # annual formation counts by zone for zone_wcrs/ccrs
    gled90d_zone_wcrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_wcrlike_yyyy_formations.pkl') 
    gled90d_zone_ccrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_ccrlike_yyyy_formations.pkl') 

    # monthly formation counts by zone for zone_wcrs/ccrs
    gled90d_zone_wcrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled90d_zone_ccrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_ccrlike_yyyy_mm_formations.pkl') 

    # GLED 180d
    # annual formation counts by zone for zone_wcrs/ccrs
    gled180d_zone_wcrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_wcrlike_yyyy_formations.pkl') 
    gled180d_zone_ccrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_ccrlike_yyyy_formations.pkl') 

    # monthly formation counts by zone for zone_wcrs/ccrs
    gled180d_zone_wcrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled180d_zone_ccrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_ccrlike_yyyy_mm_formations.pkl') 

    # merge counts  dataframes
    gled_zone_wcrlike_yyyy_formations = gled30d_zone_wcrlike_yyyy_formations + gled90d_zone_wcrlike_yyyy_formations + gled180d_zone_wcrlike_yyyy_formations
    gled_zone_ccrlike_yyyy_formations = gled30d_zone_ccrlike_yyyy_formations + gled90d_zone_ccrlike_yyyy_formations + gled180d_zone_ccrlike_yyyy_formations
    # correct year
    gled_zone_wcrlike_yyyy_formations['year'] = gled_zone_wcrlike_yyyy_formations['year']/3
    gled_zone_ccrlike_yyyy_formations['year'] = gled_zone_ccrlike_yyyy_formations['year']/3

    gled_zone_wcrlike_yyyy_mm_formations = gled30d_zone_wcrlike_yyyy_mm_formations + gled90d_zone_wcrlike_yyyy_mm_formations + gled180d_zone_wcrlike_yyyy_mm_formations
    gled_zone_ccrlike_yyyy_mm_formations = gled30d_zone_ccrlike_yyyy_mm_formations + gled90d_zone_ccrlike_yyyy_mm_formations + gled180d_zone_ccrlike_yyyy_mm_formations
    # correct year
    gled_zone_wcrlike_yyyy_mm_formations['year'] = gled_zone_wcrlike_yyyy_mm_formations['year']/3
    gled_zone_ccrlike_yyyy_mm_formations['year'] = gled_zone_ccrlike_yyyy_mm_formations['year']/3
    gled_zone_wcrlike_yyyy_mm_formations['month'] = gled_zone_wcrlike_yyyy_mm_formations['month']/3
    gled_zone_ccrlike_yyyy_mm_formations['month'] = gled_zone_ccrlike_yyyy_mm_formations['month']/3


    # save files again
    # save zone_ccr_yyyy_formations df as pickled file
    gled_zone_wcrlike_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled_zone_ccrlike_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_ccrlike_yyyy_mm_formations.pkl') 

    gled_zone_wcrlike_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_wcrlike_yyyy_formations.pkl') 
    gled_zone_ccrlike_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_ccrlike_yyyy_formations.pkl') 


    gled30d_zone_ccrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_nwa_ccrlike_eddies.pkl')
    gled30d_zone_wcrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_nwa_wcrlike_eddies.pkl') 

    gled90d_zone_ccrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_nwa_ccrlike_eddies.pkl')
    gled90d_zone_wcrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_nwa_wcrlike_eddies.pkl') 

    gled180d_zone_ccrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_nwa_ccrlike_eddies.pkl')
    gled180d_zone_wcrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_nwa_wcrlike_eddies.pkl') 

#-------------------------------------------------------------------------------------------------------------------------------

def gled_eddy_df_to_nwa_ringlike_edddies_WHOLE_NWA(path30d, path90d, path180d):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    " Function Description:                                                                                   "
    "    The ??????? function reads in GLED eddy dataset and converts it to                                   " 
    "    a pandas dataframe, filters out eddies that don't qualify as WCR-like or CCR-like, and saves         "
    "    the WCR-like and CCR-like eddy dataframes in the data/dataframes folder.                             "
    "                                                                                                         "
    " Input:                                                                                                  " 
    "    path (String)         : path to where the GLED eddy dataset is stored                                "
    "                                                                                                         "
    " Output:                                                                                                 "
    "    (None)                : saves pandas dataframe of eddy trajectories for wcr-like and ccr-like eddies "
    "                                                                                                         "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                               read in GLED eddy trajectory datasets                                     "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # read in GLED files as pandas dataframe
    gled_df30 = pd.read_json(path30d)
    gled_df90 = pd.read_json(path90d)
    gled_df180 = pd.read_json(path180d)

    # define lat/lon box of NWA, Gangopadhyay et al., 2019 bounds
    zone_lon = [[285,305],[285,290],[290,295],[295,300],[300,305]] # W
    zone_lat = [30,45] # N

    # which dataset

    Files = {'gled30d': gled_df30, 'gled90d': gled_df90,'gled180d': gled_df180}

    for whichFile in Files:
        gled_df = Files[whichFile]

        # filter to lat eddy start is in lat/lon box and net movement is westward
        bbox = []
        # checks eddy start is in lat/lon box and net movement is westward
        for i in np.arange(len(gled_df['center_lon'])):
            if ((gled_df['center_lon'][i][0]>=zone_lon[0][0]) & (gled_df['center_lon'][i][0]<=zone_lon[0][1]) & (gled_df['center_lat'][i][0] >= zone_lat[0]) & (gled_df['center_lat'][i][0] <= zone_lat[1]) & (int(gled_df['date_start'][i][0:4])<2018) & (gled_df['center_lon'][i][0] > gled_df['center_lon'][i][-1])):
                bbox.append(i)

        # nwa eddies
        gled_df = gled_df.iloc[bbox]

        # anticyclonic
        a_gled_df = gled_df[gled_df['cyc']==1]
        # cyclonic
        c_gled_df = gled_df[gled_df['cyc']==-1]

        # create pandas dataframe similar to META eddy dataframe

        # create cyclonic_type column to put in pandas DataFrame
        anticyclonic_arr = np.full(shape=len(a_gled_df.center_lon), fill_value=1)
        cyclonic_arr = np.full(shape=len(c_gled_df.center_lon), fill_value=-1)

        # to make the following lines not as long
        a_df = a_gled_df
        c_df = c_gled_df

        # change into more META type format
        a_gled_df = pd.DataFrame({'id': np.array(a_df.id), 'time' : np.array(a_df.date_start), 'duration' : np.array(a_df.duration), 'radius' : np.array(a_df.radius),'cyclonic_type': np.array(a_df.cyc), 'longitude': np.array(a_df.center_lon), 'latitude': np.array(a_df.center_lat), 'dx' :  np.array(a_df.dx), 'speed_x' : np.array(a_df.speed_x), 'dy' : np.array(a_df.dy), 'speed_y' : np.array(a_df.speed_y), 'vort' : np.array(a_df.vort), 'lavd' : np.array(a_df.lavd)})
        c_gled_df = pd.DataFrame({'id': np.array(c_df.id), 'time' : np.array(c_df.date_start), 'duration' : np.array(c_df.duration), 'radius' : np.array(c_df.radius),'cyclonic_type': np.array(c_df.cyc), 'longitude': np.array(c_df.center_lon), 'latitude': np.array(c_df.center_lat), 'dx' :  np.array(c_df.dx), 'speed_x' : np.array(c_df.speed_x), 'dy' : np.array(c_df.dy), 'speed_y' : np.array(c_df.speed_y), 'vort' : np.array(c_df.vort), 'lavd' : np.array(c_df.lavd)})

    #     #### FILTER
    #     # filter WCR like eddies
    #     wcrlike_tracks = []
    #     ind = (np.array(a_gled_df['id']))
    #     for i in np.arange(len(a_gled_df['id'])):
    #         eddy = a_gled_df[a_gled_df['id']==ind[i]]
    #         if is_wcrlike_gled(eddy, gs, bathy).all():
    #             wcrlike_tracks.append(i)

    #     # filter CCR like eddies
    #     ccrlike_tracks = []
    #     ind = (np.array(c_gled_df['id']))
    #     for i in np.arange(len(c_gled_df['id'])):
    #         eddy = c_gled_df[c_gled_df['id']==ind[i]]
    #         if is_ccrlike_gled(eddy, gs, bathy).all():
    #             ccrlike_tracks.append(i)

    #     # BY ZONE

    #     # gled_zone_wcrlike_gled_eddy_df = a_gled_df30.iloc[wcrlike_tracks]
    #     zone_wcrlike_eddies = a_gled_df.iloc[wcrlike_tracks]       
    #     zone_ccrlike_eddies = c_gled_df.iloc[ccrlike_tracks]       

        # gled_zone_wcrlike_gled_eddy_df = a_gled_df30.iloc[wcrlike_tracks]
        zone_wcrlike_eddies = a_gled_df       
        zone_ccrlike_eddies = c_gled_df  

        # save eddy file
        zone_wcrlike_eddies.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_nwa_wcrlike_eddies.pkl') 
        zone_ccrlike_eddies.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_nwa_ccrlike_eddies.pkl') 


        eddy_df = zone_wcrlike_eddies
        zone1_wcrs_arr = []
        zone2_wcrs_arr = []
        zone3_wcrs_arr = []
        zone4_wcrs_arr = []
        # checks eddy start is in lat/lon box and net movement is westward
        for i in np.arange(len(eddy_df['longitude'])):
            # zone 1
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[1][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[1][1])):
                zone1_wcrs_arr.append(i)
            # zone 2
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[2][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[2][1])):
                zone2_wcrs_arr.append(i)
            # zone 3
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[3][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[3][1])):
                zone3_wcrs_arr.append(i)
            # zone 4
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[4][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[4][1])):
                zone4_wcrs_arr.append(i)


        eddy_df = zone_ccrlike_eddies
        zone1_ccrs_arr = []
        zone2_ccrs_arr = []
        zone3_ccrs_arr = []
        zone4_ccrs_arr = []
        # checks eddy start is in lat/lon box and net movement is westward
        for i in np.arange(len(eddy_df['longitude'])):
            # zone 1
            if ((eddy_df['longitude'].iloc[i][0]>=zone_lon[1][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[1][1])):
                zone1_ccrs_arr.append(i)
            # zone 2
            elif ((eddy_df['longitude'].iloc[i][0]>=zone_lon[2][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[2][1])):
                zone2_ccrs_arr.append(i)
            # zone 3
            elif ((eddy_df['longitude'].iloc[i][0]>=zone_lon[3][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[3][1])):
                zone3_ccrs_arr.append(i)
            # zone 4
            elif ((eddy_df['longitude'].iloc[i][0]>=zone_lon[4][0]) & (eddy_df['longitude'].iloc[i][0]<=zone_lon[4][1])):
                zone4_ccrs_arr.append(i)

        zone1_wcrs = zone_wcrlike_eddies.iloc[zone1_wcrs_arr]
        zone2_wcrs = zone_wcrlike_eddies.iloc[zone2_wcrs_arr]
        zone3_wcrs = zone_wcrlike_eddies.iloc[zone3_wcrs_arr]
        zone4_wcrs = zone_wcrlike_eddies.iloc[zone4_wcrs_arr]

        # cyclones
        zone1_ccrs = zone_ccrlike_eddies.iloc[zone1_ccrs_arr]
        zone2_ccrs = zone_ccrlike_eddies.iloc[zone2_ccrs_arr]
        zone3_ccrs = zone_ccrlike_eddies.iloc[zone3_ccrs_arr]
        zone4_ccrs = zone_ccrlike_eddies.iloc[zone4_ccrs_arr]


        ## counts

        # WCR LIKE

        zones_wcrs = {'zone1':zone1_wcrs, 'zone2':zone2_wcrs, 'zone3':zone3_wcrs, 'zone4':zone4_wcrs}
        # # zones_ccrs = {'zone1':zone1_ccrs, 'zone2':zone2_ccrs, 'zone3':zone3_ccrs, 'zone4':zone4_ccrs}

        zone_wcr_yyyy_formations = gled_count_annual_formations(zone_wcrlike_eddies, 'all_zones')
        for zone in zones_wcrs:
            zone_wcr_yyyy_formations[zone] = gled_count_annual_formations(zones_wcrs[zone], zone)[zone]

        # save zone_wcr_yyyy_formations df as pickled file
        zone_wcr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_wcrlike_yyyy_formations.pkl') 

        # yearly, monthly formations wcrs
        zone_wcr_yyyy_mm_formations = gled_count_all_formations(zone_wcrlike_eddies, 'all_zones')
        for zone in zones_wcrs:
            zone_wcr_yyyy_mm_formations[zone] = gled_count_all_formations(zones_wcrs[zone], zone)[zone]
        zone_wcr_yyyy_mm_formations

        # save zone_wcr_yyyy_formations df as pickled file
        zone_wcr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_wcrlike_yyyy_mm_formations.pkl') 



        ### CCR LIKE
        zones_ccrs = {'zone1':zone1_ccrs, 'zone2':zone2_ccrs, 'zone3':zone3_ccrs, 'zone4':zone4_ccrs}

        zone_ccr_yyyy_formations = gled_count_annual_formations(zone_ccrlike_eddies, 'all_zones')
        for zone in zones_ccrs:
            zone_ccr_yyyy_formations[zone] = gled_count_annual_formations(zones_ccrs[zone], zone)[zone]

        # save zone_ccr_yyyy_formations df as pickled file
        zone_ccr_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_ccrlike_yyyy_formations.pkl') 

        # yearly, monthly formations ccrs
        zone_ccr_yyyy_mm_formations = gled_count_all_formations(zone_ccrlike_eddies, 'all_zones')
        for zone in zones_ccrs:
            zone_ccr_yyyy_mm_formations[zone] = gled_count_all_formations(zones_ccrs[zone], zone)[zone]
        zone_ccr_yyyy_mm_formations

        # save zone_ccr_yyyy_formations df as pickled file
        zone_ccr_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/'+whichFile+'_zone_ccrlike_yyyy_mm_formations.pkl') 


    # end of loop

    # load all the counts dataframes
    # GLED 30d
    # annual formation counts by zone for zone_wcrs/ccrs
    gled30d_zone_wcrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_wcrlike_yyyy_formations.pkl') 
    gled30d_zone_ccrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_ccrlike_yyyy_formations.pkl') 

    # monthly formation counts by zone for zone_wcrs/ccrs
    gled30d_zone_wcrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled30d_zone_ccrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_zone_ccrlike_yyyy_mm_formations.pkl') 

    # GLED 90d
    # annual formation counts by zone for zone_wcrs/ccrs
    gled90d_zone_wcrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_wcrlike_yyyy_formations.pkl') 
    gled90d_zone_ccrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_ccrlike_yyyy_formations.pkl') 

    # monthly formation counts by zone for zone_wcrs/ccrs
    gled90d_zone_wcrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled90d_zone_ccrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_zone_ccrlike_yyyy_mm_formations.pkl') 

    # GLED 180d
    # annual formation counts by zone for zone_wcrs/ccrs
    gled180d_zone_wcrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_wcrlike_yyyy_formations.pkl') 
    gled180d_zone_ccrlike_yyyy_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_ccrlike_yyyy_formations.pkl') 

    # monthly formation counts by zone for zone_wcrs/ccrs
    gled180d_zone_wcrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled180d_zone_ccrlike_yyyy_mm_formations = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_zone_ccrlike_yyyy_mm_formations.pkl') 

    # merge counts  dataframes
    gled_zone_wcrlike_yyyy_formations = gled30d_zone_wcrlike_yyyy_formations + gled90d_zone_wcrlike_yyyy_formations + gled180d_zone_wcrlike_yyyy_formations
    gled_zone_ccrlike_yyyy_formations = gled30d_zone_ccrlike_yyyy_formations + gled90d_zone_ccrlike_yyyy_formations + gled180d_zone_ccrlike_yyyy_formations
    # correct year
    gled_zone_wcrlike_yyyy_formations['year'] = gled_zone_wcrlike_yyyy_formations['year']/3
    gled_zone_ccrlike_yyyy_formations['year'] = gled_zone_ccrlike_yyyy_formations['year']/3

    gled_zone_wcrlike_yyyy_mm_formations = gled30d_zone_wcrlike_yyyy_mm_formations + gled90d_zone_wcrlike_yyyy_mm_formations + gled180d_zone_wcrlike_yyyy_mm_formations
    gled_zone_ccrlike_yyyy_mm_formations = gled30d_zone_ccrlike_yyyy_mm_formations + gled90d_zone_ccrlike_yyyy_mm_formations + gled180d_zone_ccrlike_yyyy_mm_formations
    # correct year
    gled_zone_wcrlike_yyyy_mm_formations['year'] = gled_zone_wcrlike_yyyy_mm_formations['year']/3
    gled_zone_ccrlike_yyyy_mm_formations['year'] = gled_zone_ccrlike_yyyy_mm_formations['year']/3
    gled_zone_wcrlike_yyyy_mm_formations['month'] = gled_zone_wcrlike_yyyy_mm_formations['month']/3
    gled_zone_ccrlike_yyyy_mm_formations['month'] = gled_zone_ccrlike_yyyy_mm_formations['month']/3


    # save files again
    # save zone_ccr_yyyy_formations df as pickled file
    gled_zone_wcrlike_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_wcrlike_yyyy_mm_formations.pkl') 
    gled_zone_ccrlike_yyyy_mm_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_ccrlike_yyyy_mm_formations.pkl') 

    gled_zone_wcrlike_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_wcrlike_yyyy_formations.pkl') 
    gled_zone_ccrlike_yyyy_formations.to_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled_zone_ccrlike_yyyy_formations.pkl') 


    # open eddy files

    gled30d_zone_ccrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_nwa_ccrlike_eddies.pkl')
    gled30d_zone_wcrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled30d_nwa_wcrlike_eddies.pkl') 

    gled90d_zone_ccrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_nwa_ccrlike_eddies.pkl')
    gled90d_zone_wcrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled90d_nwa_wcrlike_eddies.pkl') 

    gled180d_zone_ccrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_nwa_ccrlike_eddies.pkl')
    gled180d_zone_wcrlike_eddies = pd.read_pickle('/Users/elenaperez/Desktop/rings/data/dataframes/gled180d_nwa_wcrlike_eddies.pkl') 

    formation_counts_df_to_excel('gled')
    
    
#-------------------------------------------------------------------------------------------------------------------------------

def ringCensus_seasonality_lineplot(wcr_formations_df, title, ylim_min, ylim_max, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones                   "
    "    wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones                   "
    "    title (Str)                   : title of the figure, e.g. 'Map of the Northwest Atlantic'          "
    "    fig_quality (Int)             : integer of what dpi the image will be set to (e.g., 100 dpi)       "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a bar plot for distribution of formations by zone                                        "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # calculate average monthly formations
    # initialize
    wcr_monthly_std = [] # WCRS
    wcr_monthly_avg = []
    ccr_monthly_std = [] # CCRs
    ccr_monthly_avg = []
    # compute means & std
    for i in np.arange(1,13):
        wcr_monthly_std.append((wcr_formations_df[wcr_formations_df['month']==i]['all_zones']).std())
        wcr_monthly_avg.append(wcr_formations_df[wcr_formations_df['month']==i]['all_zones'].mean())
        
    # lines for standard deviation
    wcr_minus_std = np.subtract(wcr_monthly_avg,wcr_monthly_std)
    wcr_plus_std = np.add(wcr_monthly_avg,wcr_monthly_std)

    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle(title+'('+str(int(min(wcr_formations_df['year'])))+' - '+str(int(max(wcr_formations_df['year'])))+')', y=0.935, fontsize=14);

    # plot
    ax.plot(wcr_formations_df['month'].drop_duplicates(),wcr_monthly_avg,'-o',color='#EE4D4D');
    ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_plus_std, wcr_monthly_avg, alpha=0.3, color='#EE4D4D')
    ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_minus_std, wcr_monthly_avg, alpha=0.3, color='#EE4D4D')
    
    zones = ['zone_1', 'zone_2', 'zone_3', 'zone_4']

    # axes formatting
    ax.set_xlabel('Months', fontweight='bold')
    ax.set_ylabel('Average No. of Formations',fontweight='bold');
    ax.set_ylim(ylim_min,ylim_max)
    ax.set_xlim(1,12)

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#EE4D4D', lw=4), Line2D([0], [0], color='#388FEF', lw=4)]
    ax.legend(custom_lines, ['WCR'], loc='upper left');
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
