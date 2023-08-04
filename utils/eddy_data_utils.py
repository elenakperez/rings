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
    ) is_geq_100m_isobath            : returns true if eddy depth is greater than or equal to 100-m isobath
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


#-------------------------------------------------------------------------------------------------------------------------------------------
################# MIGHT DELETE AFTER I RE-RUN for Chelton 1993-2020 dataset #################
# ) !OLD! create a pandas dataframe for eddies in the nwa 

def eddy_df_old_version():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:

        
    Output:
        * saves pandas dataframe as pickle of eddies in nwa based off nwa_tracks
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # open dataset from path
    ds = nc.Dataset('/Users/elenaperez/Desktop/chatts/data/eddy_trajectory_dt_2.0_19930101_20200307.nc')

    # get length of a column (number of observations=27880804)
    lat =ds['latitude'][:].data

    # create empty matrix to store all the observations and related variables 
    matrix = np.zeros((len(lat), len(np.array(f.variables))))

    # variables that will be in the final dataframe
    # var = ['track','latitude','longitude','time','cyclonic_type']
    var = np.array(f.variables)

    # for loop to fill empty matrix with variables' values 
    counter = 0 
    for v in var:
        d = ds[v][:].data
        matrix[:, counter] = d
        counter += 1

    # create a dataframe just for time column (because np doesn't preserve datetime objects)
    df_time = f.time[:].to_dataframe()

    # create dataframe from eddy matrix
    df = pd.DataFrame(matrix, columns = var)

    # replace time column with time dataframe
    df['time'] = df_time['time']


    # filter out eddy tracks that are not in northwest atlantic
    eddy_df = df[df['track'].isin(nwa_tracks)]

    # to save a new pickled file
    # eddy_df.to_pickle("nwa_eddies.pkl") #may need to change name

    # check to see if it worked
    eddy_df

 #-------------------------------------------------------------------------------------------------------------------------------------------
################# MIGHT DELETE AFTER I RE-RUN for Chelton 1993-2020 dataset #################
# ) !OLD!

def make_eddy_df_old_version(eddy_df, eddy_tracks, time_period, bbox_lon, bbox_lat, df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy_df (DataFrame) : pandas dataframe of eddies in the larger region, e.g. northwest atlantic
        eddy_tracks (array) : tracks of eddies in the larger region, e.g. nwa_tracks = [4, 5, 7, ...]
        time_period (array) : time period, e.g. ['1994-01-01','1994-12-31']. False if no time period
        bbox_lon (array)    : upper and lower longitude bounds of box, e.g. [-75,-70] W
        bbox_lat (array)    : upper and lower latitude bounds of box, e.g. [30,45] N
        df_name (string)    : name of file to be saved 

        
    Output:
        new_df (DataFrame)  : pandas dataframe of eddies contained in the specified box and time period
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # loop through eddy tracks to see if eddy is formed in the specified box
    zone_tracks = []
    for i in eddy_tracks:
        if ((eddy_df[eddy_df['track']==i].iloc[0].longitude >= bbox_lon[0]) & (eddy_df[eddy_df['track']==i].iloc[0].longitude <= bbox_lon[1]) & (eddy_df[eddy_df['track']==i].iloc[0].latitude>= bbox_lat[0]) & (eddy_df[eddy_df['track']==i].iloc[0].latitude<= bbox_lat[1])):
            zone_tracks.append(i)
      
    # cut to specified bounding box
    new_df = eddy_df[eddy_df['track'].isin(zone_tracks)] 
    
    # cut to specified time period
    if time_period!=False:
        new_df = new_df[(new_df['time'] >= time_period[0]) & (new_df['time'] <= time_period[0])] 
    
    # save dataframe in pickled file
    new_df.to_pickle(df_name+'.pkl')

#zone_lat = [30,45]
#zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]] # 0 : whole region, 1 : zone 1, 2 : zone 2, etc.
#zone_2[(zone_2['time'] > '2013-01-01')] 


#-------------------------------------------------------------------------------------------------------------------------------------------
# ) 

def create_nwa_tracks():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:

        
    Output:
    * does not return anything but creates nwa_tracks, an np array of track ids for nwa eddies
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # create a mask for eddies that fall within the Northwest Atlantic region
    nwa_center = np.where((278<=f.longitude) & (f.longitude<=312) & (24<=f.latitude) & (f.latitude<=53))

    # indices of eddy #eddy_4 = np.where(f.track==4)4 (in the northwest atlantic)


    # create lat, lon, time coordinates for the DataSet
    # this may not be necessary if I just use np.where to mask out eddies of interest
    f_new = f.assign_coords(dict( lon=((((f.longitude + 180) % 360) - 180)), lat= f.latitude,time=f.time ))
    
    nwa_tracks=[]

    for i in range(len(nwa_center[0])):
        track = int(f.track[int(nwa_center[0][i])].values)
        if track not in nwa_tracks:
            nwa_tracks.append(track)
        
    np.save('nwa_tracks.npy', nwa_tracks) 

#-------------------------------------------------------------------------------------------------------------------------------------------
# ) returns overall mean position of the Gulf Stream 1993 - 2018

def get_gs():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
    
    Output:
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # load GS file
    gs = loadmat('/Users/elenaperez/Desktop/chatts/data/GS.mat')
    
    # save Gulf Stream data as xarray
    return (gs['lon_cell'][0][0][0] - 360),(gs['lat_cell'][0][0][0])

#-------------------------------------------------------------------------------------------------------------------------------------------
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
    gs = loadmat('/Users/elenaperez/Desktop/chatts/data/GS.mat')
    
    return (gs['lon_cell_yearly'][year-1993][0][0] - 360),(gs['lat_cell_yearly'][year-1993][0][0])

#-------------------------------------------------------------------------------------------------------------------------------------------
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
    gs = loadmat('/Users/elenaperez/Desktop/chatts/data/GS.mat')
    
    return (gs['lon_cell_monthly'][(((year-1993)*12)+month)][0][0] - 360),(gs['lat_cell_monthly'][(((year-1993)*12)+month)][0][0])


#-------------------------------------------------------------------------------------------------------------------------------------------

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
    gs = loadmat('/Users/elenaperez/Desktop/chatts/data/gs/GS_daily_CMEMS_047_50cm_contours_1993_to_nrt.mat')
    
    'convert time array to ordinal dates'
    for d in range(len(gs['time'][0])-1):
        gs['time'][0][d] = gs['time'][0][d]+date.toordinal(date(1950,1,1))
    
    'get index in lon/lat arrays from time array'
    index_from_date = np.where(gs['time'][0]==date.toordinal(date(year,month,day)))[0][0]
    
    return (gs['lon_cell'][index_from_date][0][0] - 360),(gs['lat_cell'][index_from_date][0][0])


#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns lon, lat for a single eddy

def get_eddy_formation_loc(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns lon/lat location of eddy formation as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    return eddy[eddy['observation_number']==0]['longitude'],eddy[eddy['observation_number']==0]['latitude']

#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns location (lon, lat) of demise for a single eddy

def get_eddy_demise_loc(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns lon/lat location of eddy demise as tuple

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    return eddy.iloc[-1].longitude,eddy.iloc[-1].latitude

#-------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns lifespan of eddy

def get_eddy_lifespan(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy (DataFrame)  : Pandas DataFrame with data from a single eddy 

    Output:
        (Tuple)           : returns integer equal to lifespan of eddy

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    return eddy['observation_number'].iloc[-1]

#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns true if given eddy is warm core ring, false otherwise

def is_wcr(eddy):
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

#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns true if given eddy is cold core ring, false otherwise

def is_ccr(eddy):
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


#-------------------------------------------------------------------------------------------------------------------------------------------


# ) returns true if given eddy is WCR-like feature but propagate eastward

def is_eastward_wcr(eddy, gs):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function uses mutliple criterion to filter out warm core rings from mesoscale eddies. The criteria is:
        1) is the eddy anti-cyclonic
        2) is the eddy formed north of -0.25 deg of the Gulf path
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

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day,gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day,gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon
    
    eastward_propagation = eddy_moves_east(eddy) # True if eddy moves west

    return ((eddy_lat >= gs_lat_np.T[min_index]-0.25) & (is_geq_100m_isobath(eddy) & (eddy['cyclonic_type']==1).all()) & (eddy_lat <= gs_lat_np.T[min_index]+3) & (eastward_propagation)) # if true then eddy formation is north of gs path

#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns true if given eddy is CCR-like feature but propagate eastward

def is_eastward_ccr(eddy, gs):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function uses mutliple criterion to filter out cold core rings from mesoscale eddies. The criteria is:
        1) is the eddy cyclonic
        2) is the eddy formed south of +0.25 deg of the Gulf path
        3) is the eddy formed north of 34N
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

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day,gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day,gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon
    
    eastward_propagation = eddy_moves_east(eddy) # True if eddy moves west
    
    return ((eddy_lat <= gs_lat_np.T[min_index]+0.25) & (is_geq_100m_isobath(eddy)) & (eddy_lat >= gs_lat_np.T[min_index]-3) & ((eddy['cyclonic_type']==-1).all()) & (eddy_lat >= 34) & (eastward_propagation))


#-------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------
# ) returns true if eddy moves westward

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


#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns true if eddy moves eastward

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


#-------------------------------------------------------------------------------------------------------------------------------------------

# ) returns the closest gs lon, lat to a given eddy

def closest_gs_lat(eddy,gs):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a single eddy and all the GS positions and calculates for that given eddy what is 
    the closest position of the GS.
    
    Input:
        eddy (DataFrame) : Pandas DataFrame with data from a single eddy 
        gs (array)       : array of GS positions

    Output:
        (bool)           : returns true if eddy propagates westwards, false otherwise

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    eddy_lon = get_eddy_formation_loc(eddy)[0] # lon of eddy formation
    eddy_lat = get_eddy_formation_loc(eddy)[1] # lat of eddy formation

    eddy_year = get_eddy_formation_time(eddy)[0] # year of eddy formation
    eddy_month = get_eddy_formation_time(eddy)[1] # month of eddy formation
    eddy_day = get_eddy_formation_time(eddy)[2] # day of eddy formation

    gs_lon_np = get_gs_day(eddy_year,eddy_month,eddy_day,gs)[0] # np array of gs lon
    gs_lat_np = get_gs_day(eddy_year,eddy_month,eddy_day,gs)[1] # np array of gs lat

    gs_lon_len = len(gs_lon_np) # length of gs lon array

    eddy_lon_np = np.full(gs_lon_len, eddy_lon) # populate 1D array of gs_lon len with eddy formation lon

    min_index = np.argmin(abs(gs_lon_np-eddy_lon_np)) # what index of gs lon array is the closest to eddy formation lon
    
    return gs_lon_np[min_index], gs_lat_np[min_index]


#-------------------------------------------------------------------------------------------------------------------------------------------


# ) returns true if eddy depth is greater than or equal to 100-m isobath

def is_geq_100m_isobath(eddy):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function determins if an eddy is *not* on the shelf. If its depth is greater than or equal to
    the 100-m isobath, return True since it is *not* on the shelf. If its depth is less than 100-m,
    return False because it is on the shelf.
    
    Input:
        eddy (DataFrame) : Pandas DataFrame with data from a single eddy 

    Output:
        (bool)           : returns true if eddy is deeper than 100m isobath (on shelf), returns false otherwise

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    bathy = xr.open_dataset('/Users/elenaperez/Desktop/chatts/data/nwa_bathy.nc')
    
    target_bathy = 500 # desired isobath

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

#-------------------------------------------------------------------------------------------------------------------------------------------
# ) This function loads an unedited Chelton tracks file, processes, and saves a pickled DataFrame of the 
#    eddy tracks within a specified region (x_bnds, y_bnds) and specified time period (1993-2017 otherwise)

def make_eddy_region_df(time_period, x_bnds, y_bnds, df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        time_period (array) : time period, e.g. ['1994-01-01','1994-12-31']. False if no time period
        x_bnds (array)      : western and eastern longitude bounds of box in 0-360 deg, e.g. [278,312] 
        y_bnds (array)      : southern and northern latitude bounds of box, e.g. [30,45] N
        df_name (string)    : name of file to be saved 

        
    Output:
        eddy_df (DataFrame)  : pandas dataframe of eddies contained in the specified box and time period
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # load file as xarray DataArray
    path_to_chelton_tracks = '/Users/elenaperez/Desktop/chatts/data/eddy_trajectory_dt_2.0_19930101_20200307.nc'
    eddy_da = xr.open_dataset(path_to_chelton_tracks)
    'eddy_da (DataArray)     : Xarray DataArray of global mesoscale eddies '

    # create numpy array of indices within bounding box
    region_indices = np.where((x_bnds[0]<=eddy_da.longitude) & (eddy_da.longitude<=x_bnds[1]) & (y_bnds[0]<=eddy_da.latitude) & (eddy_da.latitude<=y_bnds[1]))
    'region_indices (array)  : indices of eddies in the larger region, e.g. [274, 275, 276, ...]'

    # create numpy array of track ids within the bounding box
    for i in range(len(region_indices[0])):
            track = int(eddy_da.track[int(region_indices[0][i])].values)
            if track not in region_tracks:
                region_tracks.append(track)
    'region_tracks (array) : track ids of eddies in the larger region, e.g. [4, 5, 7, ...]'

    # save the regional track ids as numpy array
    region_tracks_filename = 'nwa_tracks'
    np.save(region_tracks_filename+'.npy', region_tracks) 

    # convert eddy_da (DataArray) to Pandas DataFrame 
    eddy_df = eddy_da.to_dataframe()
    'eddy_df (DataFrame) : pandas dataframe of eddies in specified region, e.g. northwest atlantic'

    # filter out eddy tracks that are not in specified region, e.g. Northwest Atlantic
    eddy_df = eddy_df[eddy_df['track'].isin(region_tracks)]

    # convert lon from 0-360 to -180-180
    eddy_df['longitude'] = ((((eddy_df['longitude'] + 180) % 360) - 180)) 

    # save regional eddy DataFrame in pickled file
    eddy_df.to_pickle(df_name+'.pkl')

    # make_eddy_region_df(False, [278,312], [24,54], 'nwa_eddy_df')


#-------------------------------------------------------------------------------------------------------------------------------------------
# ) This function takes a regional eddy DataFrame and cuts it down to Gangopadhyay et al., 2019's zones 1-4

def make_eddy_zone_df(eddy_df, time_period, which_zone, zone_df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy_df (DataFrame)    : pandas dataframe of eddies in specified region, e.g. Northwest Atlantic
        time_period (array)    : time period, e.g. ['1994-01-01','1994-12-31']. False if no time period
        which_zone             : zones based on Gangopadhyay et al., 2020. 0 = all zones, 1 = zone 1, etc.
        zone_df_name (string)  : name of file to be saved 

        
    Output:
        new_df (DataFrame)     : pandas dataframe of eddies contained in the specified box and time period
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    zone_lat = [30,45]
    zone_lon = [[-75,-55],[-75,-70],[-70,-65],[-65,-60],[-60,-55]]

    # loop through eddy tracks to see if eddy is formed in the region
    zone_tracks = []
    for i in np.array(eddy_df['track'].unique()):
        if ((eddy_df[eddy_df['track']==i].iloc[0].longitude >= zone_lon[which_zone][0]) & (eddy_df[eddy_df['track']==i].iloc[0].longitude <= zone_lon[which_zone][1]) & (eddy_df[eddy_df['track']==i].iloc[0].latitude>= zone_lat[0]) & (eddy_df[eddy_df['track']==i].iloc[0].latitude<= zone_lat[1])):
            zone_tracks.append(i)

    # cut to specified bounding box
    eddy_zone_df = eddy_df[eddy_df['track'].isin(zone_tracks)] 

    # cut to specified time period
    if time_period!=False:
        eddy_zone_df = eddy_zone_df[(eddy_zone_df['time'] >= time_period[0]) & (eddy_zone_df['time'] <= time_period[0])] 

    # save regional eddy DataFrame in pickled file
    eddy_zone_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+zone_df_name+'.pkl')

# # example of how to use this function
# make_eddy_zone_df(eddy_ccr_df, False, 0, 'zone_ccrs')
# make_eddy_zone_df(eddy_ccr_df, False, 1, 'zone1_ccrs')
# make_eddy_zone_df(eddy_ccr_df, False, 2, 'zone2_ccrs')
# make_eddy_zone_df(eddy_ccr_df, False, 3, 'zone3_ccrs')
# make_eddy_zone_df(eddy_ccr_df, False, 4, 'zone4_ccrs')

# make_eddy_zone_df(eddy_wcr_df, False, 0, 'zone_wcrs')
# make_eddy_zone_df(eddy_wcr_df, False, 1, 'zone1_wcrs')
# make_eddy_zone_df(eddy_wcr_df, False, 2, 'zone2_wcrs')
# make_eddy_zone_df(eddy_wcr_df, False, 3, 'zone3_wcrs')
# make_eddy_zone_df(eddy_wcr_df, False, 4, 'zone4_wcrs')

#-------------------------------------------------------------------------------------------------------------------------------------------


# ) This function takes an eddy dataFrame and trims it down to eddies that pass thru the CH box

def make_eddy_ch_df(eddy_df, time_period, new_df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy_df (DataFrame)    : pandas dataframe of eddies in specified region, e.g. Northwest Atlantic
        time_period (array)    : time period, e.g. ['1994-01-01','1994-12-31']. False if no time period
        new_df_name (string)   : name of file to be saved 
        
    Output:
        new_df (DataFrame)     : pandas dataframe of eddies contained in the specified box and time period
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # lat,lon of "the point" and CH box
    the_point = [35.14111, -74.8] 
    ch_box = [(the_point[0]-1, the_point[0]+1),(the_point[1]-1, the_point[1]+1)] 

    ch_indices = np.where((eddy_df['longitude'] >= ch_box[1][0])  & (eddy_df['longitude'] <= ch_box[1][1]) & (eddy_df['latitude'] <= ch_box[0][1]) & (eddy_df['latitude'] >= ch_box[0][0]))

    # get unique track array for CH eddies
    ch_tracks = eddy_df.iloc[ch_indices]['track'].unique()

    # cut to specified bounding box
    eddy_ch_df = eddy_df[eddy_df['track'].isin(ch_tracks)] 

    # cut to specified time period
    if time_period!=False:
        eddy_ch_df = eddy_ch_df[(eddy_ch_df['time'] >= time_period[0]) & (eddy_ch_df['time'] <= time_period[0])] 

    # save regional eddy DataFrame in pickled file
    eddy_ch_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+new_df_name+'.pkl')

# # example of how to use this function
# make_eddy_ch_df(nwa_eddies, False, 'ch_eddies')

#-------------------------------------------------------------------------------------------------------------------------------------------


# ) takes an eddy DataFrame and saves a new DataFrame of just WCRs or CCRs
def eddy_df_to_ring_df(eddy_df):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy_df (DataFrame)    : pandas dataframe of eddies
        df_name (String)       : name of df that will be saved

        
    Output:
        new_df (DataFrame)     : pandas dataframe of specified type of rings (e.g. WCRs or CCRs)
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    ### WESTWARD PROPAGATING FEATURES ###
    
    ## CCRs ##
    ccr_tracks = []
    for i in np.array(eddy_df['track'].unique()):
        eddy = eddy_df[eddy_df['track']==i]
        if ((eddy['time'].dt.year<2018).all()): # cuts off at 2018, since GS stop halfway through 2018
            if is_ccr(eddy).all():
                ccr_tracks.append(i)

    # cut eddies that aren't ccrs out of eddy df to form ring df
    eddy_ccr_df = eddy_df[eddy_df['track'].isin(ccr_tracks)]

    # save ccr df as pickled file
    eddy_ccr_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/nwa_ccr_day_df.pkl') 
        
    ## WCRs ##
    wcr_tracks = []
    for i in np.array(eddy_df['track'].unique()):
        eddy = eddy_df[eddy_df['track']==i]
        if ((eddy['time'].dt.year<2018).all()): # cuts off at 2018, since GS stop halfway through 2018
            if is_wcr(eddy).all():
                wcr_tracks.append(i)

    # cut eddies that aren't ccrs out of eddy df to form ring df
    eddy_wcr_df = eddy_df[eddy_df['track'].isin(wcr_tracks)]

    # save ccr df as pickled file
    eddy_wcr_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/nwa_wcr_day_df.pkl')
    
    
    
    ### EASTWARD PROPAGATING FEATURES ###
    
    # CCR-like features that propagate eastward
    east_ccr_tracks = []
    for i in np.array(eddy_df['track'].unique()):
        eddy = eddy_df[eddy_df['track']==i]
        if ((eddy['time'].dt.year<2018).all()): # cuts off at 2018, since GS stop halfway through 2018
            if is_eastward_ccr(eddy).all():
                east_ccr_tracks.append(i)

    # cut eddies that aren't ccrs out of eddy df to form ring df
    eddy_east_ccr_df = eddy_df[eddy_df['track'].isin(east_ccr_tracks)]

    # save ccr df as pickled file
    eddy_east_ccr_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/nwa_east_ccr_day_df.pkl') 
        
    ## WCRs ##
    east_wcr_tracks = []
    for i in np.array(eddy_df['track'].unique()):
        eddy = eddy_df[eddy_df['track']==i]
        if ((eddy['time'].dt.year<2018).all()): # cuts off at 2018, since GS stop halfway through 2018
            if is_eastward_wcr(eddy).all():
                east_wcr_tracks.append(i)

    # cut eddies that aren't ccrs out of eddy df to form ring df
    eddy_east_wcr_df = eddy_df[eddy_df['track'].isin(east_wcr_tracks)]

    # save ccr df as pickled file
    eddy_east_wcr_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/nwa_east_wcr_day_df.pkl')
    
    
    
    # # example of how to call function
    # eddy_df_to_ring_df(nwa_eddies)

#-------------------------------------------------------------------------------------------------------------------------------------------
# ) saves DataFrame of number of monthly eddy formations for an eddy DataFrame

def count_monthly_ring_formations(ring_df, ring_type, df_name):
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

    if ring_type=='wcr':
        var = ['month','wcr_formations','mean_wcr_formations']

    elif ring_type=='ccr':
        var = ['month','ccr_formations','mean_ccr_formations']   

    df_structure = np.zeros((len(month_range), len(var)))
    ring_month_count_df = pd.DataFrame(df_structure, columns = var)

    counter=0
    for month in month_range:
        monthly_avg = len((ring_df[(ring_df['time'].dt.month == month) & (ring_df['observation_number']==0)])['track'].unique())/((year_range[-1]-year_range[0])+1)
        monthly_formations = len((ring_df[(ring_df['time'].dt.month == month) & (ring_df['observation_number']==0)])['track'].unique())
        ring_month_count_df.iloc[counter]=[month, monthly_formations, monthly_avg]
        counter += 1
        
    ring_month_count_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+df_name+'.pkl')

   # example call
   # zone_wcr_month_count_df = count_monthly_ring_formations(zone_wcrs, 'wcr', 'zone_wcr_monthly_formations')

    

#-------------------------------------------------------------------------------------------------------------------------------------------
# ) saves DataFrame of number of annual ring formations for an eddy DataFrame

def count_annual_ring_formations(ring_df, ring_type, df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of annual formations for 
    that type of rings (e.g. wcr or ccr)
    
    Input:
        ring_df (DataFrame)              : pandas dataframe of rings
        ring_type (String)               : 'wcr' for warm core ring, 'ccr' for cold core ring
        
    Output:
        ring_annual_count_df (DataFrame) : pandas dataframe of annual ring formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    year_range = np.arange(min(ring_df['time'].dt.year), max(ring_df['time'].dt.year)+1)

    if ring_type=='wcr':
        var = ['year','wcr_formations']

    elif ring_type=='ccr':
        var = ['year','ccr_formations']    

    df_structure = np.zeros((len(year_range), len(var)))
    ring_annual_count_df = pd.DataFrame(df_structure, columns = var)

    counter = 0
    for i in year_range:
        annual_formations = len((ring_df[(ring_df['time'].dt.year == i) & (ring_df['observation_number']==0)])['track'].unique())
        ring_annual_count_df.iloc[counter]=[i, annual_formations]
        counter += 1
        
    ring_annual_count_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+df_name+'.pkl')
       
## example of how to call this function
# count_annual_ring_formations(eddy_ccr_df,'ccr','nwa_ccr_annual_formations')
# count_annual_ring_formations(zone_ccrs,'ccr','zone_ccr_annual_formations')


#-------------------------------------------------------------------------------------------------------------------------------------------
# ) saves DataFrame of all years, month ring counts

def count_all_ring_formations(ring_df, ring_type, df_name):
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

    if ring_type=='wcr':
        var = ['year','month','wcr_formations']

    elif ring_type=='ccr':
        var = ['year','month','ccr_formations']   

    df_structure = np.zeros((len(month_range)*len(year_range), len(var)))
    ring_all_count_df = pd.DataFrame(df_structure, columns = var)

    counter=0
    for year in year_range:
        for month in month_range:
            ring_formations = len((ring_df[(ring_df['time'].dt.month == month) & (ring_df['time'].dt.year == year) & (ring_df['observation_number']==0)])['track'].unique())
            ring_all_count_df.iloc[counter]=[year, month, ring_formations]
            counter += 1

    ring_all_count_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+df_name+'.pkl')

    
#-------------------------------------------------------------------------------------------------------------------------------------------
# ) saves merged DataFrame of ring monthly formations by zone

def merge_monthly_ring_counts():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        
        
    Output:
        * no output, but saves the merged DataFrame
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # # open ring *MONTHLY* formation counts
    zone_ccr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ccr_monthly_formations.pkl')
    zone1_ccr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_ccr_monthly_formations.pkl')
    zone2_ccr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_ccr_monthly_formations.pkl')
    zone3_ccr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_ccr_monthly_formations.pkl')
    zone4_ccr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_ccr_monthly_formations.pkl')

    zone_wcr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_wcr_monthly_formations.pkl')
    zone1_wcr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_wcr_monthly_formations.pkl')
    zone2_wcr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_wcr_monthly_formations.pkl')
    zone3_wcr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_wcr_monthly_formations.pkl')
    zone4_wcr_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_wcr_monthly_formations.pkl')
    
    ## WCRs ##
    # merge each dataframe onto the all_zones dataframe
    zone_wcr_monthly_formations.columns = ['month', 'all_zones','all_zones_mean'] # all zones
    zone_wcr_monthly_formations['zone_1'] = zone1_wcr_monthly_formations['wcr_formations'] # zone 1
    zone_wcr_monthly_formations['zone_1_mean'] = zone1_wcr_monthly_formations['mean_wcr_formations']
    zone_wcr_monthly_formations['zone_2'] = zone2_wcr_monthly_formations['wcr_formations'] # zone 2
    zone_wcr_monthly_formations['zone_2_mean'] = zone2_wcr_monthly_formations['mean_wcr_formations']
    zone_wcr_monthly_formations['zone_3'] = zone3_wcr_monthly_formations['wcr_formations'] # zone 3
    zone_wcr_monthly_formations['zone_3_mean'] = zone3_wcr_monthly_formations['mean_wcr_formations']
    zone_wcr_monthly_formations['zone_4'] = zone4_wcr_monthly_formations['wcr_formations'] # zone 4
    zone_wcr_monthly_formations['zone_4_mean'] = zone4_wcr_monthly_formations['mean_wcr_formations']

    # save dataframe
    zone_wcr_monthly_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_wcr_monthly_formations.pkl') 
    
    ## CCRs ##
    # merge each dataframe onto the all_zones dataframe
    zone_ccr_monthly_formations.columns = ['month', 'all_zones','all_zones_mean'] # all zones
    zone_ccr_monthly_formations['zone_1'] = zone1_ccr_monthly_formations['ccr_formations'] # zone 1
    zone_ccr_monthly_formations['zone_1_mean'] = zone1_ccr_monthly_formations['mean_ccr_formations']
    zone_ccr_monthly_formations['zone_2'] = zone2_ccr_monthly_formations['ccr_formations'] # zone 2
    zone_ccr_monthly_formations['zone_2_mean'] = zone2_ccr_monthly_formations['mean_ccr_formations']
    zone_ccr_monthly_formations['zone_3'] = zone3_ccr_monthly_formations['ccr_formations'] # zone 3
    zone_ccr_monthly_formations['zone_3_mean'] = zone3_ccr_monthly_formations['mean_ccr_formations']
    zone_ccr_monthly_formations['zone_4'] = zone4_ccr_monthly_formations['ccr_formations'] # zone 4
    zone_ccr_monthly_formations['zone_4_mean'] = zone4_ccr_monthly_formations['mean_ccr_formations']

    # save dataframe
    zone_ccr_monthly_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ccr_monthly_formations.pkl')
 

 # example call
# merge_monthly_ring_counts()


#-------------------------------------------------------------------------------------------------------------------------------------------
# ) returns merged DataFrame of ring formations by zone

def merge_annual_ring_counts():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        
        
    Output:
        * no output, but saves the merged DataFrame
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    ## OPEN
    zone_ccr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ccr_annual_formations.pkl')
    zone1_ccr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_ccr_annual_formations.pkl')
    zone2_ccr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_ccr_annual_formations.pkl')
    zone3_ccr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_ccr_annual_formations.pkl')
    zone4_ccr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_ccr_annual_formations.pkl')

    zone_wcr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_wcr_annual_formations.pkl')
    zone1_wcr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_wcr_annual_formations.pkl')
    zone2_wcr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_wcr_annual_formations.pkl')
    zone3_wcr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_wcr_annual_formations.pkl')
    zone4_wcr_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_wcr_annual_formations.pkl')
    
    ## WCRs ##
    # merge each dataframe onto the all_zones dataframe
    zone_wcr_annual_formations.columns = ['year', 'all_zones']
    zone_wcr_annual_formations['zone_1'] = zone1_wcr_annual_formations['wcr_formations']
    zone_wcr_annual_formations['zone_2'] = zone2_wcr_annual_formations['wcr_formations']
    zone_wcr_annual_formations['zone_3'] = zone3_wcr_annual_formations['wcr_formations']
    zone_wcr_annual_formations['zone_4'] = zone4_wcr_annual_formations['wcr_formations']

    # save dataframe
    zone_wcr_annual_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_wcr_annual_formations.pkl') 
    
    ## CCRs ##
    # merge each dataframe onto the all_zones dataframe
    zone_ccr_annual_formations.columns = ['year', 'all_zones']
    zone_ccr_annual_formations['zone_1'] = zone1_ccr_annual_formations['ccr_formations']
    zone_ccr_annual_formations['zone_2'] = zone2_ccr_annual_formations['ccr_formations']
    zone_ccr_annual_formations['zone_3'] = zone3_ccr_annual_formations['ccr_formations']
    zone_ccr_annual_formations['zone_4'] = zone4_ccr_annual_formations['ccr_formations']

    # save dataframe
    zone_ccr_annual_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ccr_annual_formations.pkl') 
    
## example of how to call the function
# merge_annual_ring_counts() 


#-------------------------------------------------------------------------------------------------------------------------------------------
# ) saves merged DataFrame of ring all years, all months formations by zone

def merge_all_ring_counts():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        
        
    Output:
        * no output, but saves the merged DataFrame
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # # open ring *MONTHLY* formation counts
    zone_ccr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ccr_all_formations.pkl')
    zone1_ccr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_ccr_all_formations.pkl')
    zone2_ccr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_ccr_all_formations.pkl')
    zone3_ccr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_ccr_all_formations.pkl')
    zone4_ccr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_ccr_all_formations.pkl')

    zone_wcr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_wcr_all_formations.pkl')
    zone1_wcr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_wcr_all_formations.pkl')
    zone2_wcr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_wcr_all_formations.pkl')
    zone3_wcr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_wcr_all_formations.pkl')
    zone4_wcr_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_wcr_all_formations.pkl')
    
    ## WCRs ##
    # merge each dataframe onto the all_zones dataframe
    zone_wcr_all_formations.columns = ['year', 'month', 'all_zones'] # all zones
    zone_wcr_all_formations['zone_1'] = zone1_wcr_all_formations['wcr_formations']
    zone_wcr_all_formations['zone_2'] = zone2_wcr_all_formations['wcr_formations'] 
    zone_wcr_all_formations['zone_3'] = zone3_wcr_all_formations['wcr_formations'] 
    zone_wcr_all_formations['zone_4'] = zone4_wcr_all_formations['wcr_formations']

    # save dataframe
    zone_wcr_all_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_wcr_all_formations.pkl') 
    
    ## CCRs ##
    # merge each dataframe onto the all_zones dataframe
    zone_ccr_all_formations.columns = ['year','month', 'all_zones'] # all zones
    zone_ccr_all_formations['zone_1'] = zone1_ccr_all_formations['ccr_formations']
    zone_ccr_all_formations['zone_2'] = zone2_ccr_all_formations['ccr_formations'] 
    zone_ccr_all_formations['zone_3'] = zone3_ccr_all_formations['ccr_formations'] 
    zone_ccr_all_formations['zone_4'] = zone4_ccr_all_formations['ccr_formations']

    # save dataframe
    zone_ccr_all_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ccr_all_formations.pkl')
 

 # example call
# merge_all_ring_counts()
#------------------------------------------------------------------------------------------------------------------------------# ) saves DataFrame of number of monthly eddy formations for an eddy DataFrame

def count_monthly_eddy_formations(eddy_df, eddy_type, df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of eddies and return a new DataFrame of the number of monthly formations 
    for that type of eddy (e.g. anticyclone or cyclone)
    
    Input:
        eddy_df (DataFrame)              : pandas dataframe of rings
        eddy_type (String)               : 'anticylonic' +1, 'cyclonic' -1
        
    Output:
        eddy_month_count_df (DataFrame)  : pandas dataframe of monthly eddy formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    year_range = np.arange(min(eddy_df['time'].dt.year), max(eddy_df['time'].dt.year)+1)
    month_range = np.arange(1,13)

    if eddy_type=='anticyclonic':
        var = ['month','anticyclonic_formations','mean_anticyclonic_formations']

    elif eddy_type=='cyclonic':
        var = ['month','cyclonic_formations','mean_cyclonic_formations']   

    df_structure = np.zeros((len(month_range), len(var)))
    eddy_month_count_df = pd.DataFrame(df_structure, columns = var)

    counter=0
    for month in month_range:
        monthly_avg = len((eddy_df[(eddy_df['time'].dt.month == month) & (eddy_df['observation_number']==0)])['track'].unique())/((year_range[-1]-year_range[0])+1)
        monthly_formations = len((eddy_df[(eddy_df['time'].dt.month == month) & (eddy_df['observation_number']==0)])['track'].unique())
        eddy_month_count_df.iloc[counter]=[month, monthly_formations, monthly_avg]
        counter += 1
        
    eddy_month_count_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+df_name+'.pkl')

   ## example call
   # count_monthly_eddy_formations(zone_aeddies, 'anticyclonic', 'zone_aeddy_monthly_formations')   

    
    
#------------------------------------------------------------------------------------------------------------------------------# ) saves DataFrame of number of annual ring formations for an eddy DataFrame

def count_annual_eddy_formations(eddy_df, eddy_type, df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of annual formations for 
    that type of rings (e.g. wcr or ccr)
    
    Input:
        eddy_df (DataFrame)              : pandas dataframe of rings
        eddy_type (String)               : 'wcr' for warm core ring, 'ccr' for cold core ring
        
    Output:
        eddy_annual_count_df (DataFrame) : pandas dataframe of annual ring formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    year_range = np.arange(min(eddy_df['time'].dt.year), max(eddy_df['time'].dt.year)+1)

    if eddy_type=='anticyclonic':
        var = ['year','anticyclonic_formations']

    elif eddy_type=='cyclonic':
        var = ['year','cyclonic_formations']    

    df_structure = np.zeros((len(year_range), len(var)))
    eddy_annual_count_df = pd.DataFrame(df_structure, columns = var)

    counter = 0
    for i in year_range:
        annual_formations = len((eddy_df[(eddy_df['time'].dt.year == i) & (eddy_df['observation_number']==0)])['track'].unique())
        eddy_annual_count_df.iloc[counter]=[i, annual_formations]
        counter += 1
        
    eddy_annual_count_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+df_name+'.pkl')
       
## example of how to call this function
# count_annual_eddy_formations(zone_ceddies,'cyclonic','zone1_ceddy_annual_formations')


#------------------------------------------------------------------------------------------------------------------------------
# ) saves DataFrame of all years, month eddy counts

def count_all_eddy_formations(eddy_df, eddy_type, df_name):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    This function takes a DataFrame of rings and return a new DataFrame of the number of monthly formations 
    for that type of ring (e.g. wcr or ccr)
    
    Input:
        eddy_df (DataFrame)              : pandas dataframe of rings
        eddy_type (String)               : 'anticyclonic' (+1), 'cyclonic' (-1)
        
    Output:
        eddy_month_count_df (DataFrame)  : pandas dataframe of monthly, yearly eddy formations
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    year_range = np.arange(min(eddy_df['time'].dt.year), max(eddy_df['time'].dt.year)+1)
    month_range = np.arange(1,13)

    if eddy_type=='anticyclonic':
        var = ['year','month','anticyclonic_formations']

    elif eddy_type=='cyclonic':
        var = ['year','month','cyclonic_formations']   

    df_structure = np.zeros((len(month_range)*len(year_range), len(var)))
    eddy_all_count_df = pd.DataFrame(df_structure, columns = var)

    counter=0
    for year in year_range:
        for month in month_range:
            eddy_formations = len((eddy_df[(eddy_df['time'].dt.month == month) & (eddy_df['time'].dt.year == year) & (eddy_df['observation_number']==0)])['track'].unique())
            eddy_all_count_df.iloc[counter]=[year, month, eddy_formations]
            counter += 1

    eddy_all_count_df.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/'+df_name+'.pkl')


# # example call
# count_all_eddy_formations(zone1_aeddies, 'anticyclonic', 'zone1_aeddy_all_formations') 



#------------------------------------------------------------------------------------------------------------------------------# ) saves merged DataFrame of eddy monthly formations by zone

def merge_monthly_eddy_counts():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        
        
    Output:
        * no output, but saves the merged DataFrame
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # # open ring *MONTHLY* formation counts
    zone_ceddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ceddy_monthly_formations.pkl')
    zone1_ceddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_ceddy_monthly_formations.pkl')
    zone2_ceddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_ceddy_monthly_formations.pkl')
    zone3_ceddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_ceddy_monthly_formations.pkl')
    zone4_ceddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_ceddy_monthly_formations.pkl')

    zone_aeddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_aeddy_monthly_formations.pkl')
    zone1_aeddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_aeddy_monthly_formations.pkl')
    zone2_aeddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_aeddy_monthly_formations.pkl')
    zone3_aeddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_aeddy_monthly_formations.pkl')
    zone4_aeddy_monthly_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_aeddy_monthly_formations.pkl')

    ## ANTICYCLONIC ##
    # merge each dataframe onto the all_zones dataframe
    zone_aeddy_monthly_formations.columns = ['month', 'all_zones','all_zones_mean'] # all zones
    zone_aeddy_monthly_formations['zone_1'] = zone1_aeddy_monthly_formations['anticyclonic_formations'] # zone 1
    zone_aeddy_monthly_formations['zone_1_mean'] = zone1_aeddy_monthly_formations['mean_anticyclonic_formations']
    zone_aeddy_monthly_formations['zone_2'] = zone2_aeddy_monthly_formations['anticyclonic_formations'] # zone 2
    zone_aeddy_monthly_formations['zone_2_mean'] = zone2_aeddy_monthly_formations['mean_anticyclonic_formations']
    zone_aeddy_monthly_formations['zone_3'] = zone3_aeddy_monthly_formations['anticyclonic_formations'] # zone 3
    zone_aeddy_monthly_formations['zone_3_mean'] = zone3_aeddy_monthly_formations['mean_anticyclonic_formations']
    zone_aeddy_monthly_formations['zone_4'] = zone4_aeddy_monthly_formations['anticyclonic_formations'] # zone 4
    zone_aeddy_monthly_formations['zone_4_mean'] = zone4_aeddy_monthly_formations['mean_anticyclonic_formations']

    # save dataframe
    zone_aeddy_monthly_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_aeddy_monthly_formations.pkl') 
    
    ## CYCLONIC ##
    # merge each dataframe onto the all_zones dataframe
    zone_ceddy_monthly_formations.columns = ['month', 'all_zones','all_zones_mean'] # all zones
    zone_ceddy_monthly_formations['zone_1'] = zone1_ceddy_monthly_formations['cyclonic_formations'] # zone 1
    zone_ceddy_monthly_formations['zone_1_mean'] = zone1_ceddy_monthly_formations['mean_cyclonic_formations']
    zone_ceddy_monthly_formations['zone_2'] = zone2_ceddy_monthly_formations['cyclonic_formations'] # zone 2
    zone_ceddy_monthly_formations['zone_2_mean'] = zone2_ceddy_monthly_formations['mean_cyclonic_formations']
    zone_ceddy_monthly_formations['zone_3'] = zone3_ceddy_monthly_formations['cyclonic_formations'] # zone 3
    zone_ceddy_monthly_formations['zone_3_mean'] = zone3_ceddy_monthly_formations['mean_cyclonic_formations']
    zone_ceddy_monthly_formations['zone_4'] = zone4_ceddy_monthly_formations['cyclonic_formations'] # zone 4
    zone_ceddy_monthly_formations['zone_4_mean'] = zone4_ceddy_monthly_formations['mean_cyclonic_formations']

    # save dataframe
    zone_ceddy_monthly_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ceddy_monthly_formations.pkl')
 

    ## example call
    # merge_monthly_eddy_counts()


    
#------------------------------------------------------------------------------------------------------------------------------
# ) returns merged DataFrame of annual eddy formations by zone

def merge_annual_eddy_counts():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        
        
    Output:
        * no output, but saves the merged DataFrame
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    ## OPEN
    zone_ceddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ceddy_annual_formations.pkl')
    zone1_ceddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_ceddy_annual_formations.pkl')
    zone2_ceddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_ceddy_annual_formations.pkl')
    zone3_ceddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_ceddy_annual_formations.pkl')
    zone4_ceddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_ceddy_annual_formations.pkl')

    zone_aeddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_aeddy_annual_formations.pkl')
    zone1_aeddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_aeddy_annual_formations.pkl')
    zone2_aeddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_aeddy_annual_formations.pkl')
    zone3_aeddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_aeddy_annual_formations.pkl')
    zone4_aeddy_annual_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_aeddy_annual_formations.pkl')

    ## ANTI-CYCLONES ##
    # merge each dataframe onto the all_zones dataframe
    zone_aeddy_annual_formations.columns = ['year', 'all_zones']
    zone_aeddy_annual_formations['zone_1'] = zone1_aeddy_annual_formations['anticyclonic_formations']
    zone_aeddy_annual_formations['zone_2'] = zone2_aeddy_annual_formations['anticyclonic_formations']
    zone_aeddy_annual_formations['zone_3'] = zone3_aeddy_annual_formations['anticyclonic_formations']
    zone_aeddy_annual_formations['zone_4'] = zone4_aeddy_annual_formations['anticyclonic_formations']

    # save dataframe
    zone_aeddy_annual_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_aeddy_annual_formations.pkl') 
    
    ## CYCLONES ##
    # merge each dataframe onto the all_zones dataframe
    zone_ceddy_annual_formations.columns = ['year', 'all_zones']
    zone_ceddy_annual_formations['zone_1'] = zone1_ceddy_annual_formations['cyclonic_formations']
    zone_ceddy_annual_formations['zone_2'] = zone2_ceddy_annual_formations['cyclonic_formations']
    zone_ceddy_annual_formations['zone_3'] = zone3_ceddy_annual_formations['cyclonic_formations']
    zone_ceddy_annual_formations['zone_4'] = zone4_ceddy_annual_formations['cyclonic_formations']

    # save dataframe
    zone_ceddy_annual_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ceddy_annual_formations.pkl') 
    
## example of how to call the function
# merge_annual_eddy_counts()   



#------------------------------------------------------------------------------------------------------------------------------
# ) saves merged DataFrame of eddy all years, all months formations by zone

def merge_all_eddy_counts():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        
        
    Output:
        * no output, but saves the merged DataFrame
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # # open ring *MONTHLY* formation counts
    zone_ceddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ceddy_all_formations.pkl')
    zone1_ceddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_ceddy_all_formations.pkl')
    zone2_ceddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_ceddy_all_formations.pkl')
    zone3_ceddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_ceddy_all_formations.pkl')
    zone4_ceddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_ceddy_all_formations.pkl')

    zone_aeddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_aeddy_all_formations.pkl')
    zone1_aeddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone1_aeddy_all_formations.pkl')
    zone2_aeddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone2_aeddy_all_formations.pkl')
    zone3_aeddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone3_aeddy_all_formations.pkl')
    zone4_aeddy_all_formations = pd.read_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone4_aeddy_all_formations.pkl')
    
    ## ANTICYCLONES ##
    # merge each dataframe onto the all_zones dataframe
    zone_aeddy_all_formations.columns = ['year', 'month', 'all_zones'] # all zones
    zone_aeddy_all_formations['zone_1'] = zone1_aeddy_all_formations['anticyclonic_formations']
    zone_aeddy_all_formations['zone_2'] = zone2_aeddy_all_formations['anticyclonic_formations'] 
    zone_aeddy_all_formations['zone_3'] = zone3_aeddy_all_formations['anticyclonic_formations'] 
    zone_aeddy_all_formations['zone_4'] = zone4_aeddy_all_formations['anticyclonic_formations']

    # save dataframe
    zone_aeddy_all_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_aeddy_all_formations.pkl') 
    
    ## CYCLONES ##
    # merge each dataframe onto the all_zones dataframe
    zone_ceddy_all_formations.columns = ['year','month', 'all_zones'] # all zones
    zone_ceddy_all_formations['zone_1'] = zone1_ceddy_all_formations['cyclonic_formations']
    zone_ceddy_all_formations['zone_2'] = zone2_ceddy_all_formations['cyclonic_formations'] 
    zone_ceddy_all_formations['zone_3'] = zone3_ceddy_all_formations['cyclonic_formations'] 
    zone_ceddy_all_formations['zone_4'] = zone4_ceddy_all_formations['cyclonic_formations']

    # save dataframe
    zone_ceddy_all_formations.to_pickle('/Users/elenaperez/Desktop/chatts/data/pd_dataframes/zone_ceddy_all_formations.pkl')

 # example call
# merge_all_eddy_counts()



#-------------------------------------------------------------------------------------------------------------------------------------------

# ) 

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#-------------------------------------------------------------------------------------------------------------------------------------------

# ) 

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

