"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Utilities to plot META eddy trajectories & statistical analyses

    1) all_eddy_tracks_map          : returns a map of NWA with ring-like eddy tracks for the whole period
    2) spatial_formations_barplot   : returns a bar plot of formation of all ring-like eddy formations
    3) seasonality_lineplot         : returns a line plot of average monthly formations (seasonality)
    4)               :
    5)               :
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# import necessary packages and functions for plotting
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as dates
import matplotlib.colors as colors 
import cartopy
import calendar
import cartopy.crs as ccrs
from os import path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from cartopy.feature import NaturalEarthFeature
import seaborn as sns
import os,imageio
import datetime


# turn off warnings
import warnings
warnings.filterwarnings("ignore")


# adds upper level to working directory, this is where the utils folder is saved
import sys
sys.path.append("..")

# import the util data functions
from utils.ring_data_utils import *

#-------------------------------------------------------------------------------------------------------------------------------
# )
def all_eddy_tracks_map(eddy_ccr_df, eddy_wcr_df, bathy, title, fig_quality):
    # gangopadhyay census bounds
    x_bnds = [-80,-55] # lon, NWA: [-82,-48]
    y_bnds = [30,46] # lat, NWA: [24,53]
    
    # map projection
    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title, fontsize=20, y=0.925)

    ## WCR-like eddies ##
    for i in np.array(eddy_wcr_df['track'].unique()):
        eddy = eddy_wcr_df[eddy_wcr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='#EE4D4D') # anticyclonic
    
    ## CCR-like eddies ##     
    for i in np.array(eddy_ccr_df['track'].unique()):
        eddy = eddy_ccr_df[eddy_ccr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='#388FEF') # cyclonic
        
    # Gulf Stream mean path 1993-2022
    mean_gs = get_gs()
    ax.plot(mean_gs[0],mean_gs[1], zorder= 10000, linewidth = 6, color='k')
    
#     # Zones
    ax.plot([-75, -75], [20, 50], transform=proj, color = 'k', zorder= 10000, alpha = 0.5)
    ax.plot([-70, -70], [20, 50], transform=proj, color = 'k', zorder= 10000, alpha = 0.5)
    ax.plot([-65, -65], [20, 50], transform=proj, color = 'k', zorder= 10000, alpha = 0.5)
    ax.plot([-60, -60], [20, 50], transform=proj, color = 'k', zorder= 10000, alpha = 0.5)
    ax.plot([-55, -55], [20, 50], transform=proj, color = 'k', zorder= 10000, alpha = 0.5)
    ax.plot([-50, -50], [20, 50], transform=proj, color = 'k', zorder= 10000, alpha = 0.5)

#     ax.vlines(x = 3, ymin = 1, ymax = 3)

    # Cape Hatteras
    ax.plot(-75.5, 35.2,'*', markerfacecolor='#03DE08', markeredgewidth=1, markeredgecolor='k', markersize=25, transform=proj, zorder=100001)

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='k', alpha=0.8, zorder=10000, linestyles='solid') #,levels=[-4000,-1000,-100] '#5E5E5E'
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#EE4D4D', lw=10), Line2D([0], [0], color='#388FEF', lw=10)] 
    ax.legend(custom_lines, ['WCR-like = 320', 'CCR-like = 225'], loc='upper left', fontsize = 17).set_zorder(10001)
#     ax.legend(custom_lines, ['Anti-cyclonic = 4,800 ', 'Cyclonic = 4,438'], loc='upper left', fontsize = 15).set_zorder(10001)

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    # declare text size
    XTEXT_SIZE = 15
    YTEXT_SIZE = 15
    # to facilitate text rotation at bottom edge, ...
    # text justification: 'ha':'right' is used to avoid clashing with map's boundary
    # default of 'ha' is center, often causes trouble when text rotation is not zero
    gl.xlabel_style = {'size': XTEXT_SIZE, 'color': 'k', 'ha':'right'}
    gl.ylabel_style = {'size':YTEXT_SIZE, 'color': 'k', 'weight': 'normal'}
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# )
def spatial_formations_barplot(ccr_formations_df, wcr_formations_df, title, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones                   "
    "    wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones                   "
    "    title (Str)             : title of the figure, e.g. 'Map of the Northwest Atlantic'                "
    "    fig_quality (Int)       : integer of what dpi the image will be set to (e.g., 100 dpi)             "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a bar plot for distribution of formations by zone                                        "
    "                                                                                                       "
    "                                                                                                       "
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    # create figure
    fig,ax = plt.subplots();
    fig.suptitle(title+'('+str(int(min(ccr_formations_df['year'])))+' - '+str(int(max(ccr_formations_df['year'])))+')', y=0.95, fontsize=14);
    fig.set_dpi(fig_quality)

    # calculate optimal width of bar plots' bars
    N = 4 # number of bar pairs i want
    ind = np.arange(N)
    width = np.min(np.diff(ind))/3 # Calculate optimal width
    
    zones = ['Zone 1', 'Zone 2', 'Zone 3', 'Zone 4']

    # plotting
    ax.bar(ind, ccr_formations_df[['zone1','zone2','zone3','zone4']].sum() , width, label='CCR-like', color='#388FEF')
    ax.bar(ind + width, wcr_formations_df[['zone1','zone2','zone3','zone4']].sum(), width, label='WCR-like', color='#EE4D4D')

    ax.set_xticks([0.15,1.15,2.15,3.15])
    ax.set_xticklabels(zones)
    ax.set_ylabel('Number of Formations');
    ax.legend();
    
    return fig, ax


#-------------------------------------------------------------------------------------------------------------------------------
# )
def seasonality_lineplot(ccr_formations_df, wcr_formations_df, title, ylim_min, ylim_max, fig_quality):
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
        ccr_monthly_std.append((ccr_formations_df[ccr_formations_df['month']==i]['all_zones']).std())
        ccr_monthly_avg.append(ccr_formations_df[ccr_formations_df['month']==i]['all_zones'].mean())
        wcr_monthly_avg.append(wcr_formations_df[wcr_formations_df['month']==i]['all_zones'].mean())
        
    # lines for standard deviation
    ccr_minus_std = np.subtract(ccr_monthly_avg,ccr_monthly_std)
    wcr_minus_std = np.subtract(wcr_monthly_avg,wcr_monthly_std)
    ccr_plus_std = np.add(ccr_monthly_avg,ccr_monthly_std)
    wcr_plus_std = np.add(wcr_monthly_avg,wcr_monthly_std)

    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle(title+'('+str(int(min(ccr_formations_df['year'])))+' - '+str(int(max(ccr_formations_df['year'])))+')', y=0.935, fontsize=14);

    # plot
    ## CCR-like
    ax.plot(ccr_formations_df['month'].drop_duplicates(),ccr_monthly_avg,'-o',color='#388FEF');
#     plt.errorbar(ccr_formations_df['month'].drop_duplicates(), ccr_monthly_avg, ccr_plus_std/2, linestyle='None', marker='^', color='#388FEF')
    ax.fill_between(ccr_formations_df['month'].drop_duplicates(), ccr_plus_std, ccr_monthly_avg, alpha=0.15, color='#388FEF')
    ax.fill_between(ccr_formations_df['month'].drop_duplicates(), ccr_minus_std, ccr_monthly_avg, alpha=0.15, color='#388FEF')
    
    ## WCR-like
    ax.plot(wcr_formations_df['month'].drop_duplicates(),wcr_monthly_avg,'-o',color='#EE4D4D');
#     plt.errorbar(wcr_formations_df['month'].drop_duplicates(), wcr_monthly_avg, wcr_plus_std/2, linestyle='None', marker='^', color='#EE4D4D')
    ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_plus_std, wcr_monthly_avg, alpha=0.15, color='#EE4D4D')
    ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_minus_std, wcr_monthly_avg, alpha=0.15, color='#EE4D4D')
    
    zones = ['zone_1', 'zone_2', 'zone_3', 'zone_4']

    # axes formatting
    ax.set_xlabel('Months', fontweight='bold')
    ax.set_ylabel('Average No. of Formations',fontweight='bold');
    ax.set_ylim(ylim_min,ylim_max)
    ax.set_xlim(1,12)

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#EE4D4D', lw=4), Line2D([0], [0], color='#388FEF', lw=4)]
    ax.legend(custom_lines, ['WCR-like', 'CCR-like'], loc='upper left');
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------
# ) 
def timeseries_lineplot(ccr_formations_df, wcr_formations_df, title, fig_quality, ylim_min, ylim_max):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones                   "
    "    wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones                   "
    "    title (Str)                   : title of the figure, e.g. 'Map of the Northwest Atlantic'          "
    "    fig_quality (Int)             : integer of what dpi the image will be set to (e.g., 100 dpi)       "
    "    ylim_min (Int)                : integer of lowest yaxis tick mar                                   "
    "    ylim_max (Int)                : integer of highest yaxis tick mar                                   "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a line plot annual formations for whole period (time-series)                             "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # compute means for pre/post-2000 formations
    min_year_ccr = int(min(ccr_formations_df['year']))
    max_year_ccr = int(max(ccr_formations_df['year']))
    min_year_wcr = int(min(wcr_formations_df['year']))
    max_year_wcr = int(max(wcr_formations_df['year']))
    pre_2000_ccr = ccr_formations_df[ccr_formations_df['year']<2000] # mean rings count pre-2000
    post_2000_ccr = ccr_formations_df[ccr_formations_df['year']>=2000] # mean ring count post-2000
    pre_2000_wcr = wcr_formations_df[wcr_formations_df['year']<2000] # mean rings count pre-2000
    post_2000_wcr = wcr_formations_df[wcr_formations_df['year']>=2000] # mean ring count post-2000
    
    # create figure & axes
    fig,ax = plt.subplots() # figsize=(10,5), for comparison to Fig. 2 Silver et al., 2021
    fig.set_dpi(fig_quality)
    fig.suptitle(title+str(int(min(ccr_formations_df['year'])))+' - '+str(int(max(ccr_formations_df['year'])))+')', y=0.935, fontsize=11);

    # WCR-like eddies
    ax.plot(wcr_formations_df['year'],wcr_formations_df['all_zones'],'-o',color='#EE4D4D')
#     ax.plot([min_year_wcr,2000],[pre_2000_wcr['all_zones'].mean(),pre_2000_wcr['all_zones'].mean()], color='r', linestyle='-'); # pre_2000 average
#     ax.plot([2000,max_year_wcr],[post_2000_wcr['all_zones'].mean(),post_2000_wcr['all_zones'].mean()], color='r', linestyle='-'); # post_2000 average

    # CCR-like eddies
    ax.plot(ccr_formations_df['year'],ccr_formations_df['all_zones'],'-o',color='#388FEF')
#     ax.plot([min_year_ccr,2000],[pre_2000_ccr['all_zones'].mean(),pre_2000_ccr['all_zones'].mean()], color='b', linestyle='-'); # pre_2000 average
#     ax.plot([2000,max_year_ccr],[post_2000_ccr['all_zones'].mean(),post_2000_ccr['all_zones'].mean()], color='b', linestyle='-'); # post_2000 average

    # axes formatting
    ax.set_xlabel('Years',fontweight='bold')
    ax.set_ylabel('# of Formations',fontweight='bold');
    ax.set_ylim(ylim_min, ylim_max) # 1993-2017
#     ax.set_xlim(1980,2020) # to compare to fig. 2 Silver et al., 2021

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#EE4D4D', lw=4), Line2D([0], [0], color='#388FEF', lw=4)]
    ax.legend(custom_lines, ['WCR-like', 'CCR-like'], loc='upper left');
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------

def ringCensus_spatial_formations_barplot(wcr_formations_df, title, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones                   "
    "    wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones                   "
    "    title (Str)                   : title of the figure                                                "
    "    fig_quality (Int)             : integer of what dpi the image will be set to (e.g., 100 dpi)       "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a bar plot for distribution of formations by zone                                        "
    "                                                                                                       "
    "                                                                                                       "
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    # create figure
    fig,ax = plt.subplots();
    fig.suptitle(title+'('+str(int(min(wcr_formations_df['year'])))+' - '+str(int(max(wcr_formations_df['year'])))+')', y=0.95, fontsize=14);
    fig.set_dpi(fig_quality)

    # calculate optimal width of bar plots' bars
    N = 4 # number of bar pairs i want
    ind = np.arange(N)
    width = np.min(np.diff(ind))/3 # Calculate optimal width
    
    zones = ['Zone 1', 'Zone 2', 'Zone 3', 'Zone 4']

    # plotting
    ax.bar(ind + width, wcr_formations_df[['zone1','zone2','zone3','zone4']].sum(), width, label='WCR', color='#EE4D4D')

    ax.set_xticks([0.15,1.15,2.15,3.15])
    ax.set_xticklabels(zones)
    ax.set_ylabel('Number of Formations');
    ax.legend();
    
    return fig, ax

#-------------------------------------------------------------------------------------------------------------------------------

# plot GLED filtered and unfiltered
def gled_all_eddy_tracks_map(eddy30d_ccr_df, eddy30d_wcr_df, eddy90d_ccr_df, eddy90d_wcr_df, eddy180d_ccr_df, eddy180d_wcr_df, bathy, title, fig_quality):
    # gangopadhyay census bounds
    x_bnds = [-80,-55] # lon, NWA: [-82,-48]
    y_bnds = [30,46] # lat, NWA: [24,53]

    # map projection
    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title, fontsize=20, y=0.925)
    

    Files = {'gled30d_wcr': eddy30d_wcr_df, 'gled30d_ccr': eddy30d_ccr_df, 'gled90d_wcr': eddy90d_wcr_df, 'gled90d_ccr': eddy90d_ccr_df, 'gled180d_wcr': eddy180d_wcr_df, 'gled180d_ccr': eddy180d_ccr_df}

    for whichFile in Files:
        if 'wcr' in whichFile:
            eddy_wcr_df = Files[whichFile]
            ## WCR-like eddies ##
            for i in np.arange(len(eddy_wcr_df)):
                plt.plot(((np.array(eddy_wcr_df.iloc[i]['longitude'])+180)%360)-180,eddy_wcr_df.iloc[i]['latitude'],color='#EE4D4D')

        else: 
            eddy_ccr_df = Files[whichFile]
            # CCR-like eddies ##  
            for i in np.arange(len(eddy_ccr_df)):
                plt.plot(((np.array(eddy_ccr_df.iloc[i]['longitude'])+180)%360)-180,eddy_ccr_df.iloc[i]['latitude'],color='#388FEF')

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#EE4D4D', lw=4), Line2D([0], [0], color='#388FEF', lw=4)] 
    ax.legend(custom_lines, ['WCR-like', 'CCR-like'], loc='upper left')

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------


def ringCensus_timeseries_lineplot(wcr_formations_df, title, fig_quality, ylim_min, ylim_max):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones                   "
    "    wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones                   "
    "    title (Str)                   : title of the figure, e.g. 'Map of the Northwest Atlantic'          "
    "    fig_quality (Int)             : integer of what dpi the image will be set to (e.g., 100 dpi)       "
    "    ylim_min (Int)                : integer of lowest yaxis tick mar                                   "
    "    ylim_max (Int)                : integer of highest yaxis tick mar                                   "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a line plot annual formations for whole period (time-series)                             "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # compute means for pre/post-2000 formations
    min_year_wcr = int(min(wcr_formations_df['year']))
    max_year_wcr = int(max(wcr_formations_df['year']))
    pre_2000_wcr = wcr_formations_df[wcr_formations_df['year']<2000] # mean rings count pre-2000
    post_2000_wcr = wcr_formations_df[wcr_formations_df['year']>=2000] # mean ring count post-2000
    
    # create figure & axes
    fig,ax = plt.subplots() # figsize=(10,5), for comparison to Fig. 2 Silver et al., 2021
    fig.set_dpi(fig_quality)
    fig.suptitle(title+str(int(min(wcr_formations_df['year'])))+' - '+str(int(max(wcr_formations_df['year'])))+')', y=0.935, fontsize=11);

    # WCR-like eddies
    ax.plot(wcr_formations_df['year'],wcr_formations_df['all_zones'],'-o',color='#EE4D4D')
#     ax.plot([min_year_wcr,2000],[pre_2000_wcr['all_zones'].mean(),pre_2000_wcr['all_zones'].mean()], color='r', linestyle='-'); # pre_2000 average
#     ax.plot([2000,max_year_wcr],[post_2000_wcr['all_zones'].mean(),post_2000_wcr['all_zones'].mean()], color='r', linestyle='-'); # post_2000 average

    # axes formatting
    ax.set_xlabel('Years',fontweight='bold')
    ax.set_ylabel('# of Formations',fontweight='bold');
    ax.set_ylim(ylim_min, ylim_max) # 1993-2017
#     ax.set_xlim(1980,2020) # to compare to fig. 2 Silver et al., 2021

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#EE4D4D', lw=4), Line2D([0], [0], color='#388FEF', lw=4)]
    ax.legend(custom_lines, ['WCR'], loc='upper left');
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------


