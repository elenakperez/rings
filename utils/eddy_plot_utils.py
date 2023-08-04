"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to plot Chelton Eddy Tracks ()

    1) eddy_map_year : returns a map of all the eddie tracks for a given year
    2) zone_bar_plot : return bar plot of formations by zone for rings
    3)               :
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
from utils.eddy_data_utils import *

#-------------------------------------------------------------------------------------------------------------------------------
# 1) takes in a pandas dataframe of eddies, specific year, and title and returns a map of the eddy tracks for that year

def eddy_map_year(eddy_df, bathy, year, title):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy_df (DataFrame) : pandas dataframe of eddies in the larger region, e.g. northwest atlantic
        year (int)          : which year you want to plot, e.g. 1994
        title (str)         : title of the figure, e.g. 'Eddy Tracks (1994)'

        
    Output:
        * returns a map of all the eddie tracks for a given year
    
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    annual_eddy = eddy_df[(eddy_df['time'].dt.year == year)]
    annual_cyclonic = eddy_df[(eddy_df['time'].dt.year == year) & (eddy_df['cyclonic_type']==-1)]
    annual_anticyclonic = eddy_df[(eddy_df['time'].dt.year == year) & (eddy_df['cyclonic_type']==+1)]

    # # define bounds for northwest Atlantic region
    # x_bnds = [-82,-48] # lon
    # y_bnds = [24,53] # lat

    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon
    y_bnds = [29,45] # lat

    # plt.plot(annual_eddy['longitude'].iloc[start:end],annual_eddy['latitude'])
    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(100)
    fig.suptitle(title+'('+str(year)+')', fontsize=25, y=0.875)

    for i in np.array(annual_eddy['track'].unique()):
        eddy = annual_eddy[annual_eddy['track']==i]
        if ((eddy['cyclonic_type']==1).all()): # if anti-cyclonic 
            ax.plot(eddy['longitude'],eddy['latitude'],color='red')

        elif ((eddy['cyclonic_type']==-1).all()): # if cyclonic 
            ax.plot(eddy['longitude'],eddy['latitude'],color='blue')

         # plot start
    #     ax.plot(eddy['longitude'].iloc[0]-360,eddy['latitude'].iloc[0],color='green',marker='o')

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')    

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------
# ) return bar plot of formations by zone for rings



def zone_bar_plot(ccr_formations_df, wcr_formations_df, title, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
     Input:
        ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones 
        wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones 
        title (Str)                   : title of the figure, e.g. 'Ring Formations Zones 1-4 (1993-2020)'
        fig_quality (Int)             : integer of what dpi the image will be set to (e.g., 100 dpi)


    Output:
        * returns a map of all the eddie tracks for a given year

     
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

    ## !OLD!
#     # data
#     anticyclones_count = [len((zone1_eddies[zone1_eddies['cyclonic_type']==+1])['track'].unique()), len((zone2_eddies[zone2_eddies['cyclonic_type']==+1])['track'].unique()), len((zone3_eddies[zone3_eddies['cyclonic_type']==+1])['track'].unique()), len((zone4_eddies[zone4_eddies['cyclonic_type']==+1])['track'].unique())]
#     cyclones_count = [len((zone1_eddies[zone1_eddies['cyclonic_type']==-1])['track'].unique()), len((zone2_eddies[zone2_eddies['cyclonic_type']==-1])['track'].unique()), len((zone3_eddies[zone3_eddies['cyclonic_type']==-1])['track'].unique()), len((zone4_eddies[zone4_eddies['cyclonic_type']==-1])['track'].unique())]
    
    # # plotting
    # ax.bar(ind, cyclones_count , width, label='Cyclones', color='#388FEF')
    # ax.bar(ind + width, anticyclones_count, width, label='Anti-cyclones', color='#EE4D4D')

    # plotting
    ax.bar(ind, ccr_formations_df[['zone_1','zone_2','zone_3','zone_4']].sum() , width, label='CCR', color='#388FEF')
    ax.bar(ind + width, wcr_formations_df[['zone_1','zone_2','zone_3','zone_4']].sum(), width, label='WCR', color='#EE4D4D')

    ax.set_xticks([0.15,1.15,2.15,3.15])
    ax.set_xticklabels(zones)
    ax.set_ylabel('Number of Formations');
    ax.legend();
    
    return fig, ax
# example call in notebook
# zone_bar_plot(zone_ccr_annual_formations, zone_wcr_annual_formations, 'Ring Formations ', 125);

#-------------------------------------------------------------------------------------------------------------------------------
# ) 

def rings_interannual_variability(ccr_formations_df, wcr_formations_df, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones 
        wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones 
        fig_quality (Int)             : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)                  : returns a time-series of interannual variability of rings
     
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
#     fig.suptitle('WCR vs. CCR Formations ('+str(int(min(ccr_formations_df['year'])))+' - '+str(int(max(ccr_formations_df['year'])))+')', y=0.935, fontsize=11);

    # anticyclones
    ax.plot(wcr_formations_df['year'],wcr_formations_df['all_zones'],'-o',color='r')
    ax.plot([min_year_wcr,2000],[pre_2000_wcr['all_zones'].mean(),pre_2000_wcr['all_zones'].mean()], color='r', linestyle='-'); # pre_2000 average
    ax.plot([2000,max_year_wcr],[post_2000_wcr['all_zones'].mean(),post_2000_wcr['all_zones'].mean()], color='r', linestyle='-'); # post_2000 average

    # cyclones
    ax.plot(ccr_formations_df['year'],ccr_formations_df['all_zones'],'-o',color='b')
    ax.plot([min_year_ccr,2000],[pre_2000_ccr['all_zones'].mean(),pre_2000_ccr['all_zones'].mean()], color='b', linestyle='-'); # pre_2000 average
    ax.plot([2000,max_year_ccr],[post_2000_ccr['all_zones'].mean(),post_2000_ccr['all_zones'].mean()], color='b', linestyle='-'); # post_2000 average

    # axes formatting
    ax.set_xlabel('Years',fontweight='bold')
    ax.set_ylabel('Number of Ring Formations',fontweight='bold');
#     ax.set_ylim(5,45)
#     ax.set_xlim(1980,2020) # to compare to fig. 2 Silver et al., 2021

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['WCR', 'CCR'], loc='upper left');
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# ) 

def eddy_interannual_variability(eddy_df, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones 
        wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones 
        fig_quality (Int)             : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)                  : returns a time-series of interannual variability of rings
     
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 

    min_year = int(min(eddy_df['year']))
    max_year = int(max(eddy_df['year']))

    pre_2000 = eddy_df[eddy_df['year']<2000] # mean rings count pre-2000
    post_2000 = eddy_df[eddy_df['year']>=2000] # mean ring count post-2000

    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle('Anti-cyclonic vs. Cyclonic Eddy Formations (1993-2020)', fontsize=11, y=0.935, fontweight='bold')

    # anticyclones
    ax.plot(eddy_df['year'].iloc[0:27],eddy_df['anticyclonic'].iloc[0:27],'-o',color='r')
    ax.plot([min_year,2000],[pre_2000['anticyclonic'].mean(),pre_2000['anticyclonic'].mean()], color='r', linestyle='-'); # pre_2000 average
    ax.plot([2000,max_year],[post_2000['anticyclonic'].mean(),post_2000['anticyclonic'].mean()], color='r', linestyle='-'); # post_2000 average

    # cyclones
    ax.plot(eddy_df['year'].iloc[0:27],eddy_df['cyclonic'].iloc[0:27],'-o',color='b')
    ax.plot([min_year,2000],[pre_2000['cyclonic'].mean(),pre_2000['cyclonic'].mean()], color='b', linestyle='-'); # pre_2000 average
    ax.plot([2000,max_year],[post_2000['cyclonic'].mean(),post_2000['cyclonic'].mean()], color='b', linestyle='-'); # post_2000 average

    # axes formatting
    ax.set_xlabel('Years',fontweight='bold')
    ax.set_ylabel('Number of Ring Formations',fontweight='bold');
# #     ax.set_ylim(5,45)
#     ax.set_xlim(1980,2020) # to compare to fig. 2 Silver et al., 2021

    return fig, ax

#-------------------------------------------------------------------------------------------------------------------------------
# ) 

def ring_monthly_formations(ccr_formations_df, wcr_formations_df, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones 
        wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones 
        fig_quality (Int)             : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)                  : returns a time-series of interannual variability of rings
     
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    
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
        
    # lines
    ccr_minus_std = np.subtract(ccr_monthly_avg,ccr_monthly_std)
    wcr_minus_std = np.subtract(wcr_monthly_avg,wcr_monthly_std)
    ccr_plus_std = np.add(ccr_monthly_avg,ccr_monthly_std)
    wcr_plus_std = np.add(wcr_monthly_avg,wcr_monthly_std)

    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle('WCR vs. CCR Monthly Formations', y=0.935, fontsize=11);

    # plot all zones mean monthly formations
#     ax.plot(ccr_formations_df['month'],ccr_formations_df['all_zones_mean'],'-o',color='blue')
#     ax.plot(wcr_formations_df['month'],wcr_formations_df['all_zones_mean'],'-o',color='r');
    ax.plot(ccr_formations_df['month'].drop_duplicates(),ccr_monthly_avg,'-o',color='blue');
#     ax.plot(ccr_formations_df['month'].drop_duplicates(),np.subtract(ccr_monthly_avg,ccr_monthly_std),color='blue')
    
    ax.plot(wcr_formations_df['month'].drop_duplicates(),wcr_monthly_avg,'-o',color='r');
    
#     # standard deivation shading
#     # ccrs
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), ccr_monthly_avg, ccr_minus_std, alpha=0.3, color='#60BEFA')
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), ccr_monthly_avg, ccr_plus_std, alpha=0.3, color='#60BEFA')
#     # wcrs
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_monthly_avg, wcr_minus_std, alpha=0.3, color='#FF938F')
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_monthly_avg, wcr_plus_std, alpha=0.3, color='#FF938F')


    zones = ['zone_1', 'zone_2', 'zone_3', 'zone_4']

#     # plot zones 1-4 mean
#     for z in zones:
#         ax.plot(wcr_formations_df['month'],wcr_formations_df[z+'_mean'],'-o',color='#FF938F')
#         ax.plot(ccr_formations_df['month'],ccr_formations_df[z+'_mean'],'-o',color='#60BEFA')

    # axes formatting
    ax.set_xlabel('Months', fontweight='bold')
    ax.set_ylabel('Average Monthly Ring Formation',fontweight='bold');
#     ax.set_ylim(0,5)

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['WCR', 'CCR']);
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------
# ) 

def zonal_ring_monthly_formations(ccr_formations_df, wcr_formations_df, fig_quality, which_zone):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones 
        wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones 
        fig_quality (Int)             : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)                  : returns a time-series of interannual variability of rings
     
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    
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
        
    # lines
    ccr_minus_std = np.subtract(ccr_monthly_avg,ccr_monthly_std)
    wcr_minus_std = np.subtract(wcr_monthly_avg,wcr_monthly_std)
    ccr_plus_std = np.add(ccr_monthly_avg,ccr_monthly_std)
    wcr_plus_std = np.add(wcr_monthly_avg,wcr_monthly_std)

    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle(which_zone+' Formations', y=0.935, fontsize=11);

    # plot all zones mean monthly formations
#     ax.plot(ccr_formations_df['month'],ccr_formations_df['all_zones_mean'],'-o',color='blue')
#     ax.plot(wcr_formations_df['month'],wcr_formations_df['all_zones_mean'],'-o',color='r');

#     ax.plot(ccr_formations_df['month'].drop_duplicates(),ccr_monthly_avg,'-o',color='blue');
#     ax.plot(wcr_formations_df['month'].drop_duplicates(),wcr_monthly_avg,'-o',color='r');
    
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), ccr_monthly_avg, ccr_minus_std, alpha=0.3, color='#60BEFA')
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), ccr_monthly_avg, ccr_plus_std, alpha=0.3, color='#60BEFA')

#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_monthly_avg, wcr_minus_std, alpha=0.3, color='#FF938F')
#     ax.fill_between(wcr_formations_df['month'].drop_duplicates(), wcr_monthly_avg, wcr_plus_std, alpha=0.3, color='#FF938F')


    zones = ['zone_1', 'zone_2', 'zone_3', 'zone_4']

#     # plot zones 1-4 mean
#     for z in zones:
#         ax.plot(wcr_formations_df['month'],wcr_formations_df[z+'_mean'],'-o',color='r')
#         ax.plot(ccr_formations_df['month'],ccr_formations_df[z+'_mean'],'-o',color='b')

    ax.plot(wcr_formations_df['month'],wcr_formations_df[which_zone+'_mean'],'-o',color='r')
    ax.plot(ccr_formations_df['month'],ccr_formations_df[which_zone+'_mean'],'-o',color='b')

    # axes formatting
    ax.set_xlabel('Months', fontweight='bold')
    ax.set_ylabel('Average Monthly Ring Formation',fontweight='bold');
    ax.set_ylim(0,2)

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['WCR', 'CCR']);
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# )

# modified function so that anti-cylcones north of the GS are "WCRs" and cyclones south of the GS are "CCRs"
def eddy_filt_map_year(eddy_df, bathy, year, title, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        eddy_df (DataFrame) : pandas dataframe of eddies in the larger region, e.g. northwest atlantic
        year (int)          : which year you want to plot, e.g. 1994
        title (str)         : title of the figure, e.g. 'Eddy Tracks (1994)'
        fig_quality (Int)             : quality of the figure (e.g. 100 dpi) 

    Output:
        * returns a map of all the zone eddy tracks for a given year
     
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    
    annual_eddy = eddy_df[(eddy_df['time'].dt.year == year)]
    annual_cyclonic = eddy_df[(eddy_df['time'].dt.year == year) & (eddy_df['cyclonic_type']==-1)]
    annual_anticyclonic = eddy_df[(eddy_df['time'].dt.year == year) & (eddy_df['cyclonic_type']==+1)]

    # # define bounds for northwest Atlantic region
    # x_bnds = [-82,-48] # lon
    # y_bnds = [24,53] # lat

    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon
    y_bnds = [29,45] # lat

    # plt.plot(annual_eddy['longitude'].iloc[start:end],annual_eddy['latitude'])
    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title+'('+str(year)+')', fontsize=25, y=0.875)
    
    # colors for lines
    colors = sns.color_palette("hls", 12)
    monthColors = ['#9ecae1','#4393c3','#2166ac','#8c510a','#bf812d','#f46d43','#d73027','#971700','#fddbc7','#d9f0d3','#5aae61','#1FC1F2']
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
#         # gs paths
#     # plot all months, years paths
#     for i in (np.arange(1993,2019)):
#         yr = i
#         for month in range(12):
#             if yr!=2018:
#                 ax.plot(get_gs_month(yr,month-1)[0],get_gs_month(yr,month-1)[1],label=months[month], color='#D8D8D8',linewidth=3); #color=monthColors[month]
#             else:
#                 if month<5:
#                     ax.plot(get_gs_month(yr,month-1)[0],get_gs_month(yr,month-1)[1],label=months[month], color='#D8D8D8',linewidth=3); #color=monthColors[month]
    
    for month in range(12):
        if year!=2018:
            ax.plot(get_gs_month(year,month-1)[0],get_gs_month(year,month-1)[1],label=months[month], color='#D8D8D8',linewidth=3); #color=monthColors[month]
        else:
            if month<5:
                ax.plot(get_gs_month(year,month-1)[0],get_gs_month(year,month-1)[1],label=months[month], color='#D8D8D8',linewidth=3); #color=monthColors[month]

    # plot gs path for given year
    ax.plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
    
    
    for i in np.array(annual_eddy['track'].unique()):
        eddy = annual_eddy[annual_eddy['track']==i]
        if (eddy['cyclonic_type']==1).all(): # if anti-cyclonic & north of the gulf stream
            ax.plot(eddy['longitude'],eddy['latitude'],color='red')
#             # plot start
#             ax.plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o')


        elif (eddy['cyclonic_type']==-1).all(): # if cyclonic & NOT north of the gulf stream (i.e. south)
            ax.plot(eddy['longitude'],eddy['latitude'],color='blue')
#             # plot start
#             ax.plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o')

       
    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  
    
    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['Anti-cyclonic', 'Cyclonic'])

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------
# ) makes a GIF of all the formations of CCRs, WCRs from 1993 to 2017, stepping through at yearly intervals

def make_gif(wcr_df, ccr_df):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
    
    Output:
        * no output, saves GIF
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    for i in (np.arange(1993,2018)):
        eddy_filt_map_zones_year(wcr_df, ccr_df, 'Ring Tracks', i)[0].savefig('//Users/elenaperez/Desktop/chatts/figures/yearly_rings/'+str(i)+'.png',bbox_inches='tight');

    images,image_file_names = [],[]
    for file_name in os.listdir('/Users/elenaperez/Desktop/chatts/figures/yearly_rings'):
        if file_name.endswith('.png'):
            image_file_names.append(file_name)       
    sorted_files = sorted(image_file_names,key=lambda x: int(os.path.splitext(x)[0]))

    for filename in sorted_files:
        images.append(imageio.imread('/Users/elenaperez/Desktop/chatts/figures/yearly_rings/'+filename))
    imageio.mimsave('/Users/elenaperez/Desktop/chatts/figures/all_years_ring_tracks.gif', images,duration=1)
    
# example call
# make_gif(zone_wcrs, zone_ccrs);


#-------------------------------------------------------------------------------------------------------------------------------
# )
def eddy_filt_map_grid(wcr_df, ccr_df, bathy, title):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
        bathy (xarray)      : dataarray of northwest atlantic bathymetry 
        title (String)      : title of the figure
    
    Output:
        fig (Figure)        : returns 5x5 grid of annual ring formations

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon, NWA: [-82,-48]
    y_bnds = [29,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()
    fig,ax = plt.subplots(nrows=5,ncols=5,sharex='col',sharey='row',constrained_layout=False,figsize=(15, 17), subplot_kw = dict(projection=proj))
    
    
    track_year = 1993
    for i in range(5):
        for j in range(5):

            
#             # plot gs path for all years
#             for year in (np.arange(1993,2018)):
#                 ax[i,j].plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
                
                    # plot all monthly paths for given year
            for month in range(12):
                if track_year!=2018:
                    ax[i,j].plot(get_gs_month(track_year,month-1)[0],get_gs_month(track_year,month-1)[1], color='grey',linewidth=3); #color=monthColors[month]
                else:
                    if month<5:
                        ax[i,j].plot(get_gs_month(track_year,month-1)[0],get_gs_month(track_year,month-1)[1], color='grey',linewidth=3); #color=monthColors[month]

#             # add zone lines
#             ax[i,j].plot([-75,29],[-75,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 1
#             ax[i,j].plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
#             ax[i,j].plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
#             ax[i,j].plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4
#         #     ax[i,j].plot([-55,29],[-55,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4.5
        
            # plot tracks *ALL YEARS*
                ## WCRs ##
            for k in np.array(wcr_df['track'].unique()):
                eddy = wcr_df[wcr_df['track']==k]
                ax[i,j].plot(eddy['longitude'],eddy['latitude'],color='#FF938F',alpha=0.1) # track
        #         ax[i,j].plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o') # formation
        #         ax[i,j].plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='red',marker='o') # demise

            annual_eddy = wcr_df[(wcr_df['time'].dt.year == track_year)]
            for l in np.array(annual_eddy['track'].unique()):
                eddy = annual_eddy[annual_eddy['track']==l]
                ax[i,j].plot(eddy['longitude'],eddy['latitude'],color='red') # track

        
            ## CCRs ##
            for m in np.array(ccr_df['track'].unique()):
                eddy = ccr_df[ccr_df['track']==m]
                ax[i,j].plot(eddy['longitude'],eddy['latitude'],color='#60BEFA', alpha=0.1) # track
        #         ax[i,j].plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o') # formation
        #         ax[i,j].plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='red',marker='o') # demise
            
            annual_eddy = ccr_df[(ccr_df['time'].dt.year == track_year)]
            for n in np.array(annual_eddy['track'].unique()):
                eddy = annual_eddy[annual_eddy['track']==n]
                ax[i,j].plot(eddy['longitude'],eddy['latitude'],color='blue') # track
                
            # axes formatting
            ax[i,j].set_title(str(track_year),fontweight="bold")
            ax[i,j].set_xlabel('Longitude')
            ax[i,j].set_ylabel('Latitude')
            ax[i,j].contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
            ax[i,j].coastlines(resolution='50m',color='gray')
            ax[i,j].set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
            ax[i,j].add_feature(cartopy.feature.LAND, color='lightgray')  

            # axes formatting
            ax[i,j].set_xlabel('Longitude')
            ax[i,j].set_ylabel('Latitude')
            ax[i,j].contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
            ax[i,j].coastlines(resolution='50m',color='gray')
            ax[i,j].set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
            ax[i,j].add_feature(cartopy.feature.LAND, color='lightgray')  

#             # custom legend
#             from matplotlib.lines import Line2D
#             custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
#             ax.legend(custom_lines, ['Warm Core Rings', 'Cold Core Rings'])

             # gridlines
            if j == 0:
                gl = ax[i,j].gridlines(crs=proj,draw_labels=True)
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.yformatter = LATITUDE_FORMATTER
                gl.xformatter = LONGITUDE_FORMATTER
                if i!=4:
                    gl.xlabels_bottom = False
            if i==5:
                gl = ax[i,j].gridlines(crs=proj,draw_labels=True)
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.ylabels_left = False
                gl.yformatter = LATITUDE_FORMATTER
                gl.xformatter = LONGITUDE_FORMATTER
            else:
                gl = ax[i,j].gridlines(crs=proj,draw_labels=True)
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xlabels_bottom = False
                gl.ylabels_left = False
                gl.yformatter = LATITUDE_FORMATTER
                gl.xformatter = LONGITUDE_FORMATTER
                
            track_year=track_year+1
    
#     fig.tight_layout()
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# ) map of the ring tracks for a given year

def ring_filt_map_zones_year(wcr_df, ccr_df, bathy, title, track_year, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
        bathy (xarray)      : dataarray of northwest atlantic bathymetry 
        title (String)      : title of the figure
        track_year (Int)    : year of interest, e.g. 1993
        fig_quality (Int)   : quality of the figure (e.g. 100 dpi) 
    
    Output:
        fig (Figure)        : returns map of ring tracks for a given year 
                                (all tracks in light color and yearly tracks in bold color)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon, NWA: [-82,-48]
    y_bnds = [29,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title+' ('+str(track_year)+')', fontsize=25, y=0.875)
    
    # colors for lines
    colors = sns.color_palette("hls", 12)
    monthColors = ['#9ecae1','#4393c3','#2166ac','#8c510a','#bf812d','#f46d43','#d73027','#971700','#fddbc7','#d9f0d3','#5aae61','#1FC1F2']
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
#     # plot gs path for all years
#     for year in (np.arange(1993,2018)):
#         ax.plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
        
        # plot all monthly paths for given year
    for month in range(12):
        if track_year!=2018:
            ax.plot(get_gs_month(track_year,month-1)[0],get_gs_month(track_year,month-1)[1], color='grey',linewidth=3); #color=monthColors[month]
        else:
            if month<5:
                ax.plot(get_gs_month(track_year,month-1)[0],get_gs_month(track_year,month-1)[1], color='grey',linewidth=3); #color=monthColors[month]

    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 1
    ax.plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
    ax.plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
    ax.plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4
#     ax.plot([-55,29],[-55,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4.5

    for i in np.array(wcr_df['track'].unique()):
        eddy = wcr_df[wcr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='#FF938F',alpha=0.3) # track
#         ax.plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o') # formation
#         ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='red',marker='o') # demise

    ## CCRs ##
    for i in np.array(ccr_df['track'].unique()):
        eddy = ccr_df[ccr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='#60BEFA', alpha=0.3) # track
#         ax.plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o') # formation
#         ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='red',marker='o') # demise


    # plot YEAR tracks
    
    ## WCRs ##
    annual_eddy = wcr_df[(wcr_df['time'].dt.year == track_year)]
    for i in np.array(annual_eddy['track'].unique()):
        eddy = annual_eddy[annual_eddy['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='red',linewidth=2) # track

    ## CCRs ##
    annual_eddy = ccr_df[(ccr_df['time'].dt.year == track_year)]
    for i in np.array(annual_eddy['track'].unique()):
        eddy = annual_eddy[annual_eddy['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='blue',linewidth=2) # track
#         ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='red',marker='o') # demise
#         ax.plot(get_eddy_formation_loc(eddy)[0],get_eddy_formation_loc(eddy)[1],color='green',marker='o') # formation


       
    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  
    
    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['Warm Core Rings', 'Cold Core Rings'])

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# ) creates a map of ring tracks for all years in the northwest atlantic

def ring_filt_map_zones(wcr_df, ccr_df, bathy, title, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
        bathy (xarray)      : dataarray of northwest atlantic bathymetry 
        title (String)      : title of the figure
        fig_quality (Int)   : quality of the figure (e.g. 100 dpi) 
    
    Output:
        fig (Figure)        : returns map of ring tracks for all years

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon, NWA: [-82,-48]
    y_bnds = [29,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title, fontsize=25, y=0.875)
    
    # colors for lines
    colors = sns.color_palette("hls", 12)
    monthColors = ['#9ecae1','#4393c3','#2166ac','#8c510a','#bf812d','#f46d43','#d73027','#971700','#fddbc7','#d9f0d3','#5aae61','#1FC1F2']
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    # plot gs path for all years
    for year in (np.arange(1993,2018)):
        ax.plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 1
    
    ## WCRs ##
    for i in np.array(wcr_df['track'].unique()):
        eddy = wcr_df[wcr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='red')

    ## CCRs ##
    for i in np.array(ccr_df['track'].unique()):
        eddy = ccr_df[ccr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='blue')
       
    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  
    
    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['Warm Core Rings', 'Cold Core Rings'])

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# ) creates a map of EDDY tracks for all years in the northwest atlantic

def eddy_filt_map_zones(eddy_df, bathy, title, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
        bathy (xarray)      : dataarray of northwest atlantic bathymetry 
        title (String)      : title of the figure
        fig_quality (Int)   : quality of the figure (e.g. 100 dpi) 
    
    Output:
        fig (Figure)        : returns map of eddy tracks for all years

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon, NWA: [-82,-48]
    y_bnds = [29,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(title, fontsize=25, y=0.875)
    
    # colors for lines
    colors = sns.color_palette("hls", 12)
    monthColors = ['#9ecae1','#4393c3','#2166ac','#8c510a','#bf812d','#f46d43','#d73027','#971700','#fddbc7','#d9f0d3','#5aae61','#1FC1F2']
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    # plot gs path for all years
    for year in (np.arange(1993,2018)):
        ax.plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 1
    
    ## anti-cyclonic ##
    for i in np.array(eddy_df['track'].unique()):
        eddy = eddy_df[eddy_df['track']==i]        
        if (eddy['cyclonic_type']==1).all(): # if anti-cyclonic & north of the gulf stream
            ax.plot(eddy['longitude'],eddy['latitude'],color='red')
        elif (eddy['cyclonic_type']==-1).all():
            ax.plot(eddy['longitude'],eddy['latitude'],color='blue') #cyclonic

       
    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  
    
    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='r', lw=4), Line2D([0], [0], color='b', lw=4)]
    ax.legend(custom_lines, ['Anti-cylonic', 'Cyclonic'])

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax


#------------------------------------------------------------------------------------------------------------------------------
# ) map of demise locations for all rings 

def demise_filt_map_zones(wcr_df, ccr_df, bathy, title, ring_type, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
        bathy (xarray)      : dataarray of northwest atlantic bathymetry 
        title (String)      : title of the figure
        ring_type (String)  : 'ccr' for CCRs, 'wcr' for WCRs, and 'all' for both
        fig_quality (Int)   : quality of the figure (e.g. 100 dpi) 
    
    Output:
        fig (Figure)        : returns map of eddy tracks for all years

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon, NWA: [-82,-48]
    y_bnds = [29,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(ring_type+title, fontsize=25, y=0.875)
    
    # colors for lines
    colors = sns.color_palette("hls", 12)
    monthColors = ['#9ecae1','#4393c3','#2166ac','#8c510a','#bf812d','#f46d43','#d73027','#971700','#fddbc7','#d9f0d3','#5aae61','#1FC1F2']
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    # plot gs path for all years
    for year in (np.arange(1993,2018)):
        ax.plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 1
    ax.plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
    ax.plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
    ax.plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4
#     ax.plot([-55,29],[-55,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4.5

    if (ring_type=='wcr') or (ring_type=='all'):
        ## WCRs ##
        for i in np.array(wcr_df['track'].unique()):
            eddy = wcr_df[wcr_df['track']==i]
    #         # track
    #         ax.plot(eddy['longitude'],eddy['latitude'],color='red')

            eddy_lifespan = get_eddy_lifespan(eddy)

#     #         # formation
#             if (eddy_lifespan <= 150):
#                 ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='green',marker='o',label='Longlived Rings')
#             else:
#                 ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='orange',marker='o',label='Longlived Rings')

            # demise    
            if (eddy_lifespan <= 150):
                ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='green',marker='o',label='Longlived Rings')
            else:
                ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='orange',marker='o',label='Longlived Rings')
     
    if (ring_type=='ccr') or (ring_type=='all'):
        ## CCRs ##
        for i in np.array(ccr_df['track'].unique()):
            eddy = ccr_df[ccr_df['track']==i]
    #         # track
    #         ax.plot(eddy['longitude'],eddy['latitude'],color='blue')

            eddy_lifespan = get_eddy_lifespan(eddy)

#     #         # formation
#             if (eddy_lifespan <= 150):
#                 ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='green',marker='o',label='Longlived Rings')
#             else:
#                 ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='orange',marker='o',label='Longlived Rings')

            # demise    
            if (eddy_lifespan <= 150):
                ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='green',marker='o',label='Longlived Rings')
            else:
                ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='orange',marker='o',label='Longlived Rings')      
    
    
    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  
    
    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='g', lw=4), Line2D([0], [0], color='orange', lw=4)]
    ax.legend(custom_lines, ['Shortlived WCRs', 'Longlived WCRs'],loc='upper left')
#     ax.legend()

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax


#------------------------------------------------------------------------------------------------------------------------------
# ) map of formation locations for all rings

def formation_filt_map_zones(wcr_df, ccr_df, bathy, title, ring_type, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        wcr_df (DataFrame)  : pandas dataframe of warm core rings
        ccr_df (DataFrame)  : pandas dataframe of cold core rings
        bathy (xarray)      : dataarray of northwest atlantic bathymetry 
        title (String)      : title of the figure
        ring_type (String)  : 'ccr' for CCRs, 'wcr' for WCRs, and 'all' for both
        fig_quality (Int)   : quality of the figure (e.g. 100 dpi) 
    
    Output:
        fig (Figure)        : returns map of eddy tracks for all years

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # gangopadhyay census bounds
    x_bnds = [-85,-55] # lon, NWA: [-82,-48]
    y_bnds = [29,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    fig.suptitle(ring_type+title, fontsize=25, y=0.875)
    
    # colors for lines
    colors = sns.color_palette("hls", 12)
    monthColors = ['#9ecae1','#4393c3','#2166ac','#8c510a','#bf812d','#f46d43','#d73027','#971700','#fddbc7','#d9f0d3','#5aae61','#1FC1F2']
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    # plot gs path for all years
    for year in (np.arange(1993,2018)):
        ax.plot(get_gs_year(year)[0],get_gs_year(year)[1],label=year, color='gray', linewidth=2);
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 1
    ax.plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
    ax.plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
    ax.plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4
#     ax.plot([-55,29],[-55,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 4.5

    if (ring_type=='wcr') or (ring_type=='all'):
        ## WCRs ##
        for i in np.array(wcr_df['track'].unique()):
            eddy = wcr_df[wcr_df['track']==i]
    #         # track
    #         ax.plot(eddy['longitude'],eddy['latitude'],color='red')

            eddy_lifespan = get_eddy_lifespan(eddy)

    #         # formation
            if (eddy_lifespan <= 150):
                ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='green',marker='o',label='Longlived Rings')
            else:
                ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='orange',marker='o',label='Longlived Rings')

    #         # demise    
    #         if (eddy_lifespan <= 150):
    #             ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='green',marker='o',label='Longlived Rings')
    #         else:
    #             ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='orange',marker='o',label='Longlived Rings')
     
    if (ring_type=='ccr') or (ring_type=='all'):
        ## CCRs ##
        for i in np.array(ccr_df['track'].unique()):
            eddy = ccr_df[ccr_df['track']==i]
    #         # track
    #         ax.plot(eddy['longitude'],eddy['latitude'],color='blue')

            eddy_lifespan = get_eddy_lifespan(eddy)

    #         # formation
            if (eddy_lifespan <= 150):
                ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='green',marker='o',label='Longlived Rings')
            else:
                ax.plot(eddy['longitude'].iloc[0], eddy['latitude'].iloc[0], color='orange',marker='o',label='Longlived Rings')

    #         # demise    
    #         if (eddy_lifespan <= 150):
    #             ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='green',marker='o',label='Longlived Rings')
    #         else:
    #             ax.plot(eddy['longitude'].iloc[-1], eddy['latitude'].iloc[-1], color='orange',marker='o',label='Longlived Rings')
            

       
    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  
    
    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='g', lw=4), Line2D([0], [0], color='orange', lw=4)]
    ax.legend(custom_lines, ['Shortlived WCRs', 'Longlived WCRs'],loc='upper left')
#     ax.legend()

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax


#----------------------------------------------------------------------------------------------------------------------------
# ) 3x4 grid of histograms that show seasonal distribution of ring formations

def hist_ring_forms_by_month(ring_df, title, ring_type, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
     Input:
        ring_df (DataFrame)   : pandas dataframe of formation counts for a specific ring type
        title (Str)           : title of the figure, e.g. 'Ring Formations Zones 1-4 (1993-2020)'
        ring_type (Str)       : 'ccr' for CCRs, 'wcr' for WCRs
        fig_quality (Int)     : integer of what dpi the image will be set to (e.g., 100 dpi)


    Output:
        * returns grid of histograms that show distributions of number of formations vs. frequency of formation

     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    fig,ax = plt.subplots(nrows=4,ncols=3,constrained_layout=False,figsize=(15, 17))

    # figure
    fig.set_dpi(fig_quality)
    fig.suptitle(title, fontsize=20,fontweight='bold', y=1,x=0.54)
    fig.subplots_adjust(hspace=0.1, wspace=0.1)

    month = 1
    for i in range(4):
        for j in range(3):
            if ring_type=='ccr':
                ax[i,j].hist((ring_df[ring_df['month']==month])['all_zones'],color='blue',bins=np.arange(1,11));
            elif ring_type=='wcr':
                ax[i,j].hist((ring_df[ring_df['month']==month])['all_zones'],color='red',bins=np.arange(1,11));


            #axes
            datetime_object = datetime.datetime.strptime(str(month), "%m")
            month_name = datetime_object.strftime("%b")
            ax[i,j].set_title(month_name)
            ax[i,j].set_ylim(0,14)
            month+=1

    fig.tight_layout()
    
    return fig,ax

#------------------------------------------------------------------------------------------------------------------------------
# ) plots seasonal variability of formations by zones for WCRs and CCRs
def ring_seasonal_var_by_zone(ccr_formations_df, wcr_formations_df, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        ccr_formations_df (DataFrame) : pandas dataframe of cold core rings in all zones 
        wcr_formations_df (DataFrame) : pandas dataframe of warm core rings in all zones 
        fig_quality (Int)             : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)                  : returns a time-series of interannual variability of rings
     
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 

    fig,ax = plt.subplots(nrows=2,ncols=2,sharex='col',sharey='row',constrained_layout=False)

    # figure
    fig.set_dpi(fig_quality)
    # fig.suptitle('WCR vs. CCR Monthly Formations by Zone', fontsize=20,fontweight='bold', y=1,x=0.54)
    fig.subplots_adjust(hspace=0.1, wspace=0.1)

    zones = ['zone_1', 'zone_2', 'zone_3', 'zone_4']
    ctr = 0
    for i in range(2):
        for j in range(2):
            ax[i,j].plot(wcr_formations_df['month'],wcr_formations_df[zones[ctr]+'_mean'],'-o',color='r',label=zones[ctr])
            ax[i,j].plot(ccr_formations_df['month'],ccr_formations_df[zones[ctr]+'_mean'],'-o',color='b',label=zones[ctr])
    #         ax[i,j].set_title(zones[ctr])

            ax[i,j].legend()

    #         axes formatting
            if i==1:
                ax[i,j].set_xlabel('Months', fontweight='bold')

            fig.text(0.04, 0.5, 'Average Monthly Ring Formation', va='center', rotation='vertical')

            ctr +=1 

            ax[i,j].set_ylim(0,1.5)
            
    return fig,ax


#------------------------------------------------------------------------------------------------------------------------------
# )

# ) 
def lifespan_bar_chart_zone(short_wcr_df, long_wcr_df, short_ccr_df, long_ccr_df, title, ring_type, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        ccr_df (DataFrame)         : pandas dataframe of cold core rings in all zones 
        wcr_df (DataFrame)         : pandas dataframe of warm core rings in all zones 
        title (String)             : title of the figure
        ring_type (String)         : type of ring (or eddy) to be plotted (e.g., 'ccr', 'wcr', 'all')
        fig_quality (Integer)      : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)               : returns a comparative time-series of interannual variability of rings
     
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    
    # create figure
    fig,ax = plt.subplots();
    fig.suptitle(title, y=0.95, fontsize=14);
    fig.set_dpi(fig_quality)
    
    # calculate optimal width of bar plots' bars
    N = 4 # number of bar pairs i want
    ind = np.arange(N)
    width = np.min(np.diff(ind))/3 # Calculate optimal width
    
    zones = ['Zone 1', 'Zone 2', 'Zone 3', 'Zone 4']
    
    # plotting
    # CCR
    if ring_type=='ccr':
        ax.bar(ind, short_ccr_df.iloc[0] , width, label='Short-lived', color='#388FEF')
        ax.bar(ind + width, long_ccr_df.iloc[0], width, label='Long-lived', color='blue')

    if ring_type=='wcr':
        ax.bar(ind, short_wcr_df.iloc[0], width, label='Short-lived', color='#EE4D4D')
        ax.bar(ind + width, long_wcr_df.iloc[0], width, label='Long-lived', color='orange')
    
    ax.set_xticks([0.15,1.15,2.15,3.15])
    ax.set_xticklabels(zones)
    ax.set_ylabel('Number of Formations');
    ax.legend();

    
    return fig, ax

#------------------------------------------------------------------------------------------------------------------------------
# ) compare clark chart annual formations to chelton track annual formations
def compare_clark_chelton(clark_ccr, clark_wcr, chelton_eddies, chelton_ccr, chelton_wcr, title, ring_type, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Input:
        clark_ccr (DataFrame)          : pandas dataframe of Clark Chart cold core rings in all zones 
        clark_wcr (DataFrame)          : pandas dataframe of Clark Chart warm core rings in all zones 
        chelton_eddies (DataFrame)     : pandas dataframe of Chelton eddies in all zones
        chelton_ccr (Dataframe)        : pandas dataframe of Chelton cold core rings in all zones
        chelton_wcr (Dataframe)        : pandas dataframe of Chelton warm core rings in all zones
        title (String)                 : title of the figure
        ring_type (String)             : type of ring (or eddy) to be plotted (e.g., 'ccr', 'wcr', 'all')
        fig_quality (Integer)          : quality of the figure (e.g. 100 dpi) 
        
    Output:
        fig (Figure)                   : returns a comparative time-series of interannual variability of rings
     
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""     
    
    # 1) create figure & axes
    fig,ax = plt.subplots() # figsize=(10,5), for comparison to Fig. 2 Silver et al., 2021
    fig.set_dpi(fig_quality)
    fig.suptitle(title, y=0.935, fontsize=11);

    # 2) plot
    # CCRs or cyclones
    if (ring_type=='ccr') or (ring_type=='all'):
        ## CLARK RINGS ##
        ax.plot(clark_ccr['year'],clark_ccr['all_zones'],'-o',color='g')

        # create twin axes for chelton eddies
        ax2=ax.twinx()

        ## CHELTON EDDIES ##
        ax2.plot(chelton_eddies['year'],chelton_eddies['cyclonic'],'-o',color='k')

        ## CHELTON RINGS ##
        ax.plot(chelton_ccr['year'],chelton_ccr['all_zones'],'-o',color='orange')
    
    # WCRs or anti-cyclones
    elif (ring_type=='wcr') or (ring_type=='all'):
        ## CLARK RINGS ##
        ax.plot(clark_wcr['year'],clark_wcr['all_zones'],'-o',color='g')

        # create twin axes for chelton eddies
        ax2=ax.twinx()

        ## CHELTON EDDIES ##
        ax2.plot(chelton_eddies['year'],chelton_eddies['anticyclonic'],'-o',color='k')

        ## CHELTON RINGS ##
        ax.plot(chelton_wcr['year'],chelton_wcr['all_zones'],'-o',color='orange')


    # 3) axes formatting
    ax.set_xlabel('Years',fontweight='bold')
    ax.set_ylabel('Number of Ring Formations',fontweight='bold');
    ax2.set_ylabel('Number of Eddy Formations', fontweight='bold')
#     ax.set_ylim(5,45)
#     ax.set_xlim(1980,2020) # to compare to fig. 2 Silver et al., 2021

    # 4) custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='g', lw=4), Line2D([0], [0], color='orange', lw=4)]
    ax.legend(custom_lines, ['Chelton Eddies', 'Clark Rings', 'Chelton Rings'], loc='upper left');
    
    return fig,ax



#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#------------------------------------------------------------------------------------------------------------------------------
# )

# def ():
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     Input:

        
#     Output:
    
     
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


