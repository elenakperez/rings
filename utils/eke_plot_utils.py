"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Utilities to plot CMEMS geostrophic velocity data
 
    1) plot_region_masks_example           : returns a figure of the NWA's regional masks and GS path for Jan. 1993
    2) plot_regionalArea                   : returns a figure of area for a specified region (Slope, Gulf Stream, or Sargasso)
    3) plot_monthly_eke_ts                 : returns a timeseries of monthly EKE for regions (Slope, GS, Sargasso, NWA)
    4) plot_annual_eke_ts                  : returns a timeseries of annual EKE for regions (Slope, GS, Sargasso, NWA)
    5) plot_speed                          : returns a spatial plot of the speed of NWA for a given period (pre, post, or diff)
    6) plot_ke                             : returns a spatial plot of the KE of NWA for a given period (pre, post, or diff)
    7) plot_gs_eke_wcr_ts                  : returns a dual timeseries plot of EKE vs. WCR formations
    

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

# import the util functions
from utils.ring_data_utils import * # for plotting Gulf Stream paths

# turn off warnings
import warnings
warnings.filterwarnings("ignore")


#-------------------------------------------------------------------------------------------------------------------------------
# 1) 
def plot_region_masks_example(geo_vels, bathy_nwa, fig_quality, title):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    geo_vels (DataArray)       : xarray DataArray of geostrophic velocities and masks                  "
    "    bathy_nwa (DataArray)      : bathymetry file for plotting                                          "
    "    fig_quality (Int)          : desired quality of the figure                                         "
    "    title (String)             : desired title of the figure                                           "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a map of the regional masks and Gulf Stream position for Jan. 1993                       "
    "                                                                                                       "
    "                                                                                                       "
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    
    x_bnds = [-85,-55]
    y_bnds = [30,45]

    # map projection
    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)
    # fig.suptitle(title, fontsize=20, y=0.925)

    # uncomment to plot an example regional mask
    (geo_vels["mask_slope"][0]).plot(add_colorbar=False, cmap='winter'); # Slope Sea mask for Jan. 1993
    (geo_vels["mask_gs"][0]).plot(add_colorbar=False, cmap='summer', label='Gulf Stream'); # Gulf Stream mask for Jan. 1993
    (geo_vels["mask_sag"][0]).plot(add_colorbar=False, cmap='autumn'); # Sargasso Sea mask for Jan. 1993

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy_nwa.lon,bathy_nwa.lat,bathy_nwa.z,levels=[-4000,-1000,-100],colors='k', alpha=0.8, zorder=10000, linestyles='solid') #,levels=[-4000,-1000,-100] '#5E5E5E'
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray');

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='royalblue', lw=10), Line2D([0], [0], color='mediumseagreen', lw=10), Line2D([0], [0], color='darkorange', lw=10)] 
    ax.legend(custom_lines, ['Slope Sea', 'Gulf Stream', 'Sargasso Sea'], loc='upper left', fontsize = 17).set_zorder(10001)

    ax.plot(get_gs_month(1993,1)[0],get_gs_month(1993,1)[1], color='k',linewidth=4);
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------
# 2)
def plot_regionalArea(regionArea_monthly, regionArea_annual, region, fig_quality, c_monthly, c_annual):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    regionArea_monthly (DataArray)   : xarray DataArray of monthly areas of region                     "
    "    regionArea_annual (DataArray)    : xarray DataArray of annual areas of region                      "
    "    region (String)                  : string of desired region to plot                                "
    "    fig_quality (Int)                : desired quality of the figure                                   "
    "    c_monthly (String)               : desired color for monthly timeseries                            " 
    "    c_annual (String)                : desired color for annual timeseries                             " 
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a map of the regional masks and Gulf Stream position for Jan. 1993                       "
    "                                                                                                       "
    "                                                                                                       "
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" 
    
    title='Area of the '+region
    
    # create figure
    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle(title, y=0.935, fontsize=11);
    
    # plot
    ax.plot(regionArea_monthly['time'],regionArea_monthly, color=c_monthly) # monthly
    ax.plot(regionArea_annual['time'], regionArea_annual, color=c_annual, linewidth=2) # annual 
    
    # axes formatting
    ax.set_xlabel('Time', fontweight='bold')
    ax.set_ylabel('Area [m^2]',fontweight='bold');
    
    return fig,ax

#-------------------------------------------------------------------------------------------------------------------------------
# 3)
def plot_monthly_eke_ts(nwa_eke, slope_eke, gs_eke, sag_eke, title, fig_quality, ymin, ymax):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    nwa_eke (DataArray)       : xarray DataArray NWA monthly EKE values                                "
    "    slope_eke (DataArray)     : xarray DataArray Slope Sea monthly EKE values                          "
    "    gs_eke (DataArray)        : xarray DataArray Gulf Stream monthly EKE values                        "
    "    sag_eke (DataArray)       : xarray DataArray Sargasso Sea monthly EKE values                       "
    "    title (String)            : title of the plot                                                      "
    "    fig_quality (Int)         : desired quality of the figure                                          "
    "    ymin (Int)                : set y-axis minimum value                                               "
    "    ymax (Int)                : set y-axis maximum value                                               "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle(title, y=0.935, fontsize=11);

    # plot time series
    ax.plot(slope_eke.time, slope_eke, label='Slope', color='#abd9e9'); 
    ax.plot(gs_eke.time, gs_eke, label='GS', color='#d7191c');
    ax.plot(sag_eke.time, sag_eke, label='Sargasso',  color='#fdae61');
    ax.plot(nwa_eke.time, nwa_eke, label='NWA', color='#2c7bb6')


    # axes formatting
    ax.set_ylabel(r"Monthly EKE $\left[\frac{J}{m^3}\right]$")
    ax.set_ylim([ymin,ymax])
    ax.set_xlabel('Time');
    ax.legend(loc='upper left', ncol=4);

    ax.set_xlim([slope_eke.time[0],slope_eke.time[-1]]);

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#2c7bb6', lw=4), Line2D([0], [0], color='#abd9e9', lw=4), Line2D([0], [0], color='#d7191c', lw=4), Line2D([0], [0], color='#fdae61', lw=4)]
    ax.legend(custom_lines, ['NWA', 'Slope', 'GS', 'Sargasso'], loc='upper left', ncols=4);
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# 4)
def plot_annual_eke_ts(nwa_eke, slope_eke, gs_eke, sag_eke, title, fig_quality, ymin, ymax):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    nwa_eke (DataArray)       : xarray DataArray NWA annual EKE values                                 "
    "    slope_eke (DataArray)     : xarray DataArray Slope Sea annual EKE values                           "
    "    gs_eke (DataArray)        : xarray DataArray Gulf Stream annual EKE values                         "
    "    sag_eke (DataArray)       : xarray DataArray Sargasso Sea annual EKE values                        "
    "    title (String)            : title of the plot                                                      "
    "    fig_quality (Int)         : desired quality of the figure                                          "
    "    ymin (Int)                : set y-axis minimum value                                               "
    "    ymax (Int)                : set y-axis maximum valu                                                "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # create figure
    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)
    fig.suptitle(title, y=0.935, fontsize=11);

    years = np.arange(1994,2018,1)

    # plot time series
    ax.plot(slope_eke.time, slope_eke, label='Slope', color='#abd9e9', linewidth=2); 
    ax.plot(gs_eke.time, gs_eke, label='GS', color='#d7191c', linewidth=2);
    ax.plot(sag_eke.time, sag_eke, label='Sargasso',  color='#fdae61', linewidth=2);
    ax.plot(nwa_eke.time, nwa_eke, label='NWA', color='#2c7bb6', linewidth=2)

    # axes formatting
    ax.set_ylabel(r"Annual EKE $\left[\frac{J}{m^3}\right]$")
    ax.set_ylim([ymin,ymax])
    ax.set_xlabel('Time');
    ax.set_xlim([nwa_eke.time[0],nwa_eke.time[-1]]);
    ax.legend(loc='upper left', ncol=4);

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#2c7bb6', lw=4), Line2D([0], [0], color='#abd9e9', lw=4), Line2D([0], [0], color='#d7191c', lw=4), Line2D([0], [0], color='#fdae61', lw=4)]
    ax.legend(custom_lines, ['NWA', 'Slope', 'GS', 'Sargasso'], loc='upper left', ncols=4);
    
    return fig,ax
    
#-------------------------------------------------------------------------------------------------------------------------------
# 5)

def plot_speed(speed_pre, speed_post, period, bathy_nwa, fig_quality, vmin_val, vmax_val):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    speed_pre (DataArray)     : xarray DataArray of speed values (1993–1999) in the NWA, units [m/s]   "
    "    speed_post (DataArray)    : xarray DataArray of speed values (2000–2017) in the NWA, units [m/s]   "
    "    period (String)           : String indicating period, e.g. 'pre', 'post', or 'difference'          "
    "    bathy_nwa (DataArray)     : xarray DataArray bathymetry of Northwest Atlantic                      "
    "    fig_quality (Int)         : desired quality of the figure                                          "
    "    vmin_val (Int)            : set minimum value for color shading                                    "
    "    vmax_val (Int)            : set maximum value for color shading                                    "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # gangopadhyay census bounds
    x_bnds = [-80,-55] # lon, NWA: [-82,-48]
    y_bnds = [30,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)

    # plot speed
    if period=='pre':
        speed_pre.plot(cmap='Spectral_r', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
    elif period=='post':
        speed_post.plot(cmap='Spectral_r', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
    elif period=='diff':
        (speed_post-speed_pre).plot(cmap='bwr', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
       
    # plot Gulf Stream reference path (1993–2022)
    ax.plot(get_gs()[0],get_gs()[1],transform=proj, linewidth=3, color='k')
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='k',linestyle='-') # zone 1
    ax.plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
    ax.plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
    ax.plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-'); # zone 4

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy_nwa.lon,bathy_nwa.lat,bathy_nwa.z,levels=[-4000,-1000,-100],colors='k', transform = proj) #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    if period=='pre':
        ax.legend(['Speed (1993–1999)'], loc='upper left', fontsize=13);
    elif period=='post':
        ax.legend(['Speed (2000–2017)'], loc='upper left', fontsize=13);
    elif period=='diff':
        ax.legend(['Speed Difference'], loc='upper left', fontsize=13);

    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# 6)
def plot_ke(ke_pre, ke_post, period, bathy_nwa, fig_quality, vmin_val, vmax_val):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ke_pre (DataArray)     : xarray DataArray of KE values (1993–1999) in the NWA, units [J/m^3]       "
    "    ke_post (DataArray)    : xarray DataArray of KE values (2000–2017) in the NWA, units [J/m^3]       "
    "    period (String)           : String indicating period, e.g. 'pre', 'post', or 'difference'          "
    "    bathy_nwa (DataArray)     : xarray DataArray bathymetry of Northwest Atlantic                      "
    "    fig_quality (Int)         : desired quality of the figure                                          "
    "    vmin_val (Int)            : set minimum value for color shading                                    "
    "    vmax_val (Int)            : set maximum value for color shading                                    "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # gangopadhyay census bounds
    x_bnds = [-80,-55] # lon, NWA: [-82,-48]
    y_bnds = [30,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)

    # plot speed
    if period=='pre':
        ke_pre.plot(cmap='Spectral_r', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
    elif period=='post':
        ke_post.plot(cmap='Spectral_r', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
    elif period=='diff':
        (ke_post-ke_pre).plot(cmap='bwr', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
       
    # plot Gulf Stream reference path (1993–2022)
    ax.plot(get_gs()[0],get_gs()[1],transform=proj, linewidth=3, color='k')
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='k',linestyle='-') # zone 1
    ax.plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
    ax.plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
    ax.plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-'); # zone 4

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy_nwa.lon,bathy_nwa.lat,bathy_nwa.z,levels=[-4000,-1000,-100],colors='k', transform = proj) #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    if period=='pre':
        ax.legend(['KE (1993–1999)'], loc='upper left', fontsize=13);
    elif period=='post':
        ax.legend(['KE (2000–2017)'], loc='upper left', fontsize=13);
    elif period=='diff':
        ax.legend(['KE Difference'], loc='upper left', fontsize=13);

    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# 7)
def plot_eke(eke_pre, eke_post, period, bathy_nwa, fig_quality, vmin_val, vmax_val):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    eke_pre (DataArray)       : xarray DataArray of EKE values (1993–1999) in the NWA, units [J/m^3]   "
    "    eke_post (DataArray)      : xarray DataArray of EKE values (2000–2017) in the NWA, units [J/m^3]   "
    "    period (String)           : String indicating period, e.g. 'pre', 'post', or 'difference'          "
    "    bathy_nwa (DataArray)     : xarray DataArray bathymetry of Northwest Atlantic                      "
    "    fig_quality (Int)         : desired quality of the figure                                          "
    "    vmin_val (Int)            : set minimum value for color shading                                    "
    "    vmax_val (Int)            : set maximum value for color shading                                    "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # gangopadhyay census bounds
    x_bnds = [-80,-55] # lon, NWA: [-82,-48]
    y_bnds = [30,45] # lat, NWA: [24,53]

    proj = ccrs.PlateCarree()

    # create figure 
    fig,ax = plt.subplots(subplot_kw = dict(projection=proj),figsize=(10,7))
    fig.set_dpi(fig_quality)

    # plot speed
    if period=='pre':
        eke_pre.plot(cmap='Spectral_r', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
    elif period=='post':
        eke_post.plot(cmap='Spectral_r', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
    elif period=='diff':
        (eke_post-eke_pre).plot(cmap='bwr', transform = proj, ax = ax, 
                          cbar_kwargs = dict(orientation="vertical", fraction = 0.023, pad = 0.1), vmin = vmin_val, vmax = vmax_val)       
       
    # plot Gulf Stream reference path (1993–2022)
    ax.plot(get_gs()[0],get_gs()[1],transform=proj, linewidth=3, color='k')
    
    # add zone lines
    ax.plot([-75,29],[-75,45],transform=proj,linewidth=3,color='k',linestyle='-') # zone 1
    ax.plot([-70,29],[-70,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 2
    ax.plot([-65,29],[-65,45],transform=proj,linewidth=3,color='b',linestyle='-') # zone 3
    ax.plot([-60,29],[-60,45],transform=proj,linewidth=3,color='b',linestyle='-'); # zone 4

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy_nwa.lon,bathy_nwa.lat,bathy_nwa.z,levels=[-4000,-1000,-100],colors='k', transform = proj) #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    if period=='pre':
        ax.legend(['EKE (1993–1999)'], loc='upper left', fontsize=13);
    elif period=='post':
        ax.legend(['EKE (2000–2017)'], loc='upper left', fontsize=13);
    elif period=='diff':
        ax.legend(['EKE Difference'], loc='upper left', fontsize=13);

    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
# 8)
def plot_gs_eke_wcr_ts(ga, var, regionName, wcrs, fig_quality):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ga (DataArray)            : geostrophic velocity anomalies DataArray (EKE, KE, speed)              "
    "    var (String)              : which variable to plot, e.g. 'EKE_gs'                                  "
    "    regionName (String)       : region from 'var', e.g. 'Gulf Stream'                                  "
    "    wcrs (DataFrame)          : DataFrame of WCR formations                                            "
    "    fig_quality (Int)         : desired quality of the figure                                          "                                              "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    fig,ax1 = plt.subplots()
    fig.set_dpi(fig_quality)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    years = np.arange(1993,2018,1)

    # plot Gulf Stream EKE on axis 1
    ax1.plot(years, ga[var].mean(("longitude", "latitude")).resample(time = "1Y").mean(), label='GS', color = '#d7191c', linewidth=2)
    ax1.set_ylabel(r"Annual EKE $\left[\frac{J}{m^3}\right]$")
    ax1.set_xlabel('time')
    ax1.set_ylim([70,155]);

    # plot WCR formations on axis 2
    ax2.scatter(years, wcrs['all_zones'], color='#2c7bb6', linewidth=2)
    ax2.plot(years, wcrs['all_zones'], color='#2c7bb6')

    # axes formatting
    ax2.set_ylabel("Annual formations of WCRs")

    ax1.legend(loc='upper left');
    ax1.set_xlim([1993,2017]);

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#d7191c', lw=4), Line2D([0], [0], color='#2c7bb6', lw=4)]
    ax1.legend(custom_lines, [regionName,'WCRs'], loc='upper left');
    
    return fig,ax1,ax2


#-------------------------------------------------------------------------------------------------------------------------------
# 9)
def plot_leadlag_gs_wcr(ga, var, regionName, wcrs, fig_quality, ymin, ymax):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "                                                                                                       "
    " Input:                                                                                                "
    "    ga (DataArray)            : geostrophic velocity anomalies DataArray (EKE, KE, speed)              "
    "    var (String)              : which variable to plot, e.g. 'EKE_gs'                                  "
    "    regionNAme (String)       : region from 'var', e.g. 'Gulf Stream'                                  "
    "    wcrs (String)             : DataFrame of WCR formations                                            "
    "    fig_quality (Int)         : desired quality of the figure                                          "
    "    ymin (Int)                : set y-axis minimum value                                               "
    "    ymax (Int)                : set y-axis maximum value                                               "
    "                                                                                                       "
    " Output:                                                                                               "
    "    * returns a timeseries plot of monthly EKE for each region (Slope, GS, Sargasso, NWA)              "
    "                                                                                                       "
    "                                                                                                       "
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # calculations before plotting
    years = np.arange(1993,2018,1)
    def normalize(x):
        return (x - np.mean(x)) / np.std(x)

    x1 = ga[var].mean(("longitude", "latitude")).resample(time = "1Y").mean().values
    x2 = wcrs['all_zones'].values

    lags, correlations, _, _ = plt.xcorr(normalize(x1), normalize(x2));
    plt.close()
    
    # create figure
    fig,ax = plt.subplots()
    fig.set_dpi(fig_quality)

    # plot
    ax.plot(lags, correlations, color='#d7191c')
    ax.plot(lags[10:], correlations[10:], color='#2c7bb6')

    # axes formatting
    ax.set_ylabel("Correlation")
    ax.set_xlabel('Time lag (years)')
    ax.legend(loc='upper left');
    ax.set_xlim([-10,10]);
    ax.set_xticks(np.arange(-10, 11))
    ax.grid()
    ax.set_ylim([ymin,ymax])
    ax.vlines(0,-1,1,color='k')

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#d7191c', lw=4), Line2D([0], [0], color='#2c7bb6', lw=4)]
    ax.legend(custom_lines, ['Gulf Stream EKE \nleads WCR formations       ', 'Gulf Stream EKE \nlags WCR formations        '], loc='upper left', ncols=2);

    return fig, ax


#-------------------------------------------------------------------------------------------------------------------------------
# 10) example function
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











