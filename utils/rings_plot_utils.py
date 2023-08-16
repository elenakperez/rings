"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Utilities to plot META eddy trajectories & statistical analyses

    1)  all_eddy_tracks_map : returns a map of NWA with eddy tracks for the whole period
    2)                      : 
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
# )
def all_eddy_tracks_map(eddy_wcr_df, eddy_ccr_df, bathy, title, fig_quality):
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
        ax.plot(eddy['longitude'],eddy['latitude'],color='#F42119') # anticyclonic
    
    ## CCR-like eddies ##     
    for i in np.array(eddy_ccr_df['track'].unique()):
        eddy = eddy_ccr_df[eddy_ccr_df['track']==i]
        ax.plot(eddy['longitude'],eddy['latitude'],color='#114DCE') # cyclonic

    # axes formatting
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.contour(bathy.lon,bathy.lat,bathy.z,levels=[-4000,-1000,-100],colors='gray') #,levels=[-4000,-1000,-100]
    ax.coastlines(resolution='50m',color='gray')
    ax.set_extent([x_bnds[0],x_bnds[1],y_bnds[0],y_bnds[1]], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, color='lightgray')  

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#F42119', lw=4), Line2D([0], [0], color='#114DCE', lw=4)]
    ax.legend(custom_lines, ['WCR-like', 'CCR-like'], loc='upper left')

    # gridlines
    gl = ax.gridlines(crs=proj,draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    
    return fig,ax


#-------------------------------------------------------------------------------------------------------------------------------
















