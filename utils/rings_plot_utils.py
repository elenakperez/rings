"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Utilities to plot META eddy trajectories & statistical analyses

    1)  : 
    2)  : 
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