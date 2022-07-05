import numpy as numpy
import matplotlib.pyplot as plt
import pandas as pd

import shapely.geometry as sg
import shapely.ops as so
from descartes import PolygonPatch

from astropy.coordinates import SkyCoord
import astropy.constants as c
import astropy.units as u

import matplotlib
matplotlib.rcParams.update({'font.size':20})