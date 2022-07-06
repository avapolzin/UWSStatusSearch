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


### Dragonfly UW observing status ###
UWobs = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/UWobs.txt', sep = '\t')
coords = SkyCoord(UWobs.RA, UWobs.Dec, unit = ['hr', 'deg'])


### Mowla et al. (2019) -- https://ui.adsabs.harvard.edu/abs/2019ApJ...872L..13M/abstract ###
def get_Rad(mass):
	r_p = 8.6*u.kpc
	M_p = 10**10.2 * u.M_sun
	a = 0.17
	b = 0.50
	d = 6   
	
	r_80 = r_p * (mass/M_p)**a * ((1/2)*(1 + (mass/M_p)**d))**((b-a)/d)

	r_vir = r_80/0.047
	return r_vir


### Updated Nearby Galaxies Catalog -- https://ui.adsabs.harvard.edu/abs/2013AJ....145..101K/abstract ###
LV_ = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/LV_nearbygalaxiescatalog.txt', header = 0, sep = '\s+', comment = '#')
# will assume a mass-to-light of 1 for K-band
LV_coords = SkyCoord(LV_['RA'].values, LV_['Dec'].values, unit = ['hr', 'deg'])

### adapted from the NASA-Sloan Atlas -- http://nsatlas.org ###


h=0.7
z_conv = (c.c/(70*(u.km/u.s)/u.Mpc)).decompose().to(u.Mpc).value

##########
# * * * * 
##########


def survey(highlight_fields = False, save_path = False, **kwargs):

	"""
	Plot overall survey status.

	Parameters:
		highlight_fields (list, str): Default is False. If not False, list of fields to 
			mark in plot of survey observing status.
		save_path (str, path): Default is False, in which case figure is not saved. 
			Otherwise, feed save a path to save the plot.
		**kwargs: Figure keyword arguments to feed to matplotlib.

	Returns:
		fig (obj): The figure object.
		ax (obj): The subplot axes object. 

	"""
	g = Github()
	repo = g.get_repo('avapolzin/UWSStatusSearch')
	commits = repo.get_commits(path='UWobs.txt')
	last_update = commits[0].commit.committer.date.strftime('%b %d %Y')

	fig= plt.figure(figsize=(((max(UWcoords.ra.deg) - min(UWcoords.ra.deg))/(max(UWcoords.dec.deg) - min(UWcoords.dec.deg)))*10, 15), **kwargs)
	ax = plt.subplot(projection='mollweide')
	ax.grid(color='gray', ls='dashed')
	ax.scatter(UWcoords.ra.radian[UWobs.Nobs.values == 0] - np.deg2rad(180), UWcoords.dec.radian[UWobs.Nobs.values == 0], 
		s = 130., c='lavenderblush', label = r'N$_{obs}$ = 0')
	ax.scatter(UWcoords.ra.radian[np.logical_and((UWobs.Nobs.values > 0),(UWobs.Nobs.values < 100))] - np.deg2rad(180), 
	           UWcoords.dec.radian[np.logical_and((UWobs.Nobs.values > 0),(UWobs.Nobs.values < 100))], s = 130., c='goldenrod', 
	           label = r'0 < N$_{obs}$ < 100')
	ax.scatter(UWcoords.ra.radian[UWobs.Nobs.values >= 100] - np.deg2rad(180), UWcoords.dec.radian[UWobs.Nobs.values >= 100], 
		s = 130., c='seagreen', label = r'N$_{obs}$ >= 100')

	ax.text(np.deg2rad(230), np.deg2rad(-75), 
		s = 'Survey is %i%% completed as of %s.'%((len(UWobs[UWobs.Nobs >= 100])/len(UWobs))*100, last_update), fontsize = 12)

	if highlight_fields:
		for i in highlight_fields:
			highlight_ra = UWcoords.ra.radian[UWobs.ID == i] - np.deg2rad(180)
			highlight_dec = UWcoords.dec.radian[UWobs.ID == i]
			ax.scatter(highlight_ra, highlight_dec, marker = 'x', s = 130, c = 'k')
			ax.text(highlight_ra - np.deg2rad(1), highlight_dec + np.deg2rad(5), 
				s = '%s'%i, color = 'k', fontsize = 15)

	tick_labels = np.array([30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
	ax.set_xticklabels(tick_labels) 
	ax.legend(loc = 'best')
	ax.set_xlabel('Right Ascension (deg)')
	ax.set_ylabel('Declination (deg)')

	plt.show()
	if save_path:
		fig.savefig(save_path + 'UWS_observed_fields_%s.png'%last_update.replace(' ', ''), 
			bbox_inches = 'tight')

	return fig, ax


def field(field, highlight_coords = False, save_path = False, virial = False, 
	dist_limits = [0, 20], verbose = True, **kwargs):

	"""
	Plot single UW survey field, including crossmatch of previously ID'd objects from 
	other surveys.

	Parameters:
		field (str): UW field to plot.
		highlight_coords (list str): Default is False. If not False, list of coordinates to 
			mark in plot of survey observing status. RA should be in hours 
			and Dec in degrees.
		save_path (str, path): Default is False, in which case figure is not saved. 
			Otherwise, feed save a path to save the plot.
		virial (bool): Default is False. If True, will plot virial radii of crossmatched 
			galaxies (within distance limits) that overlap with the field. Will also read
			out percentage of field covered by virial radii.
		dist_limits (list, range): Default is [0, 20], which limits distance for 
			crossmatched galaxies to those < 20 Mpc. The format is [lower limit, upper limit]
			in Mpc. Because of the catalogs in use, the maximum distance that can be 
			accomodated is 
		verbose (bool): Default is True. If True, will print names, locations, inferred 
			stellar masses, and distances of crossmatched galaxies within the field.
		**kwargs: Figure keyword arguments to feed to matplotlib.

	Returns:
		fig (obj): The figure object.
		ax (obj): The subplot axes object. 

	"""










	return fig, ax

