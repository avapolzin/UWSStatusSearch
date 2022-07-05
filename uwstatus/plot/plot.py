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


UWobs = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/UWobs.txt', sep = '\t')
coords = SkyCoord(UWobs.RA, UWobs.Dec, unit = ['hr', 'deg'])

def survey(highlight_fields = False, save_path = False, **kwargs):

	"""
	Plot overall survey status.

	Parameters:
		highlight_fields (list, str): Default is False. If not none, list of fields to 
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
		ax.scatter()

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
	verbose = False, ):





