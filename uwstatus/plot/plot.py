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
NSA_table = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/nsa_v0_1_2_abridged.txt', header = 0, sep = '\t')
h=0.7
z_conv = (c.c/(70*(u.km/u.s)/u.Mpc)).decompose().to(u.Mpc).value

##########
# * * * * 
##########


def survey(highlight_fields = False, text_label = True, save_path = False, **kwargs):

	"""
	Plot overall survey status.

	Parameters:
		highlight_fields (list, str): Default is False. If not False, list of fields to 
			mark in plot of survey observing status.
		text_label (bool): Default is True. Label highlighted fields with 'UWXXXX' on plot.
			If highlight_fields = False, this parameter is ignored.
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
			if text_label:
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
			crossmatched galaxies to those <= 20 Mpc. The format is inclusive for 
			[lower limit, upper limit] in Mpc. Because of the catalogs in use, the 
			maximum distance that can be accomodated is ~239 Mpc.
		verbose (bool): Default is True. If True, will print names, locations, inferred 
			stellar masses, and distances of crossmatched galaxies within the field.
		**kwargs: Figure keyword arguments to feed to matplotlib.

	Returns:
		fig (obj): The figure object.
		ax (obj): The subplot axes object. 

	"""
	field_ = sg.box(coords[UWobs.ID == field].ra.deg - 1.5, coords[UWobs.ID == field].dec.deg - 1, 
		coords[UWobs.ID == field].ra.deg + 1.5, coords[UWobs.ID == field].dec.deg + 1)
	virial_ = []

	deg_cut = 20 #radius to include in the virial calc
	
	## Nearby Galaxies Catalog
	RA_LV = LV_coords.ra.deg[np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= deg_cut), 
		(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= deg_cut)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
	Dec_LV = LV_coords.dec.deg[np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= deg_cut), 
		(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= deg_cut)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
	dist_LV = LV_['dist_Mpc'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= deg_cut), 
		(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= deg_cut)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
	loglum_LV = LV_['loglumK'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= deg_cut), 
		(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= deg_cut)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
	loglum_LV.reset_index(drop = True)
	
	## NASA Sloan ATLAS
	dist_ = z_conv*NSA_table['ZDIST'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= deg_cut), 
	   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= deg_cut)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
	RA_ = NSA_table['RA'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= deg_cut), 
	   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= deg_cut)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
	Dec_ = NSA_table['DEC'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= deg_cut), 
	   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= deg_cut)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
	
	masses_ = h**2 * NSA_table['MASS'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= deg_cut), 
	   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= deg_cut)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
	
	fig, ax = plt.subplots(1, 1, figsize = (12, 8))

	if virial:
		if len(dist_) > 0:
			for j in range(len(dist_)):
				ang_rad = (get_Rad(masses_[j]*u.Msun)/(dist_[j]*u.Mpc)).decompose()*u.rad.to(u.deg)
				# test = sg.Point(RA_[j], Dec_[j]).buffer(ang_rad)
				virial_.append(sg.Point(RA_[j], Dec_[j]).buffer(ang_rad))

		if len(loglum_LV) > 0:
			for k in range(len(loglum_LV)):
				mass_LV = 10**loglum_LV.values[k]
				ang_rad = (get_Rad(mass_LV*u.Msun)/(dist_LV.values[k]*u.Mpc)).decompose()*u.rad.to(u.deg)
				virial_.append(sg.Point(RA_LV[k], Dec_LV[k]).buffer(ang_rad))
	
		area = field.intersection(so.cascaded_union(virial_)).area
	
	
		for m in virial_:
			ax.add_patch(PolygonPatch(m, fc = 'gray', ec = 'gray', alpha = 0.6, zorder = 0))

		ax.set_title(str(UWobs.ID[UWobs.ID == field]) + ', %.2f percent covered'%(area/6 * 100))
	
	ax.scatter(RA_LV[np.isfinite(loglum_LV.values) & (loglum_LV.values >= 10)], Dec_LV[np.isfinite(loglum_LV.values) & (loglum_LV.values >= 10)], 
		color = 'k', marker = 'X', label = r'$M_* >= 10^{10} \, M_\odot$', s = 150)
	ax.scatter(RA_[(dist_ <= 20) & (masses_ >= 10**10)], Dec_[(dist_ <= 20) & (masses_ >= 10**10)], color = 'red', marker = 'X', s = 150)

	ax.scatter(RA_LV[np.isfinite(loglum_LV.values) & (loglum_LV.values <= 10) & (loglum_LV.values >= 7)], 
		Dec_LV[np.isfinite(loglum_LV.values) & (loglum_LV.values <= 10) & (loglum_LV.values >= 7)], 
		color = 'k', marker = 'x', label = r'$M_* <= 10^{10} \, M_\odot$ and $M_* >= 10^{7} \, M_\odot$', s = 150)
	ax.scatter(RA_[(masses_ <= 10**10) & (masses_ >= 10**7)], 
		Dec_[(masses_ <= 10**10) & (masses_ >= 10**7)], 
		color = 'k', marker = 'x', s = 150)

	ax.scatter(RA_LV[~np.isfinite(loglum_LV.values)], Dec_LV[~np.isfinite(loglum_LV.values)], 
		color = 'k', label = r'$M_* <= 10^{7} \, M_\odot$', s = 75)
	ax.scatter(RA_LV[np.isfinite(loglum_LV.values) & (loglum_LV.values <= 7)], Dec_LV[np.isfinite(loglum_LV.values) & (loglum_LV.values <= 7)], 
		color = 'k', s = 75)
	ax.scatter(RA_[masses_ <= 10**7], 
		Dec_[masses_ <= 10**7], 
		color = 'k', s = 75)

	ax.legend(loc = 'best', fontsize = 10)

	ax.set_xlim(coords[UWobs.ID == field].ra.deg - 1.5, coords[UWobs.ID == field].ra.deg + 1.5)
	ax.set_xlabel('RA (deg)')
	ax.set_ylim(coords[UWobs.ID == field].dec.deg - 1, coords[UWobs.ID == field].dec.deg + 1)
	ax.set_ylabel('Dec (deg)')
	ax.invert_xaxis()
	
	fig.savefig(save_path + str(UWobs.ID[UWobs.ID == field]) + '.png', bbox_inches = 'tight')
	plt.close()


	if verbose:
		## Nearby Galaxies Catalog
		name_LV = LV_['name'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= 1.5), 
			(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= 1)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
		RA_LV = LV_coords.ra.deg[np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= 1.5), 
			(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= 1)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
		Dec_LV = LV_coords.dec.deg[np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= 1.5), 
			(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= 1)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
		dist_LV = LV_['dist_Mpc'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= 1.5), 
			(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= 1)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
		loglum_LV = LV_['loglumK'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - LV_coords.ra.deg) <= 1.5), 
			(abs(coords[UWobs.ID == field].dec.deg - LV_coords.dec.deg) <= 1)) & (LV_['dist_Mpc'] >= dist_limits[0]) & (LV_['dist_Mpc'] <= dist_limits[1])]
		loglum_LV.reset_index(drop = True)
		masses_LV = 10**loglum_LV.values #assumes K-band M/L = 1
		
		## NASA Sloan ATLAS
		name_ = NSA_table['IAUNAME'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= 1.5), 
		   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= 1)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
		dist_ = z_conv*NSA_table['ZDIST'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= 1.5), 
		   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= 1)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
		RA_ = NSA_table['RA'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= 1.5), 
		   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= 1)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
		Dec_ = NSA_table['DEC'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= 1.5), 
		   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= 1)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]
		
		masses_ = h**2 * NSA_table['MASS'][np.logical_and((abs(coords[UWobs.ID == field].ra.deg - NSA_table['RA']) <= 1.5), 
		   (abs(coords[UWobs.ID == field].dec.deg - NSA_table['DEC']) <= 1)) & (z_conv*NSA_table['ZDIST'] >= dist_limits[0]) & (z_conv*NSA_table['ZDIST'] <= dist_limits[1])]

		gals = {'Name':np.concatenate((name_LV, name_)), 'RA':np.concatenate((RA_LV, RA_)), 
		'Dec':np.concatenate((Dec_LV, Dec_)), 'Dist(Mpc)':np.concatenate((dist_LV, dist_)), 
		'Mass(Msun)':np.concatenate((masses_LV, masses_))}

		df = pd.DataFrame.from_dict(gals)

		print('\n\n')
		print('The known galaxies in this field are:\n')
		print(df.to_string())
		print('\n\n')



	return fig, ax

