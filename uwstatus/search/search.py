import pandas as pd
import datetime
from github import Github
from astropy.coordinates import SkyCoord

### Dragonfly UW observing status ###
UWobs = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/UWobs.txt', sep = '\t')
coords = SkyCoord(UWobs.RA, UWobs.Dec, unit = ['hr', 'deg'])

##########
# * * * * 
##########

def status(search, type = 'field', verbose = True):

	"""
	Relays the observing status of a UW field, including dates observed.

	Parameters:
		search (str): Depending on the type, the field name (type "field"), 
			object name (type "name"), or object coordinates (type "coord").
			Default is "field". If type = "coord", RA must be in hours and Dec 
			in degrees.
		type (str): Options are "field", "name", and "coord", which defines the 
			what the search string is.
		verbose (bool): Allows or suppresses printed information. Default is True. 

	Returns:
		Nframes (int): Number of frames taken for the relevant field.
		obs_dates (list, str): Lists the dates on which the relevant field was observed.
	"""

	if type == 'field':
		result = UWobs.loc[UWobs.ID == search]
		if len(result) == 0:
			raise Exception('Not a valid UW field name. Field IDs are listed UWXXX.')
		field = result.ID.values[0]
		Nframes = result.Nobs.values[0]
		obs_dates = result.obs_dates.values[0]
		field_center = SkyCoord(result.RA, result.Dec, unit = ['hr', 'deg'])

		if verbose:

			if Nframes >= 100:
				status = 'been fully observed'
			if Nframes == 0:
				status = 'not been observed'
			if (Nframes < 100) and (Nframes > 0): 
				status = 'been partially observed'

			g = Github()
			repo = g.get_repo('avapolzin/UWSStatusSearch')
			commits = repo.get_commits(path='UWobs.txt')

			last_update = commits[0].commit.committer.date.strftime('%b %d %Y')

			print('\n\n')
			print('%s is centered at %s (in deg: %s).'%(search, str(field_center.to_string('hmsdms')).strip('[]\'\''), str(field_center.to_string('decimal')).strip('[]\'\'')))
			print('\n\n')
			print('Observation status, last updated %s:\n'%last_update)
			print('%s has %s (%i frames). \n'%(field, status, Nframes))
			print('Nights observed: %s'%(obs_dates))
			print('\n')

	if type == 'name':
		c = SkyCoord.from_name(search)

		result = UWobs.loc[(abs(coords.ra.deg - c.ra.deg) <= 1.5) & (abs(coords.dec.deg - c.dec.deg) <= 0.9)]
		if len(result) ==0:
			raise Exception('Search falls outside of the UWS footprint.')

		field = result.ID.values[0]
		Nframes = result.Nobs.values[0]
		obs_dates = result.obs_dates.values[0]

		if verbose:

			print('\n\n')
			print('%s is located in %s'%(search, field))
			print('\n')

			if Nframes >= 100:
				status = 'been fully observed'
			if Nframes == 0:
				status = 'not been observed'
			if (Nframes < 100) and (Nframes > 0): 
				status = 'been partially observed'

			g = Github()
			repo = g.get_repo('avapolzin/UWSStatusSearch')
			commits = repo.get_commits(path='UWobs.txt')

			last_update = commits[0].commit.committer.date.strftime('%b %d %Y')


			print('Observation status, last updated %s:\n'%last_update)
			print('%s has %s, with %i frames. \n'%(field, status, Nframes))
			print('Nights observed: %s'%(obs_dates))
			print('\n')

	if type == 'coord':
		c = SkyCoord(search, unit = ['hr', 'deg'])

		result = UWobs.loc[(abs(coords.ra.deg - c.ra.deg) <= 1.5) & (abs(coords.dec.deg - c.dec.deg) <= 0.9)]
		if len(result) ==0:
			raise Exception('Search falls outside of the UWS footprint.')

		field = result.ID.values[0]
		Nframes = result.Nobs.values[0]
		obs_dates = result.obs_dates.values[0]

		if verbose:

			print('\n\n')
			print('%s is located in %s'%(search, field))
			print('\n')

			if Nframes >= 100:
				status = 'been fully observed'
			if Nframes == 0:
				status = 'not been observed'
			if (Nframes < 100) and (Nframes > 0): 
				status = 'been partially observed'

			g = Github()
			repo = g.get_repo('avapolzin/UWSStatusSearch')
			commits = repo.get_commits(path='UWobs.txt')

			last_update = commits[0].commit.committer.date.strftime('%b %d %Y')


			print('Observation status, last updated %s:\n'%last_update)
			print('%s has %s, with %i frames. \n'%(field, status, Nframes))
			print('Nights observed: %s'%(obs_dates))
			print('\n')

	return Nframes, obs_dates


def namesearch(search, verbose = True):

	"""
	Matches object to UW field and relays observing status of that field.

	Parameters:
		search (str): The common object name. Used to search for coordinates, 
			which are matched to a UW field.
		verbose (bool): Allows or suppresses printed information. Default is True. 

	Returns:
		field (str): Name of UW field object falls into.
		Nframes (int): Number of frames taken for the relevant field.
		obs_dates (list, str): Lists the dates on which the relevant field was observed.
	"""

	c = SkyCoord.from_name(search)

	result = UWobs.loc[(abs(coords.ra.deg - c.ra.deg) <= 1.5) & (abs(coords.dec.deg - c.dec.deg) <= 0.9)]
	if len(result) ==0:
		raise Exception('Search falls outside of the UWS footprint.')

	field = result.ID.values[0]
	Nframes = result.Nobs.values[0]
	obs_dates = result.obs_dates.values[0]

	if verbose:

		print('\n\n')
		print('%s is located in %s'%(search, field))
		print('\n')

		if Nframes >= 100:
			status = 'been fully observed'
		if Nframes == 0:
			status = 'not been observed'
		if (Nframes < 100) and (Nframes > 0): 
			status = 'been partially observed'

		g = Github()
		repo = g.get_repo('avapolzin/UWSStatusSearch')
		commits = repo.get_commits(path='UWobs.txt')

		last_update = commits[0].commit.committer.date.strftime('%b %d %Y')


		print('Observation status, last updated %s:\n'%last_update)
		print('%s has %s, with %i frames. \n'%(field, status, Nframes))
		print('Nights observed: %s'%(obs_dates))
		print('\n')


	return field, Nframes, obs_dates

def coordsearch(search, verbose = True):

	"""
	Matches coordinates to a UW field and relays the observing status of that field.

	Parameters:
		search (str): Coordinates with RA in hours and Dec in degrees.
		verbose (bool): Allows or suppresses printed information. Default is True. 

	Returns:
		field (str): Name of UW field coordinates fall into.
		Nframes (int): Number of frames taken for the relevant field.
		obs_dates (list, str): Lists the dates on which the relevant field was observed.
	"""

	c = SkyCoord(search, unit = ['hr', 'deg'])

	result = UWobs.loc[(abs(coords.ra.deg - c.ra.deg) <= 1.5) & (abs(coords.dec.deg - c.dec.deg) <= 0.9)]
	if len(result) ==0:
		raise Exception('Search falls outside of the UWS footprint.')

	field = result.ID.values[0]
	Nframes = result.Nobs.values[0]
	obs_dates = result.obs_dates.values[0]

	if verbose:

		print('\n\n')
		print('%s is located in %s'%(search, field))
		print('\n')

		if Nframes >= 100:
			status = 'been fully observed'
		if Nframes == 0:
			status = 'not been observed'
		if (Nframes < 100) and (Nframes > 0): 
			status = 'been partially observed'

		g = Github()
		repo = g.get_repo('avapolzin/UWSStatusSearch')
		commits = repo.get_commits(path='UWobs.txt')

		last_update = commits[0].commit.committer.date.strftime('%b %d %Y')


		print('Observation status, last updated %s:\n'%last_update)
		print('%s has %s, with %i frames. \n'%(field, status, Nframes))
		print('Nights observed: %s'%(obs_dates))
		print('\n')

	return field, Nframes, obs_dates


