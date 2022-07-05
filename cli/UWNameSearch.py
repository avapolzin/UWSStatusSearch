import pandas as pd
from astropy.coordinates import SkyCoord
import datetime
import sys
from github import Github

search = str(sys.argv[1])

UWobs = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/UWobs.txt', sep = '\t')

coords = SkyCoord(UWobs.RA, UWobs.Dec, unit = ['hr', 'deg'])

c = SkyCoord.from_name(search)

result = UWobs.loc[(abs(coords.ra.deg - c.ra.deg) <= 1.5) & (abs(coords.dec.deg - c.dec.deg) <= 0.9)]
if len(result) ==0:
	raise Exception('Search falls outside of the UWS footprint.')

field = result.ID.values[0]
Nframes = result.Nobs.values[0]
obs_dates = result.obs_dates.values[0]

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