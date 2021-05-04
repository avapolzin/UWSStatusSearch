import pandas as pd
import datetime
import sys
from github import Github

search = str(sys.argv[1])

UWobs = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSCoordSearch/main/UWobs.txt', sep = '\t')

result = UWobs.loc[UWobs.ID == search]
if len(result) == 0:
	raise Exception('Not a valid UW field name. Field IDs are listed UWXXX.')
field = result.ID.values[0]
Nframes = result.Nobs.values[0]
obs_dates = result.obs_dates.values[0]

if Nframes >= 100:
	status = 'been fully observed'
if Nframes == 0:
	status = 'not been observed'
if (Nframes < 100) and (Nframes > 0): 
	status = 'been partially observed'

g = Github()
repo = g.get_repo('avapolzin/UWSCoordSearch')
commits = repo.get_commits(path='UWobs.txt')

last_update = commits[0].commit.committer.date.strftime('%b %d %Y')

print('\n\n')
print('Observation status last updated %s:\n'%last_update)
print('%s has %s, with %i frames. \n'%(field, status, Nframes))
print('Nights observed: %s'%(obs_dates))
print('\n')
