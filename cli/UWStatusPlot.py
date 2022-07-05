import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
import pandas as pd
from astropy.coordinates import SkyCoord
import numpy as np
import datetime
import sys
from github import Github

UWobs = pd.read_csv('https://raw.githubusercontent.com/avapolzin/UWSStatusSearch/main/UWobs.txt', sep = '\t')

UWcoords = SkyCoord(UWobs.RA, UWobs.Dec, unit = ['hr', 'deg'])

g = Github()
repo = g.get_repo('avapolzin/UWSStatusSearch')
commits = repo.get_commits(path='UWobs.txt')
last_update = commits[0].commit.committer.date.strftime('%b %d %Y')

fig= plt.figure(figsize=(((max(UWcoords.ra.deg) - min(UWcoords.ra.deg))/(max(UWcoords.dec.deg) - min(UWcoords.dec.deg)))*10, 15))
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

tick_labels = np.array([30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
ax.set_xticklabels(tick_labels) 
ax.legend(loc = 'best')
ax.set_xlabel('Right Ascension (deg)')
ax.set_ylabel('Declination (deg)')

plt.show()
fig.savefig('UWS_observed_fields_%s.png'%last_update.replace(' ', ''), bbox_inches = 'tight')
