# UWSCoordSearch
Match object/coordinates to UW field + see observation status.



Crossmatching of catalogs (just the NASA-Sloan Atlas and the Updated Nearby Galaxy Catalog for now -- please open an issue if you believe there are additional catalogs that would be suited to crossmatching with the survey) is available in `uwstatus.plotting.field()`.

If you use this package or the other scripts in this repository (described further below) in a publication or presentation, please add a footnote linking to _https://github.com/avapolzin/UWSStatusSearch_ and consider adding this software to your acknowledgments. If you do use the single field crossmatching function, please cite [Karachentsev, Makarov, & Kaisina (2013)](https://ui.adsabs.harvard.edu/abs/2013AJ....145..101K/abstract) and reference the [NASA-Sloan Atlas](http://nsatlas.org/data). If you use the functionality to view overlapping virial radii from crossmatched sources, please also cite [Mowla et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...872L..13M/abstract).

***
## Inside _UWSStatusSearch/cli_
Quick scripts that can be run as is, aliased, or made executable and put in the path for easy access.

Not all functionality is available with the scripts in this directory, since they're based on an older version of these tools. For instance, the only plotting option is to show the status of the whole survey. To cobble together the other plotting functionality, you can see the scripts that inspired that functionality at [here](https://github.com/avapolzin/DFUltrawide).

***UWCoordSearch.py*** takes coordinates (i.e. python UWCoordSearch.py '10:00:30.03 +02:08:59.47') and returns the matching UW field. 

***UWNameSearch.py*** takes an object name (i.e. python UWNameSearch.py NGC1052) and returns the matching UW field.

Both UWCoordSearch and UWNameSearch return the observation status of the relevant field. ***UWObsStatus.py*** bypasses the search element and checks the status of a single field (i.e. python UWObsStatus.py UW208), returning the number and date of observations at the time of last update.

***UWStatusPlot.py*** makes it easy to view and save a map of all of the UWS fields and their observation status. It takes no arguments (i.e., python UWStatusPlot.py) and saves the map to your working directory.


Dependencies:
  - PyGithub
  - Astropy
  - Pandas
  - Datetime
  - Sys
  - Matplotlib
  - Numpy




