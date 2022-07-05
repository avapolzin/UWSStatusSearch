# UWSCoordSearch
Match object/coordinates to UW field + see observation status.



Crossmatching of catalogs (just the NASA-Sloan Atlas and the Updated Nearby Galaxy Catalog for now -- please open an issue if you believe there are additional catalogs that would be suited to crossmatching with the survey) is available in `uwstatus.plotting.field()`.

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




