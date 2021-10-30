# UWSCoordSearch
Match object/coordinates to UW field + see observation status.

Quick scripts that can be run as is or can be aliased for easy access.



*UWCoordSearch.py* takes coordinates (i.e. python UWCoordSearch.py '10:00:30.03 +02:08:59.47') and returns the matching UW field. 

*UWNameSearch.py* takes an object name (i.e. python UWNameSearch.py NGC1052) and returns the matching UW field.

Both UWCoordSearch and UWNameSearch return the observation status of the relevant field. *UWObsStatus.py* bypasses the search element and checks the status of a single field (i.e. python UWObsStatus.py UW208), returning the number and date of observations at the time of last update.

*UWStatusPlot.py* makes it easy to view and save a map of all of the UWS fields and their observation status. It takes no arguments (i.e., python UWStatusPlot.py) and saves the map to your working directory.


Dependencies:
  - PyGithub
  - Astropy
  - Pandas
  - Datetime
  - Sys
  - Matplotlib
  - Numpy




