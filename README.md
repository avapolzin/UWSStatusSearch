# UWSCoordSearch
Match object/coordinates to UW field + see observation status.

Quick scripts that can be run as is or can be aliased for easy access.



UWCoordSearch.py takes coordinates (i.e. python UWCoordSearch.py '10:00:30.03 +02:08:59.47') and returns the matching UW field. 

UWNameSearch.py takes an object name (i.e. python UWNameSearch.py NGC1052) and returns the matching UW field.

Both UWCoordSearch and UWNameSearch return the observation status of the relevant field. UWObsStatus.py bypasses the search element and checks the status of a single field (i.e. python UWObsStatus.py UW208), returning the number and date of observations at the time of last update.


Dependencies:
  - PyGithub
  - Astropy
  - Pandas
  - Datetime
  - Sys


![UWS_observed_fields](https://user-images.githubusercontent.com/29441772/124369043-848ad300-dc35-11eb-9881-c00e384337d2.png)


![UWS_observed+reduced_fields](https://user-images.githubusercontent.com/29441772/124369045-86549680-dc35-11eb-9944-5d1ce44d7028.png)



