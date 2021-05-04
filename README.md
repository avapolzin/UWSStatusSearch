# UWSCoordSearch
Match object/coordinates to UW field + see observation status.

Quick scripts that can be run as is or can be aliased for easy access.

UWObsStatus.py takes a single field as an argument (i.e. UWObsStatus UW208) and returns the number and date of observations at the time of last update.

UWCoordSearch.py takes coordinates (i.e. 10:00:30.03 +02:08:59.47) or an object name and a flag (i.e. -n NGC1052) and returns the matching UW field. The verbose flag (-v) will also return the number and date of observations.
