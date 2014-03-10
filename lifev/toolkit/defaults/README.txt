This folder contains default solver settings for LifeV Cardiac Toolkit.

DO NOT MODIFY THE FILES IN THIS FOLDER UNLESS YOU ARE SURE WHAT YOU ARE DOING.

If you wish to modify the settings of your particular solver, it is better to copy the particular solver setting to your local and provide it to the executable
as a command line parameter. For example, the local file MyMonoDomainSettings.xml would be passed as follows:

  lvelectrom -m MyModel -f MyMonoDomainSettings.xml -o OutputFolder

where the 
