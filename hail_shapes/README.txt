This software package enables the creation of MESH swath shapefiles for three specified thresholds on a daily, realtime basis, as well as for offline analyses.

###### SETUP ######
The main driver script is 'dailyHailSwath.py'.
'logging_def.py' and 'shapefile.py' are supplemental libraries. Keep them all in the same directory.
This software has been tested in python2.7 (miniconda distribution -- https://conda.io/miniconda.html)

You may also need some other libraries, which should be readily downloaded and installed using the 'conda' command in miniconda:
- shapely
- numpy
- pygrib
- netCDF4 (only need for testing purposes)
See https://conda.io/miniconda.html for installing the miniconda distribution and these extra libraries.

###### EXAMPLES ######
Example for realtime:
[user@machine]$ python dailyHailSwath.py -o RTshapes -l RTshapes/rt.log -r

Example for offline:
[user@machine]$ python dailyHailSwath.py -fp /data/MRMS_MESH_Max_1440min_*_20180413_12* -o shapes 

Example for realtime, non-default thresholds, and non-default 12Z->12Z accumulation (this will contour 1", 2", and 3" MESH for a 15Z->15Z accumulation):
[user@machine]$ python dailyHailSwath.py -o RTshapes -t 1.0 2.0 3.0 -hr 15

Example for offline simplification of shapes:
[user@machine]$ python dailyHailSwath.py -fp /data/MRMS_MESH_Max_1440min_*_20180413_12* -o shapes -s 0.01

###### FULL HELP ######
To see the full help, run from the command line:
[user@machine ~]$ python dailyHailSwath.py -h
usage: dailyHailSwath.py [-h] [-f FILE] [-fp FILEPATTERN]
                         [-t THRESHOLDS THRESHOLDS THRESHOLDS] [-o OUTDIR]
                         [-hr HOUR] [-r] [-s SIMPLIFY] [-l LOGFILE]
                         [-rn READ_NETCDF]

This program will make shapefiles for daily MRMS MESH tracks in realtime, or
create shapefiles for MRMS MESH tracks for custom timeframes in offline mode.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  [optional] Full path and filename to MRMS grib2 file,
                        for realtime. Default = 'http://mrms.ncep.noaa.gov/dat
                        a/2D/MESH/MRMS_MESH_Max_1440min.latest.grib2.gz
  -fp FILEPATTERN, --filepattern FILEPATTERN
                        [optional] Full filepattern to process local archived
                        files. E.g., '/home/user/MRMS_MESH_Max_1440min_*_20180
                        310_2000*grib2', or
                        '/home/user/MRMS_MESH_*_20180310-*grib2'. In this way,
                        whether you have a 1440min accumulation or MESH grib2s
                        every 2 min you can create the shapefiles.
  -t THRESHOLDS THRESHOLDS THRESHOLDS, --thresholds THRESHOLDS THRESHOLDS THRESHOLDS
                        [optional] Thresholds to contour MESH. Default = 0.5,
                        1.0, 2.5 (inches)
  -o OUTDIR, --outdir OUTDIR
                        [optional] Output directory. Default = ${PWD}
  -hr HOUR, --hour HOUR
                        [optional] For realtime mode, provide the UTC hour
                        that you want the 24-hr accumulation (e.g., 06).
                        Default = 12
  -r, --realtime        [optional] Run in realtime mode. This means it will
                        run indefinitely. If not set, you should specify a
                        filepattern with -fp.
  -x, --execute_once    [optional] Execute once, but for realtime data. This
                        can be used in conjunction with crontab.
  -s SIMPLIFY, --simplify SIMPLIFY
                        [optional] Use this simplification factor to reduce
                        the number of vertices. Try 0.0025 first, then
                        increase the factor if you need to. You should not go
                        over a factor of 0.5.
  -l LOGFILE, --logfile LOGFILE
                        [optional] Logfile to write output. Will always write
                        to screen regardless.
  -rn READ_NETCDF, --read_netcdf READ_NETCDF
                        [optional] Read in this netcdf file (an accumulated
                        grid). This is for testing.

Example for realtime: [user@machine]$ python dailyHailSwath.py -o RTshapes -l
RTshapes/rt.log -r

