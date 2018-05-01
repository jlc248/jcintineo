import time, sys, os
import numpy as np
import pygrib as pg
import argparse
import logging
from logging.config import dictConfig
from logging_def import logging_def
from glob import glob
from subprocess import call,Popen,PIPE
import shapefile
from datetime import datetime, timedelta
from shapely.geometry import Polygon
import shapely
#############################################################################################################
def read_data(file, var_name, atts=None, global_atts=None):
  from netCDF4 import Dataset
  FUNCTION_NAME = 'READ_DATA'
  logging.info(FUNCTION_NAME + ': Process started. Received ' + file)

  if(file[-3:] == '.gz'):
    gzip=1
    strsize = len(file)
    unzipped_file = file[0:len(file)-3]
    the_file = unzipped_file
    s = call(['gunzip','-f',file])
  else:
    gzip=0
    the_file = file

  try:
    openFile = Dataset(the_file,'r')
  except OSError as err:
    logging.error(FUNCTION_NAME + ': ' + str(err))
    return -1

  data = openFile.variables[var_name][:]
  #get variable attributes
  if(isinstance(atts,dict)):
    for att in atts:
      if(att in openFile.variables[var_name].ncattrs()):
        atts[att] = getattr(openFile.variables[var_name],att) #atts is modified and returned
  #get global attributes
  if(isinstance(global_atts,dict)):
    for gatt in global_atts:
      if(gatt in openFile.ncattrs()):
        global_atts[gatt] = getattr(openFile,gatt) #global_atts is modified and returned
  openFile.close()

  #--gzip file after closing?
  if (gzip):
    s = call(['gzip','-f',unzipped_file])

  logging.info(FUNCTION_NAME + ': Got new ' + var_name + '. Process ended.')
  return data
#############################################################################################################
def max_filter(data2D):
  FUNCTION_NAME = 'max_filter'
  logging.info(FUNCTION_NAME + ': Process started.')

  ny,nx = np.shape(data2D)
  filtered_data = np.zeros((ny,nx))

  for ii in range(1,ny-1):
    for jj in range(1,nx-1):
      filtered_data[ii,jj] = np.max(data2D[ii-1:ii+1,jj-1:jj+1])

  logging.info(FUNCTION_NAME + ': Process ended.')

  return filtered_data
#############################################################################################################
def trace_object(data2D,pixelFirst,maxsize=100000):
  FUNCTION_NAME = 'trace_object'

  #data2D is the 2D grid, with missing being < 0
  #pixelFirst is a 2D-index vector in the grid, with a starting point for an object to trace
  #***it must be on the boundary of an object***
  #returns indices in data2D for the boundary of the object

  firstPixVal = data2D[pixelFirst[0],pixelFirst[1]]
  #tempI = min(data2D[pixelFirst[0]-1:pixelFirst[0]+2,pixelFirst[1]])
  if(data2D[pixelFirst[0]-1,pixelFirst[1]] != firstPixVal or data2D[pixelFirst[0]+1,pixelFirst[1]] != firstPixVal):
    directionFirst = [1,0]
  else:
    #tempI = min(data2D[pixelFirst[0],pixelFirst[1]-1:pixelFirst[1]+2])
    if(data2D[pixelFirst[0],pixelFirst[1]-1] != firstPixVal or data2D[pixelFirst[0],pixelFirst[1]+1] != firstPixVal):
      directionFirst = [0,1]
    else:
      logging.warning(FUNCTION_NAME + ': No start pixel!')
      return -1

  counter = 0

  listCONTOUR = np.zeros((maxsize,2),dtype=int) - 999
  listCONTOUR[counter,:] = pixelFirst
  listDIRECTION = np.zeros((maxsize,2),dtype=int) - 999
  listDIRECTION[counter,:] = directionFirst
  counter += 1

  checkLOOPEND = 0

  ROTmatrix90 = [[0,1],[-1,0]]
  ROTmatrix180 = [[-1,0],[0,-1]]

  tempPIX = pixelFirst
  tempDIR = directionFirst

  while(checkLOOPEND == 0):
    tempDIR = np.dot(ROTmatrix180,tempDIR) # ROTmatrix180
    tempPIX = np.add(tempPIX,tempDIR)
    if(data2D[tempPIX[0],tempPIX[1]] >= 1):
      listCONTOUR[counter,:] = tempPIX
      listDIRECTION[counter,:] = tempDIR
      counter += 1
    else:
      tempDIR = np.dot(ROTmatrix90,tempDIR) # ROTmatrix90
      tempPIX = np.add(tempPIX,tempDIR)
      if(data2D[tempPIX[0],tempPIX[1]] >= 1):
        listCONTOUR[counter,:] = tempPIX
        listDIRECTION[counter,:] = tempDIR
        counter += 1
      else:
        tempDIR = np.dot(ROTmatrix90,tempDIR) # ROTmatrix90
        tempPIX = np.add(tempPIX,tempDIR)
        if(data2D[tempPIX[0],tempPIX[1]] >= 1):
          listCONTOUR[counter,:] = tempPIX
          listDIRECTION[counter,:] = tempDIR
          counter += 1
        else:
          tempPIX = np.add(tempPIX,tempDIR)
          if(data2D[tempPIX[0],tempPIX[1]] >= 1):
            listCONTOUR[counter,:] = tempPIX
            listDIRECTION[counter,:] = tempDIR
            counter += 1
          else:
            tempDIR = np.dot(ROTmatrix90,tempDIR) # ROTmatrix90
            tempPIX = np.add(tempPIX,tempDIR)
            if(data2D[tempPIX[0],tempPIX[1]] >= 1):
              listCONTOUR[counter,:] = tempPIX
              listDIRECTION[counter,:] = tempDIR
              counter += 1
            else:
              tempPIX = np.add(tempPIX,tempDIR)
              if(data2D[tempPIX[0],tempPIX[1]] >= 1):
                listCONTOUR[counter,:] = tempPIX
                listDIRECTION[counter,:] = tempDIR
                counter += 1
              else:
                tempDIR = np.dot(ROTmatrix90,tempDIR) # ROTmatrix90
                tempPIX = np.add(tempPIX,tempDIR)
                if(data2D[tempPIX[0],tempPIX[1]] >= 1):
                  listCONTOUR[counter,:] = tempPIX
                  listDIRECTION[counter,:] = tempDIR
                  counter += 1
                else:
                  tempPIX = np.add(tempPIX,tempDIR)
                  if(data2D[tempPIX[0],tempPIX[1]] >= 1):
                    listCONTOUR[counter,:] = tempPIX
                    listDIRECTION[counter,:] = tempDIR
                    counter += 1
                  else:
                    tempDIR = np.dot(ROTmatrix90,tempDIR) # ROTmatrix90
                    tempPIX = np.add(tempPIX,tempDIR)
                    if(data2D[tempPIX[0],tempPIX[1]] >= 1):
                      listCONTOUR[counter,:] = tempPIX
                      listDIRECTION[counter,:] = tempDIR
                      counter += 1
                    else:
                      tempDIR = np.dot(ROTmatrix90,tempDIR) # ROTmatrix90
                      tempPIX = np.add(tempPIX,tempDIR)
                      if(data2D[tempPIX[0],tempPIX[1]] >= 1):
                        listCONTOUR[counter,:] = tempPIX
                        listDIRECTION[counter,:] = tempDIR
                        counter += 1

    #==loop ends when *second* pixel is revisited, entered from the same direction as it was entered on its first visit
    if(counter >= 3 and (tempPIX[0] == listCONTOUR[1,0] and tempPIX[1] == listCONTOUR[1,1]) and (tempDIR[0] == listDIRECTION[1,0] and tempDIR[1] == listDIRECTION[1,1])):
      checkLOOPEND = 1
  #==remove last pixel (repeat of 2nd pixel...keep the repeated first for shapefiles)...also remove missing slots
  listCONTOUR = listCONTOUR[0:counter-1,:]

  return tuple([listCONTOUR[:,0],listCONTOUR[:,1]]) #==return a tuple of inds so that it is readily used in other 2D grids

#############################################################################################################
def hysteresis(data,upper_thresh=0.0001,lower_thresh=0.0001):
  FUNCTION_NAME = 'hysteresis'

  logging.info(FUNCTION_NAME + ': Process started')

  ny,nx = np.shape(data)
  radius = 1
  label = 1
  label_grid = np.zeros((ny,nx),dtype=int)

  #fix thresholds
  if(upper_thresh < lower_thresh):
    logging.warning(FUNCTION_NAME + ': ' + str(upper_thresh) + ' is < ' + str(lower_thresh) + '...fixing')
    temp = lower_thresh
    lower_thresh = upper_thresh
    upper_thresh = temp
  
  logging.info(FUNCTION_NAME + ': lower threshold = ' + str(lower_thresh) + '; upper threshold = ' + str(upper_thresh))
 
  for ii in range(2,ny-2):
    for jj in range(2,nx-2):
      if(data[ii,jj] >= upper_thresh and label_grid[ii,jj] == 0):
        label_grid[ii,jj]=label
        #initiate lists
        yList = [ii]
        xList = [jj]
        #--- make sure lists aren't empty and we're not in the first row or column
        if(xList[0] == 0 or yList[0] == 0 or xList[0] == nx-2 or yList[0] == ny-2): continue
        while(len(xList) > 0):
          if(xList[0] == 0 or yList[0] == 0 or xList[0] == nx-2 or yList[0] == ny-2):
            xList = xList[1:]
            yList = yList[1:]
          else:
            for nn in range(-radius+yList[0],radius+yList[0]+1):
              for mm in range(-radius+xList[0],radius+xList[0]+1):
                if(data[nn,mm] >= lower_thresh and label_grid[nn,mm] == 0):
                  label_grid[nn,mm] = label
                  yList.append(nn)
                  xList.append(mm)
            #remove the (first) x and y from the lists
            xList = xList[1:]
            yList = yList[1:]
        #endwhile
        label += 1 #--update label

  logging.info(FUNCTION_NAME + ': Process ended')

  return label_grid
#############################################################################################################
def contour_grid(label_grid,simplify=-1):
  
  FUNCTION_NAME = 'contour_grid'
  logging.info(FUNCTION_NAME + ': Process started. Tracing objects.')

  if(simplify > 0):
    if(simplify > 0.5):
      simplify = 0.5
    elif(simplify < 0.0025):
      simplify = 0.0025
    logging.info(FUNCTION_NAME + ': Using simplification factor of ' + str(simplify))

  contours = {}

  uniq_labels = np.unique(label_grid)
  uniq_labels = uniq_labels[np.where(uniq_labels > 0)]

  for lab in uniq_labels:
    good_ind = np.where(label_grid == lab)
    first_pixel = [good_ind[0][0],good_ind[1][0]]
    inds = trace_object(label_grid,first_pixel)

    if(len(inds[0]) < 4): continue

    if(simplify > 0):
      poly = Polygon([(x,y) for x,y in zip(inds[0],inds[1])])
      h = poly.simplify(simplify,preserve_topology=False)
      if(isinstance(h, shapely.geometry.polygon.Polygon)):
        #this line is taking the output of simplify (python 'array' obj), and creating a tuple of numpy 'ndarrays' of type int
        inds = (np.array(np.asarray(h.exterior.coords.xy[0])).astype(np.int),np.array(np.asarray(h.exterior.coords.xy[1])).astype(np.int))
    
    contours[str(int(lab))] = inds

  logging.info(FUNCTION_NAME + ': Process ended')

  return contours
#############################################################################################################
def write_shapes(thresh_contours,lats,lons,outdir,start_dt,end_dt):
  FUNCTION_NAME = 'write_shapes'
  logging.info(FUNCTION_NAME + ': Process started')

  for t in thresh_contours:
    if(thresh_contours[t]): #if not empty
      w = shapefile.Writer(shapefile.POLYGON)
      w.field('LABEL','N')
      w.field('THRESH','N',decimal=2)
      for label in thresh_contours[t]:
        good_lats = lats[thresh_contours[t][label]]
        good_lons = lons[thresh_contours[t][label]]
        coords = [[[lon,lat] for lon,lat in zip(good_lons,good_lats)]]
        w.record(label,t)
        w.poly(parts=coords, shapeType=shapefile.POLYGON)

      try: #trying to write out shapefiles
        w.save(outdir + ('_').join(['MESH',str(t).replace('.','p'),'in','accum',start_dt,end_dt]))
      except (OSError,IOError) as err:
        logging.error(FUNCTION_NAME + ': ' + str(err))
      else:
        logging.info(FUNCTION_NAME + ': Wrote out ' + outdir + ('_').join(['MESH',str(t).replace('.','p'),'in','accum',start_dt,end_dt])+'.[shp,shx,dbf].')

  logging.info(FUNCTION_NAME + ': Process ended.')
#############################################################################################################
def main():
  parser = argparse.ArgumentParser()
  parser.description = "This program will make shapefiles for daily MRMS MESH tracks in realtime, or create shapefiles for MRMS MESH tracks for custom timeframes in offline mode."
  parser.add_argument("-f","--file",help="[optional] Full path and filename to MRMS grib2 file, for realtime. Default = 'http://mrms.ncep.noaa.gov/data/2D/MESH/MRMS_MESH_Max_1440min.latest.grib2.gz",type=str,nargs=1,default=['http://mrms.ncep.noaa.gov/data/2D/MESH_Max_1440min/MRMS_MESH_Max_1440min.latest.grib2.gz'])
  parser.add_argument("-fp","--filepattern",help="[optional] Full filepattern to process local archived files. E.g., '/home/user/MRMS_MESH_Max_1440min_*_20180310_2000*grib2', or '/home/user/MRMS_MESH_*_20180310-*grib2'. In this way, whether you have a 1440min accumulation or MESH grib2s every 2 min you can create the shapefiles.",type=str,nargs='*')
  parser.add_argument("-t","--thresholds",help="[optional] Thresholds to contour MESH. Default = 0.5, 1.0, 2.0 (inches)",type=str,nargs=3,default=[0.5,1.0,2.0])
  parser.add_argument("-o","--outdir",help="[optional] Output directory. Default = ${PWD}",type=str,nargs=1,default=[os.environ['PWD']])
  parser.add_argument("-hr","--hour",help="[optional] For realtime mode, provide the UTC hour that you want the 24-hr accumulation (e.g., 06). Default = 12",type=str,nargs=1,default=["12"])
  parser.add_argument("-r","--realtime",help="[optional] Run in realtime mode. This means it will run indefinitely. If not set, you should specify a filepattern with -fp.",action="store_true")
  parser.add_argument("-x","--execute_once",help="[optional] Execute once, but for realtime data. This can be used in conjunction with crontab.",action="store_true")
  parser.add_argument("-s","--simplify",help="[optional] Use this simplification factor to reduce the number of vertices. Try 0.0025 first, then increase the factor if you need to. You should not go over a factor of 0.5.",type=np.float, nargs=1)
  parser.add_argument("-l","--logfile",help="[optional] Logfile to write output. Will always write to screen regardless.", type=str, nargs=1)
  parser.add_argument("-rn","--read_netcdf",help="[optional] Read in this netcdf file (an accumulated grid). This is for testing.",type=str,nargs=1)
  parser.epilog = "Example for realtime: [user@machine]$ python dailyHailSwath.py -o RTshapes -l RTshapes/rt.log -r"
  args = parser.parse_args()

  logfile = args.logfile[0] if(args.logfile) else None
  netcdf_file = args.read_netcdf[0] if(args.read_netcdf) else None
  filepattern = args.filepattern if(args.filepattern) else []
  simpfact = args.simplify[0] if(args.simplify) else -1
  thresholds = [np.float(x) for x in args.thresholds]

  dailyHailSwath(file=args.file[0],\
                 filepattern=filepattern,\
                 thresh=thresholds,\
                 outdir=args.outdir[0],\
                 logfile=logfile,\
                 utchour=args.hour[0],\
                 realtime=args.realtime,\
                 execute_once=args.execute_once,\
                 simplify=simpfact,\
                 netcdf_file=netcdf_file)
#############################################################################################################
def dailyHailSwath(file=file,\
                   thresh=[0.5,1.0,2.5],\
                   outdir=os.environ['PWD'],\
                   logfile=None,\
                   utchour='12',\
                   filepattern=[],\
                   realtime=False,\
                   execute_once=False,\
                   simplify=-1,\
                   **kwargs):

  FUNCTION_NAME = 'dailyHailSwath'

  #checking outdir
  if(outdir[-1] != '/'): outdir += '/'
  if(not(os.path.isdir(outdir))):
    try:
      os.makedirs(outdir)
    except (IOError,OSError) as err:
      print("Can't make " + outdir + "!")
      sys.exit()

  #define the logger
  logDict = logging_def(logfile = logfile)
  dictConfig(logDict) 

  logging.info(FUNCTION_NAME + ': Process started')

  #fill utchour
  utchour = utchour.zfill(2)
  
  #sort thresholds. Make sure they are >= 0 and <= 10
  if(len(thresh) != 3):
    logging.error(FUNCTION_NAME + ': You can only do exactly three (3) thresholds.')
    sys.exit()
  thresh = np.sort(thresh)
  try: #making sure they are numbers
    test1 = thresh[0]/1; test2 = thresh[1]/1; test3 = thresh[2]/1
  except TypeError as err:
    logging.error(FUNCTION_NAME + ': ' + str(err))
    sys.exit()
  #sorry, you can only do 0 <= MESH < 10
  if(thresh[0] < 0 or thresh[0] >= 10 or thresh[1] < 0 or thresh[1] >= 10 or thresh[2] < 0 or thresh[2] >= 10):
    logging.error(FUNCTION_NAME + ': One of your thresholds is out of range (0 <= threshold < 10)')
    sys.exit()
  

  netcdf_file = kwargs['netcdf_file'] if('netcdf_file' in kwargs) else None

  #archive mode
  if(not(realtime) and not(execute_once)):
    logging.info(FUNCTION_NAME + ': Starting archived case')
    ####this block is testing####
    if(netcdf_file):
      #read the netcdf
      atts={'Latitude':0,'Longitude':0,'LatGridSpacing':0.01,'LonGridSpacing':0.01}
      mesh_accum = read_data(netcdf_file,'MESH_Max_1440min',global_atts=atts) * 0.0393701

      ny,nx = np.shape(mesh_accum)
      NW_lat = atts['Latitude']; NW_lon = atts['Longitude']; dy = atts['LatGridSpacing']; dx = atts['LonGridSpacing']
      lats = np.zeros((ny,nx)); lons = np.zeros((ny,nx))
      tmpdt = datetime.strptime(os.path.basename(netcdf_file)[0:15],'%Y%m%d-%H%M%S')
      end_dt = tmpdt.strftime('%Y%m%d-%H%M%S')
      start_dt = (tmpdt - timedelta(days=1)).strftime('%Y%m%d-%H%M%S')
      for ii in range(ny):
        lats[ii,:] = NW_lat - ii * dy
      for jj in range(nx):
        lons[:,jj] = NW_lon + jj * dx
    else:
      mesh_files = filepattern 
      if(len(mesh_files) == 1 and '*' in filepattern[0]):
        logging.error(FUNCTION_NAME + ': No files exist for ' + filepattern[0])
        sys.exit()
      mesh_files = np.sort(mesh_files)
      #e.g.,MRMS_MESH_00.50_20170511-120437.grib2
      #e.g.,MRMS_MESH_Max_1440min_latest.grib2
      #e.g.,MRMS_MESH_Max_1440min_00.50_20180404-123832.grib2

      nfiles = len(mesh_files)
      for mm in range(nfiles):
        logging.info(FUNCTION_NAME + ': Processing ' + mesh_files[mm])
        #gunzip if necessary
        if(mesh_files[mm][-3:] == '.gz'):
          s = call(['gunzip','-f',mesh_files[mm]])
          grib2 = mesh_files[mm][:-3]
        else:
          grib2 = mesh_files[mm]
        #read grib2
        grb = pg.open(grib2)
        #get data
        msg = grb[1]  #should only be 1 'message' in grib2 file
  
        #get the start and end time info
        if(nfiles == 1):
          if(int(msg['parameterName']) == 34): #a 1440min file
            end_dt = (msg.validDate).strftime('%Y%m%d-%H%M%S')
            start_dt = (datetime.strptime(end_dt,'%Y%m%d-%H%M%S') - timedelta(days=1)).strftime('%Y%m%d-%H%M%S')
          elif(int(msg['parameterName']) == 28): #a 2-min file...I don't know why you would want to make shapes for a 2-min increment, but oh well.
            end_dt = (msg.validDate).strftime('%Y%m%d-%H%M%S')
            start_dt = (msg.validDate).strftime('%Y%m%d-%H%M%S')
          else:
            logging.error(FUNCTION_NAME + ': Unknown parameterName:'+str(msg['parameterName']))
        else:
          if(mm == 0):
            if(int(msg['parameterName']) == 34): #a 1440min file
              start_dt_tmp = (msg.validDate).strftime('%Y%m%d-%H%M%S')
              start_dt = (datetime.strptime(start_dt_tmp,'%Y%m%d-%H%M%S') - timedelta(days=1)).strftime('%Y%m%d-%H%M%S')
            elif(int(msg['parameterName']) == 28): # a 2-min file
              start_dt = (msg.validDate).strftime('%Y%m%d-%H%M%S')
            else:
              logging.error(FUNCTION_NAME + ': Unknown parameterName:'+str(msg['parameterName']))
          elif(mm == (nfiles - 1)):
            end_dt = (msg.validDate).strftime('%Y%m%d-%H%M%S')

        mesh_data = msg.values * 0.0393701 #convert from mm to in
        ny,nx = np.shape(mesh_data)
        if('lats' not in locals()): lats = (grb[1]['latitudes']).reshape(ny,nx)
        if('lons' not in locals()): lons = (grb[1]['longitudes']).reshape(ny,nx) - 360.
        #close grib2
        grb.close()
        #add to accumulation grid...after creating mesh_accum grid, if necessary
        if('mesh_accum' not in locals()): mesh_accum = np.zeros((ny,nx),dtype=np.float) - 999.
        mesh_accum = np.maximum(mesh_accum,mesh_data)

    mesh_accum = max_filter(mesh_accum)

    #make label grids for the different thresholds, using region growing
    thresh1_label = hysteresis(mesh_accum,upper_thresh=thresh[0],lower_thresh=thresh[0])
    thresh2_label = hysteresis(mesh_accum,upper_thresh=thresh[1],lower_thresh=thresh[1])
    thresh3_label = hysteresis(mesh_accum,upper_thresh=thresh[2],lower_thresh=thresh[2])

    #make contours based on label grids
    thresh_contours = {}

    logging.info(FUNCTION_NAME + ': Contouring threshold1(' + str(thresh[0]) + ') labels')
    thresh_contours[thresh[0]] = contour_grid(thresh1_label,simplify=simplify) 
    
    logging.info(FUNCTION_NAME + ': Contouring threshold2(' + str(thresh[1]) + ') labels')
    thresh_contours[thresh[1]] = contour_grid(thresh2_label,simplify=simplify)
 
    logging.info(FUNCTION_NAME + ': Contouring threshold3(' + str(thresh[2]) + ') labels')
    thresh_contours[thresh[2]] = contour_grid(thresh3_label,simplify=simplify)

    #write shapefile
    write_shapes(thresh_contours,lats,lons,outdir,start_dt,end_dt)


  else:
  #now starting infinite loop for realtime processing
    logging.info(FUNCTION_NAME + ': Starting realtime processing for 24 accumulations from ' + utchour + 'Z to ' + utchour + 'Z.')
    lastprocess = 0
    while(True):
      #get time
      now = datetime.now()
      datetimestr = now.strftime('%Y%m%d-%H%M%S')
      dt_hr = now.strftime('%H')
      nowsecs = int(now.strftime('%s'))
    
      #get realtime file. First try to remove the old one, if present 
      if(os.path.isfile(file)):
        try:
          logging.info(FUNCTION_NAME + ': Removing old ' + file)
          os.remove(file)
        except (OSError,IOError) as err:
          logging.error(FUNCTION_NAME + ': ' + str(err))
 
      logging.info(FUNCTION_NAME + ': Fetching ' + file)
      status = call(['wget',file])
      if(status != 0):
        logging.error(FUNCTION_NAME + ': Something went wrong fetching ' + file)
        logging.info(FUNCTION_NAME + ': Sleeping 30 min')
        time.sleep(30*60) 
        continue

      logging.info(FUNCTION_NAME + ': Processing ' + os.path.basename(file))
      #gunzip if necessary
      thegrib = os.path.basename(file)
      if(thegrib[-3:] == '.gz'):
        s = call(['gunzip','-f',thegrib])
        grib2 = thegrib[:-3]
      else:
        grib2 = thegrib
      #read grib2
      grb = pg.open(grib2)
      #get data
      msg = grb[1]  #should only be 1 'message' in grib2 file
      #check time. If we have the correct hour and we last processed over 2 hours ago, then go ahead and process
      if((msg.validDate.strftime('%H') == utchour and (nowsecs - lastprocess) > 3600*2) or execute_once):

        end_dt = msg.validDate
        end_dtSTR = msg.validDate.strftime('%Y%m%d-%H%M%S')
        start_dt = (end_dt - timedelta(days=1)).strftime('%Y%m%d-%H%M%S')

        mesh_accum = msg.values * 0.0393701
        ny,nx = np.shape(mesh_accum)
        lats = (grb[1]['latitudes']).reshape(ny,nx)
        lons = (grb[1]['longitudes']).reshape(ny,nx) - 360.
        #close grib2
        grb.close()

        mesh_accum = max_filter(mesh_accum)

        #make label grids for the different thresholds, using region growing
        thresh1_label = hysteresis(mesh_accum,upper_thresh=thresh[0],lower_thresh=thresh[0])
        thresh2_label = hysteresis(mesh_accum,upper_thresh=thresh[1],lower_thresh=thresh[1])
        thresh3_label = hysteresis(mesh_accum,upper_thresh=thresh[2],lower_thresh=thresh[2])

        #make contours based on label grids
        thresh_contours = {}

        logging.info(FUNCTION_NAME + ': Contouring threshold1(' + str(thresh[0]) + ') labels')
        thresh_contours[thresh[0]] = contour_grid(thresh1_label,simplify=simplify)

        logging.info(FUNCTION_NAME + ': Contouring threshold2(' + str(thresh[1]) + ') labels')
        thresh_contours[thresh[1]] = contour_grid(thresh2_label,simplify=simplify)

        logging.info(FUNCTION_NAME + ': Contouring threshold3(' + str(thresh[2]) + ') labels')
        thresh_contours[thresh[2]] = contour_grid(thresh3_label,simplify=simplify)

        #write shapefile
        write_shapes(thresh_contours,lats,lons,outdir,start_dt,end_dtSTR)

        #reset process time
        lastprocess = nowsecs

      else:
        logging.info(FUNCTION_NAME + ': Hour = ' + msg.validDate.strftime('%H') + 'Z; waiting for ' + utchour + 'Z')
        grb.close()  #close grib


      logging.info(FUNCTION_NAME + ': Removing ' + grib2)
      try:
        os.remove(grib2)
      except (OSError,IOError) as err:
        logging.warning(FUNCTION_NAME + ': ' + str(err))

      if(execute_once):
        logging.info(FUNCTION_NAME + ': Executed once and now exiting.')
        sys.exit()

      logging.info(FUNCTION_NAME + ': Sleeping 30 min')
      time.sleep(30*60)
    #endwhile

  logging.info(FUNCTION_NAME + ': Process ended.')

#############################################################################################################
if __name__ == "__main__":
  main()
