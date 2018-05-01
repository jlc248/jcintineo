import numpy as np
import shapefile as sp
import json
import sys
import os
from glob import glob

#20150424-25  20160410-11  20160412-13  20160509-10  20160511-12  20170228-0301  20170508-09  20170611-12
#20150619-20  20160411-12  20160508-09  20160510-11  20170225-26  20170326-27    20170518-19  20170629-30

dir = 'TESTTEST'

dates = ['20150424','20150619','20160410','20160411','20160412','20160508','20160509','20160510','20160511',\
         '20170225','20170228','20170326','20170508','20170518','20170611','20170629']
dates=['20180403']
for date in dates:
  print(date)
  files = glob(dir + '/*accum_' + date + '*.shp')

  outDict = {}
  outDict['type'] = 'FeatureCollection'
  outDict['features'] = []

  for file in files:

    sf = sp.Reader(file)

    shapes = sf.shapes()
    fields = sf.fields
    recs = sf.records()

    sRecs = sf.shapeRecords()

    for rec in sRecs:
  
      label = rec.record[0]
      threshold = rec.record[1]

      if(threshold < 1):
        color = '255 165 0'
        BWIDTH='3'
      elif(threshold < 2):
        color = '255 0 0'
        BWIDTH='3'
      else:
        color = '255 0 255'
        BWIDTH='4'


      info='Label: ' + str(label) + ' Threshold: ' + str(threshold) + ' in.'

      coords = [[[np.round(p[0],2),np.round(p[1],2)] for p in rec.shape.points]]

      feature = {'type':'Feature','properties':{'WIDTH':BWIDTH,'INFO':info,'COLOR':color,'BCOLOR':color,'BOPACITY':'100','OPACITY':'0'},\
                 'geometry':{'type':'Polygon','coordinates':coords}}

      outDict['features'].append(feature)


  with open(dir + '/MESHshapes_'+date+'_120000.force.json','w') as out:
    json.dump(outDict,out,indent=2)
