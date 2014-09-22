# Extracts lat-lons for all london hire bicycle stations.  Importing https in
# python requires ssl and is thus fraught with OS system differences and
# unreliable routines in urllib. To circumvent problems, this routine can be run
# by simply downloading the source of the webpage at
# https://web.barclayscyclehire.tfl.gov.uk/maps/
# which should be titled "maps.htm". The python script then extracts the lats
# and lons of all stations from this html source.

from bs4 import BeautifulSoup
import re
f = open ('../data/maps.htm', 'r')
page = f.read ()
f.close ()
soup = BeautifulSoup (page)
list = soup.findAll ('script')
list = [i for i in list if i.find (text=True)]
list = [i for i in list if
        i.find(text=True).encode('utf-8').find('stationMarker') >= 0] [0]
# list is a single bs4 object, which is then converted to a string
list = list.find (text=True).encode ('utf-8') 
list = list.split ('station=')
fname = '../data/station_latlons.txt'
f = open (fname, 'w')
outs = ('id', 'lat', 'long', 'name')
f.write ('id, lat, long, name\n')
count = 0
minstn = 9999
maxstn = 0
for li in list:
    if li.find ("name") > 0:
        # This regex extracts all text NOT contained within double quotes, which
        # thus contains all field names of the javascript source
        linames = re.findall ('(?:^|")([^"]*)(?:$|")', li)
        linames = [filter (None, re.findall ('([^{,:]*)', i))[0] for i in linames]
        # And this regex extract all text within double quotes, which are the
        # corresonding field values.
        lidat = re.findall ('"([^"]*)"', li)
        # f is written as comma-delimited, so remove commas from station names
        lidat = [i.replace (',', ' ') for i in lidat]
        indx = [linames.index (i) for i in outs]
        if len (indx) == len (outs):
            for i in indx:
                f.write (lidat [i])
                if i != indx [-1]:
                    f.write (',')

            count = count + 1
            f.write ('\n')
            # Update max & min station numbers
            stnum = int (lidat [linames.index ('id')])
            if stnum < minstn:
                minstn = stnum
            elif stnum > maxstn:
                maxstn = stnum

f.close ()
print count, 'stations in [', minstn, ',', maxstn, '] written to ', fname
