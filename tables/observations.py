"""
Creates a table observations from the fits files arranged into the default file structure, as used by GALAH

Authors: Janez Kos (jkos@physics.usyd.edu.au), Sanjib Sharma, Ghayandi de Silva, Sarah Martell, GALAH collaboration

Copyright (C) 2015  Janez Kos

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



HOW TO RUN THIS CODE:

1. Make sure you have the following modules installed:
- os
- pyfits
- ephem
- sys
- numpy
- psycopg2
- argparse

2. Install the postgres sql server

3. Create a new postgres user and a database

4. Set the username and database name in the main()

5. Run with the only argument being the path to the folder in which the ususal folder structure of the GALAH data is.

The code will produce a sql database, its dump and a csv table.
"""

from os import walk, system, path
import psycopg2 as mdb
import pyfits
import numpy as np
import ephem
import sys
import argparse

def lbr2xyz(l,b,r=1.0):
#    dtor=np.pi/180.0
	if l==None: return 999
	if b==None: return 999
	l=np.radians(l)
	b=np.radians(b)
	return [r*np.cos(b)*np.cos(l),r*np.cos(b)*np.sin(l),r*np.sin(b)]

def angsep(l1,b1,l2,b2):
	if l1==None or b1==None or l2==None or b2==None: return 999
	x1,y1,z1=lbr2xyz(l1,b1)
	x2,y2,z2=lbr2xyz(l2,b2)
	return np.degrees(2*np.arcsin(np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/2.0))

def label2obsid(label):
    status=0
    for temp in label.split(','):
        if temp.split('=')[0].strip() == 'omversion':
            status=1    
    if status == 1:
        result=[]
        for i,temp in enumerate(label.split(',')):
            x=temp.split('=')
            if len(x)>1:
                y=x[1].strip()
            else:
                y=x[0].strip()
            if (i==1) or (i==4):
                y=int(y)
            result.append(y)
        if len(result)==4:
            result.append(result[1])
            result.append(0)
        if len(result)<6:
            #print label
            result=['',-1,'','',-1]
    else:
        result=['',-1,'','',-1]
    return result[1],result[4]

def read_comments(comments,ymd,run):
	"""
	will read comments file  and find obstatus
	"""
	run=int(run)
	data=[]
	try:
		comments=comments+"/"+str(ymd)+"/comments/"+"comments_%s.txt" % (str(ymd))
		data=np.loadtxt(comments, dtype=np.str, delimiter=',')
		runn=data[:,0]
	except:
		return None
	obstatus=data[:,1]
	ind= np.where(runn==str(run))
	if len(ind[0])>=1: ind=ind[0]
	else: return None
	return obstatus[ind][0]

def read_comments2(comments,ymd,run):
	"""
	will read comments file and find dome flats
	"""
	run=int(run)
	data=[]
	try:
		comments=comments+"/"+str(ymd)+"/comments/"+"comments_%s.txt" % (str(ymd))
		data=np.loadtxt(comments, dtype=np.str, delimiter=',')
		runn=data[:,0]
	except:
		return 0

	cmt=data[:,4]
	ind= np.where(runn==str(run))
	if len(ind[0])>=1: ind=ind[0]
	else: return 0
	if cmt[ind][0].find('dome')!=-1 and cmt[ind][0].find('flat')!=-1:
		return 1
	else:
		return 0

def check_channels(fits):
	channels=[1,2,3,4]
	folders=['ccd_1', 'ccd_2', 'ccd_3', 'ccd_4']

	this_channel=int(fits[-10])
	this_folder=fits.split("/")[-2]
	channels = [x for x in channels if x != this_channel]
	folders=[x for x in folders if x != this_folder]
	
	for i,j in zip(channels, folders):
		fits_tmp=fits.split('/')
		fits_tmp[-2]=j
		b= bytearray(fits_tmp[-1])
		b[-10] = str(i)
		fits_tmp[-1]=str(b)
		to_check= '/'.join(fits_tmp)

		try:
			h0=pyfits.getheader(to_check,0)
		except:
			return 1

		return 0

def moon_distance(utdate, utstart, ra, de):
	obs = ephem.Observer()
	obs.lon = '149.0658'
	obs.lat = '-31.2769'
	obs.elevation = 1164
	obs.date = "/".join(utdate.split(":"))+" "+utstart
	moon= ephem.Moon()
	moon.compute(obs)

	star = ephem.FixedBody()
	star._ra = np.radians(ra)
	star._dec = np.radians(de)
	star.compute(obs)

	sep=ephem.separation(moon,star)
	if sep< ephem.degrees('30:0:0.0'): 
		return 1
	else: return 0

def saturated(image):
	f = pyfits.open(image)
	scidata = f[0].data
	f.close()

	if len(np.where(scidata[0:-1, 1500:2500]>65000)[0])>5000: return 1
	else: return 0

def check_signal(image,ndfclass):
	if ndfclass!='MFOBJECT': return 0
	else:
		f = pyfits.open(image)
		scidata = f[0].data
		f.close()

		if np.average(scidata[0:-1, 1500:2500])<400: return 1
		else: return 0

def sun_elev(utdate, utstart):
	obs = ephem.Observer()
	obs.lon = '149.0658'
	obs.lat = '-31.2769'
	obs.elevation = 1164
	obs.date = "/".join(utdate.split(":"))+" "+utstart
	sun= ephem.Sun()
	sun.compute(obs)

	if sun.alt<ephem.degrees('-8:0:0.0'): return 0
	else: return 1

def check_exptime(ndfclass, exposed, dirname):
	exposed=float(exposed)
	if ndfclass=='BIAS' and exposed>0.1: return 1
	elif ndfclass=='MFARC' and exposed>400.0: return 1
	elif ndfclass=='MFOBJECT' and exposed>3000.0: return 1
	elif ndfclass=='MFFLX' and exposed>20.0 and dirname<20140400: return 1
	elif ndfclass=='MFFLX' and exposed>200.0 and dirname>20140400: return 1
	else: return 0


def read_files(i, fullname):
	#things you read from the image name
	ymd=int(fullname[-4])
	dirname=ymd
	ccd=fullname[-2][-1]
	fitsfile=fullname[-1]
	run=fullname[-1][6:10]

	#key will be used as a unique primary key
	key=(dirname*10000+int(run))*10+int(ccd)
	runkey=dirname*10000+int(run)

	#check if filesize is ok
	fsize=path.getsize(i)
	if fsize<33000000:
		return key, runkey, dirname, fitsfile, run, -1, ccd, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0, 64, 'NULL', 'NULL', 'NULL'

	#things you read from the header
	h0=pyfits.getheader(i,0)

	utmjd=h0["UTMJD"]
	utdate=h0["UTDATE"]
	utstart=h0["UTSTART"]
	try:#because some commissioning data has no NDFCLASS
		ndfclass=h0["NDFCLASS"]
	except:
		ndfclass='NONE'

	meanra=h0["MEANRA"]
	meandec=h0["MEANDEC"]
	zenith=h0["ZDSTART"]
	airmass=1.0/np.cos(float(zenith)/180.0*np.pi)
	exposed=h0["EXPOSED"]
	obstype=h0["OBSTYPE"]
	instrument=h0["INSTRUME"]
	try:#plate is sometimes 'invalid' and sometimes it cant be found at all. Set to 'invalid' in the later case
		plate=h0["SOURCE"]
	except:
		plate='invalid'
	if plate=='invalid':
		plate=-1
	else:
		plate=plate.split()[1]
	cfg_file=h0["CFG_FILE"]
	obj_name=h0["OBJECT"]
	try:
		maskstate=h0["SLITMASK"]
	except:
		maskstate='IN'

	h1=pyfits.getheader(i,'STRUCT.MORE.FIBRES')
	cenra=h1["CENRA"]
	cendec=h1["CENDEC"]
	cenra=np.degrees(cenra)
	cendec=np.degrees(cendec)

	#things you deduce from header
	if h0.__contains__('STD_NAME'):
		stdname=h0['STD_NAME']
		fieldid=-2
		obsid=-2
	elif angsep(np.degrees(h1['CENRA']),np.degrees(h1['CENDEC']),h0['MEANRA'],h0['MEANDEC'])>0.05:
		stdname='NULL'
		fieldid=-2
		obsid=-2
	else:
		stdname='NULL'
		fieldid,obsid=label2obsid(h1['LABEL'])
		if h1.__contains__('FILENAME'):
			if h0['CFG_FILE'].startswith('gf0_'):
				x=h0['CFG_FILE'].split('_')
				obsid=int(x[1])
			elif h0['CFG_FILE'].startswith('gf1_'):
				x=h0['CFG_FILE'].split('_')
			elif h0['CFG_FILE'].startswith('kepler_'):
				x=h0['CFG_FILE'].split('_')
				obsid=fieldid+10000
			elif h0['CFG_FILE'].startswith('k2c'):
				x=h0['CFG_FILE'].split('_')
			else:
				fieldid=-1
				obsid=-1
			if (fieldid>=0)and(obsid<0):
				raise RuntimeError('fieldid>=0 and obsid<0')
		else:
			raise RuntimeError('some prob')
			obsid=-1
			fieldid=-1

	#things you read from the fits table:
	hdulist = pyfits.open(i)
	fdata=hdulist["STRUCT.MORE.FIBRES"].data
	tp= fdata["TYPE"].tolist()
	n_of_targets=tp.count('P')
	hdulist.close()

	#things you read from comments file:
	obstatus=read_comments('/'.join(i.split("/")[0:-4]),dirname,run)
	if obstatus==None: obstatus=1 # obstatus = 1 unles comments file says otherwise

	#ndfclass_updated is just the copy of the ndfclass from the header, unless the following applies
	ndfclass_updated=ndfclass

	#==============FIX NDFCLASS====================
	#for twilight flats:
	if ndfclass=='MFOBJECT' and sun_elev(utdate, utstart) and (int(exposed) in [120, 240, 480]) and angsep(float(meanra),float(meandec),float(cenra),float(cendec))>0.1:
		ndfclass_updated='SFLAT'

	#for dome flats:
	if obstype=='DOME':
		ndfclass_updated='DOME'
	read_comments2('/'.join(i.split("/")[0:-4]),dirname,run)
	if read_comments2('/'.join(i.split("/")[0:-4]),dirname,run):
		ndfclass_updated='DOME'

	#==============SET OBSTATUS====================

	#set obstatus to 0 in the following cases:
	if ndfclass not in ['BIAS', 'DARK', 'MFFFF', 'MFARC', 'MFOBJECT', 'MFFLX', 'MFSKY', 'SFLAT']: obstatus=0

	#set obstatus if BIAS has a long exposure
	if ndfclass=='BIAS' and exposed>=0.1: obstatus=0

	#set obstatus if BIAS has too much signal:
	if ndfclass=='BIAS':
		f = pyfits.open(i)
		scidata = f[0].data
		mean=np.mean(scidata[1500:2500, 1500:2500])

		if mean>400.0: obstatus=0#if mean signal in bias is more than 400, it is invalid
		f.close()

	#set obstatus=0 if plate cannot be found in the header
	if plate==-1:
		obstatus=0

	#set obstatus=0 if AAOmega data leaked into GALAH directories
	if instrument=='AAOMEGA-2dF': obstatus=0


	#==============SET QFLAG====================
	qflag=0# by default everything is ok

	#check if data exists in all four chanels. If not, qflag+=1
	if check_channels(i): qflag+=1

	#check distance from the moon. If distance is <30 deg, qflag+=2
	if moon_distance(utdate, utstart,meanra, meandec): qflag+=2

	#check if there are a lot of saturated pixels. If number of saturated pixels is > 20 000, qflag+=4
	if saturated(i): qflag+=4

	#check if there is low signal in images of regular fields. This could mean cloudy weather.
	if check_signal(i,ndfclass): qflag+=8

	#check if the slit mask was in:
	if maskstate=='IN': qflag+=16

	#check if exposure time seem peculiar:
	if check_exptime(ndfclass, exposed, dirname): qflag+=32

	#if filesize is not ok
	#qflag is set to 64 when checking the filesize and breaking out of this function

	#==============SET OCLASS====================
	bensby=[(140712, 56), (140712, 57), (140712, 58), (140712, 59), (140712, 60), (140114, 50), (140114, 51), (140114, 52), (140114, 47), (140114, 48), (140114, 49), (140114, 53), (140114, 54), (140114, 55), (140113, 47), (140114, 48), (140114, 49),  (140113, 47),  (140113, 48),  (140113, 49),  (140113, 53),  (140113, 54),  (131119, 17),  (131119, 18),  (131119, 19),  (140112, 33),  (140112, 34),  (140111, 44),  (131219, 24),  (131219, 25),  (131219, 26), (140114, 58), (140114, 59), (140114, 60), (140114, 61), (140114, 62), (140114, 63), (140807, 17), (140807, 18), (140807, 21), (140807, 22), (140807, 25), (140807, 26), (140807, 27), (140807, 28), (140807, 31), (140807, 32), (140713, 39), (140712, 53), (140712, 54), (140712, 55)]
	telluric=[(131114, 57), (131114, 58), (131114, 59), (131114, 66), (131114, 67), (131114, 68), (131114, 72), (131114, 73), (131114, 74), (141202, 44), (141202, 45), (141202, 46), (141202, 49), (141202, 50), (141202, 51), (141202, 52), (141202, 53), (141202, 54), (150206, 31), (141202, 55), (150205, 31), (150205, 39), (150205, 47), (150205, 55), (140118, 40), (140118, 41), (140111, 42), (140111, 43), (150204, 34), (150204, 35), (150204, 36), (131114, 60), (131114, 61), (131114, 62), (131114, 63), (131114, 64), (131114, 65), (131114, 69), (131114, 70), (131114, 71)]
	if cfg_file.find('benchmark')!=-1 or cfg_file.find('Benchmark')!=-1 or cfg_file.find('BENCHMARK')!=-1 or (ndfclass=='MFFLX' and angsep(float(meanra),float(meandec),float(cenra),float(cendec))>0.1 and int(exposed)<600):
		oclass='benchmark'
	elif cfg_file.find('bensby')!=-1 or cfg_file.find('Bensby')!=-1 or cfg_file.find('BENSBY')!=-1 or ((int(dirname), int(run)) in bensby):
		oclass='bensby'
	elif cfg_file.find('telluric')!=-1 or cfg_file.find('Telluric')!=-1 or cfg_file.find('TELLURIC')!=-1 or ((int(dirname), int(run)) in telluric):
		oclass='telluric'
	elif ndfclass=='MFOBJECT':
		oclass='survey'
	else:
		oclass='NULL'


	#assign cob_id=-1 for now. Will set tis later
	cob_id=-1
	#print key, runkey, dirname, fitsfile, run, cob_id, ccd, plate, cfg_file, meanra,meandec, ndfclass, fieldid, obsid, exposed, obstatus
	return key, runkey, dirname, fitsfile, run, cob_id, ccd, plate, cfg_file, meanra, meandec, cenra,cendec, utmjd, ndfclass, ndfclass_updated, fieldid, obsid, exposed, stdname, obstatus, qflag, oclass, airmass, n_of_targets






def main(dbname, user, motherfolder,add,night):
	con=mdb.connect("dbname=%s user=%s" % (dbname, user))
	cur = con.cursor()


	
	#create a table
	if add==False and night==None:
		cur.execute("DROP TABLE IF EXISTS observations")# be careful! will delete the previous table
		cur.execute("CREATE TABLE observations(runccd_id bigint primary key)")
		cur.execute("ALTER TABLE observations ADD run_id bigint")
		cur.execute("ALTER TABLE observations ADD dirname int")
		cur.execute("ALTER TABLE observations ADD fitsfile varchar(16)")
		cur.execute("ALTER TABLE observations ADD run smallint")
		cur.execute("ALTER TABLE observations ADD cob_id bigint")
		cur.execute("ALTER TABLE observations ADD ccd smallint")
		cur.execute("ALTER TABLE observations ADD plate smallint")
		cur.execute("ALTER TABLE observations ADD cfg_file varchar(256)")
		cur.execute("ALTER TABLE observations ADD meanra float8")
		cur.execute("ALTER TABLE observations ADD meandec float8")
		cur.execute("ALTER TABLE observations ADD cenra float8")
		cur.execute("ALTER TABLE observations ADD cendec float8")
		cur.execute("ALTER TABLE observations ADD utmjd float8")
		cur.execute("ALTER TABLE observations ADD ndfclass varchar(16)")
		cur.execute("ALTER TABLE observations ADD ndfclass_updated varchar(16)")
		cur.execute("ALTER TABLE observations ADD fieldid varchar(16)")
		cur.execute("ALTER TABLE observations ADD obsid varchar(16)")
		cur.execute("ALTER TABLE observations ADD exposed float8")
		cur.execute("ALTER TABLE observations ADD std_name varchar(256)")
		cur.execute("ALTER TABLE observations ADD obstatus smallint") #obstatus 0 or 1
		cur.execute("ALTER TABLE observations ADD qflag bigint")#quality flag. 0 if everything is ok. use bitmasking for everything else
		cur.execute("ALTER TABLE observations ADD oclass varchar(16)")
		#cur.execute("ALTER TABLE observations ADD airmass float8")
		#cur.execute("ALTER TABLE observations ADD targets int")
	elif add==False and night!=None:
		cur.execute("delete from observations where dirname=%s" % (night))
	else:
		pass

	print add



	#a folder in which all the observations are. Let us call it the motherfolder
	#it is a command line argument
	if motherfolder[-1]!="/": motherfolder=motherfolder+"/"


	#create a list of all fits files found in the motherfolder
	fits=[]
	imagepaths=[]

	for root, dirs, files in walk(motherfolder, followlinks=False):
		print " + making filelist for "+root
		for i in files:
			if night==None:
				if i[-5:]=='.fits': 
					if add==True:#check if there is any new data to add
						cur.execute("select count(*) from observations where dirname=%s and ccd=%s and fitsfile='%s'" % (root.split("/")[-3], root.split("/")[-1][-1], i))
						if cur.fetchall()[0][0]==0:
							imagepaths.append(root+"/"+i)
							fits.append(i)
						else: 
							pass
					else:#default
						imagepaths.append(root+"/"+i)
						fits.append(i)
			elif night!=None and root.split("/")[-3]==night:#if only one night must be processed
				if i[-5:]=='.fits': 
					if add==True:
						cur.execute("select count(*) from observations where dirname=%s and ccd=%s and fitsfile='%s'" % (root.split("/")[-3], root.split("/")[-1][-1], i))
						if cur.fetchall()[0][0]==0:
							imagepaths.append(root+"/"+i)
							fits.append(i)
						else: 
							pass
					else:
						imagepaths.append(root+"/"+i)
						fits.append(i)
			else:
				pass

	#do everything for each file seperatly:
	n=0
	N=len(imagepaths)
	for i in imagepaths:
		n+=1
		print " + %s/%s creating the table for file" % (n,N), i
		fullname=i.split("/")
		if fullname[-2][:4]=='ccd_' and fullname[-1][-5:]==".fits" and len(fullname[-1])==15 and fullname[-3]=='data' and len(fullname[-4])==6: #check if the filename looks fine

			#read parameters of each fits file and set everything
			key, runkey, dirname, fitsfile, run, cob_id, ccd, plate, cfg_file, meanra, meandec, cenra,cendec, utmjd, ndfclass, ndfclass_updated, fieldid, obsid, exposed, stdname, obstatus, qflag, oclass, airmass, n_of_targets=read_files(i, fullname)
		

			#write everything into the database
			#print key, dirname, fitsfile
			cur.execute("insert into observations values (%s, %s, %s, '%s', %s, %s, %s, %s, '%s', %s, %s, %s, %s, %s, '%s', '%s', '%s', '%s', %s, '%s', %s, %s, '%s')" % (key, runkey, dirname, fitsfile, run, cob_id, ccd, plate, cfg_file, meanra, meandec, cenra,cendec, utmjd, ndfclass, ndfclass_updated, fieldid, obsid, exposed, stdname, obstatus, qflag, oclass))
		else:
			pass

	con.commit()

	cur.execute("select distinct(dirname) from observations")
	ymds=cur.fetchall()
	
	#now calculate the cob_ids
	for ymd in ymds:
		cur.execute("select distinct(run), ndfclass_updated, meanra, meandec, plate from observations where dirname=%s and ndfclass_updated!='MFSKY' order by run" % ymd)
		ids=cur.fetchall()

		ra_p=-1
		dec_p=-1
		plate_p=-1
		ndfclass_p='none'

		for i in ids:
			if i[1]=='BIAS' or i[1]=='DOME': 
				#set cob_id to -1 for biases
				cur.execute("update observations set cob_id=-1 where dirname=%s and run=%s" % (ymd[0], i[0]))
			else:
				cob_id_tmp=int(ymd[0])*10000+int(i[0])
				if angsep(i[2],i[3],ra_p,dec_p)>0.05 or i[4]!=plate_p: change=1# if the telescope moves more than 0.05 deg, this is a new field
				else: change=0
				if i[1]=='SFLAT' and i[4]==plate_p and i[1]==ndfclass_p: change=0#if we are making twilight flats, don't change the cob_id

				ra_p=i[2]
				dec_p=i[3]
				plate_p=i[4]
				ndfclass_p=i[1]
				if change:
					cob_id=cob_id_tmp
				#print cob_id
				cur.execute("update observations set cob_id=%s where dirname=%s and run=%s" % (cob_id, ymd[0], i[0]))


	con.commit()

if __name__=="__main__":
	#computer/user specific:
	dbname='hermes_master'#name of the postgresql database. Must be created beforehand
	user='jkos'#username for the database.

	#parse the arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("motherfolder", help="path to the folder where galah data is")
	parser.add_argument("--add", help="add to database rather than rewrite the whole database", action="store_true")
	parser.add_argument("--one_night", help="add to database only one night of data. Use yymmdd to specify the night")
	args = parser.parse_args()

	main(dbname, user, args.motherfolder, args.add, args.one_night)
	system("pg_dump -O %s -t observations > observations.sql" % (dbname))
	system("psql %s -c \"COPY (SELECT * from observations order by dirname,run,ccd) TO stdout DELIMITER ',' CSV HEADER\" > observations.csv" % (dbname))
