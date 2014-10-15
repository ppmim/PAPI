################################################################################
#
# PANICtool
#
# dataset.py
#
# Author : jmiguel@iaa.es
#
# History:
# --------
# 15/05/2008  : created  
# 14/04/2009  : added ra,dec fields 
# 25/05/2009  : added object field to DB
# 31/03/2010  : added source as a file_list containing the list files
################################################################################

# Import requered modules
#from pysqlite2 import dbapi2 as sqlite
import sqlite3 as sqlite
import datahandler
import os
import sys
import math
import fileinput
from optparse import OptionParser


###### Enable logging
from misc.paLog import log

__docformat__ = "restructuredtext"  

############################################################
class DataSet(object):  

    """
    Class used to define a data set of frames load from a local directory or a Data Base  
    """

    ############################################################
    TABLE_COLUMNS="(id, run_id, ob_id, ob_pat, expn, nexp, filename, date, \
                    ut_time, mjd, type, filter, texp, ra, dec, object, detector_id)"
                
    # Maximum seconds (10min=600secs aprox) of temporal distant allowed between 
    # two consecutive frames (1/86400.0)*10*60
    MAX_MJD_DIFF = 6.95e-3
    
    # Maximum arc seconds (10 arcmin=600arcsecs aprox.) of spatial distant 
    # allowed between two consecutive frames.
    MAX_RA_DEC_DIFF = 600
    
    # Maximum number or files allowed in a non 'OT' sequence (filter grouped)
    MAX_NFILES = 50
    ############################################################

    def __init__( self , source, instrument=None ):
        """
        Initialize the object.
        
        Parameters
        ----------
        source : str
            Can be a 'directory' name, a 'filename' containing the
            list file or python list havind the files of the DataSet

        instrument: str
            Name of the source instrument of the data files. It must match
            with the FITS keyword 'INSTRUME'; whether the keyword does not
            exist, the file is inserted into DB.
        """
        self.con = None #connection
        self.source = source
        self.id = 0
        
        if instrument!=None:
            self.instrument = instrument.lower()
        else:
            self.instrument = None
            

    ############################################################    
    def createDB (self):
        """
        Create the dataset table
        """
        
        self.con = sqlite.connect(":memory:")
        cur = self.con.cursor()
        cur.execute("create table dataset " + DataSet.TABLE_COLUMNS  )
        self.id = 0

    ############################################################    
    def load( self , source=None ):

        """
        Load the source for files and insert them into the dataset DB
        """

        log.debug("Loading DB ...")
        
        if source==None: source = self.source
        
        # 1. Load the source
        if type(source)==type(list()): contents = source 
        elif os.path.isdir(source):
            log.debug("Loadding Source Directory %s" %source)
            contents = [os.path.join(source, file) for file in os.listdir(source)]
        elif os.path.isfile(source):
            log.debug("Loadding Source File %s" %source)
            contents = [line.replace( "\n", "") for line in fileinput.input(source)]
        else:
            log.error("Error, DB input source not supported !!")
            raise Exception("Error, DB input source not supported")

        # 2. Insert loaded data into in memory DB
        #    -Load and check the FITS file
        #    -Insert into 'dataset' table a new row with data from FITS file
        for file in contents:
            try:
                self.insert( file )
            except Exception,e:
                log.error("Error while inserting file %s " %file)
                #raise
                continue
                        
    ############################################################
    def insert ( self, filename ):
        """
        Insert new FITS file into dateset

        Parameters
        ----------
        filename : str
            input filename to insert into the dataset

        Returns
        -------
        True if all was successful, otherwise False
        """

        log.debug("Inserting file %s into dataset" % filename)
        if filename==None: return False
        
        try:
            if self.GetFileInfo(filename)!=None:
                log.error("File %s not inserted, it is already in Database."%filename)
                return False
        except Exception,e:
            log.exception( "Unexpected error reading FITS file %s" %filename )
            raise e
        
        try:
            fitsf = datahandler.ClFits ( filename )
        except Exception,e:
            log.exception( "Unexpected error reading FITS file %s" %filename )
            raise e
        
        data = (self.id, fitsf.runID, 
                fitsf.obID, fitsf.obPat, fitsf.pat_expno, fitsf.pat_noexp, 
                filename, fitsf.date_obs, fitsf.time_obs, fitsf.mjd, fitsf.type, 
                fitsf.filter, fitsf.exptime, fitsf.ra, fitsf.dec, fitsf.object,
                fitsf.detectorID)
        
        print "dataDB_tuple = ", data
         
        cur = self.con.cursor()

        # Check instrument id
        if (self.instrument==None or 
            (self.instrument==fitsf.getInstrument().lower())):
            try:
                cur.execute("insert into dataset" + DataSet.TABLE_COLUMNS +
                            "values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", data)
                self.con.commit()
    
            except sqlite.DatabaseError,e:
                self.con.rollback()
                log.exception("error inserting into DB")
                raise e
    
            self.id+=1
            return True
            
        else:
            log.error("File %s does not match instrument"%filename)
            return False

    
    ############################################################
    def delete ( self, filename, date=None ):

        """
        Delete a row file from the dataset table

        Parameters
        ----------
        filename: str
            file to be deleted from dataset
        
        Returns
        -------
        0 if all was successful, otherwise <0
          
        """
        
        ###
        if date==None:
            s_date="date>=?"
            date=""
        else:
            s_date="date=?"
        ###    
        try:
            cur = self.con.cursor()
            s_select="delete from dataset where filename=? and %s" %(s_date)
            cur.execute(s_select, (filename, date))
            self.con.commit()
            
        except sqlite.DatabaseError:
            self.con.rollback()
            log.exception("Error in DataSet.delete function")
            raise

        return 0

    ############################################################        
    def clearDB( self ):
        """Clear all DB rows"""
        try:
            cur = self.con.cursor()
            cur.execute("delete from dataset where filename like '%'", (""))
            self.con.commit()
            
        except sqlite.DatabaseError:
            self.con.rollback()
            log.exception("Error in DataSet.clearDB function")
            raise
            
    ############################################################        
    def GetFile( self, detectorId, type, texp, filter, date, runId=None):
        """
        Return the filename/s (with path) from the data frames with specified fields. 
        We can ask for any type of file (calib, science, master calib, reduced, ...)

        Parameters
        ----------
        
        detectorId : str
        type :  str
        texp :  float
        filter :   str
        date : str
        runId : str

        Returns
        -------
        
        A list with the filenames that match the specified fields, otherwise None
        """

        try:
            #The run id may no be specified
            if runId==None or runId=='*':
                s_run_id = "run_id like '%'"
                runId = ""
            else:
                s_run_id = "run_id=?"
            
            #The master flat does not have a 'texp' requirement 
            if type=='MASTER_DOME_FLAT' or type=='MASTER_SKY_FLAT':
                s_texp="texp>? and texp<?"
                ROUND = 10000
            else:
                s_texp = "texp>? and texp<?"
                ROUND = 0.5  # We do not need an accurate value !

            #The master dark does not have a 'filter' requirement
            if type=='MASTER_DARK' or filter=="ANY":
                s_filter = "filter like '%'"
                filter = ""
            else:
                s_filter="filter=?"

            s_select = "select filename from dataset where detector_id=? and  type=? and %s and %s and date=? and %s" %(s_texp,s_filter, s_run_id)
            print s_select
            cur = self.con.cursor()
            #cur.execute("select filename from dataset where detector_id=? and  type=? and texp>? and texp<?  and filter=? and date=? and run_id=?",
            #                 (detectorId, type, texp-ROUND, texp+ROUND, filter, date, runId))
            cur.execute(s_select,(detectorId, type, texp-ROUND, texp+ROUND, filter, date, runId))
            
            rows = cur.fetchall()
            if len(rows)==0:
                # Any match
                return None
            elif len(rows)>1:
                for row in rows:
                    print row #only for debug
                return rows
            else:
                return rows
                
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetFile function...")
            raise
            return None

    ############################################################        
    def GetFiles( self, detectorId='ANY', type='ANY', texp=-1, filter='ANY', 
                  mjd=55000, ra=0, dec=0, delta_pos=360*3600/2, 
                  delta_time=9999999, runId=None):
        """
        Return the filenames (with path) from the data set with specified 
        fields. We can ask for any type of file (calib, science, master calib, 
        reduced, ...). The files are ascending sorted by the MJD. 

        Parameters
        ----------

          detectorId : str
              Detector type
          type : srt
              Image type to look for
          texp : float
              Exposition time to look for
          filter : str
              Image filter used
          mjd : float
              Modified julian date of observation (days)
          ra : float
              Right ascension coordinate
          dec : float
              Declination coordinate
          delta_pos : float
              arcsec for ra/dec search box
          delta_time: float
              secs for data-time obs search time window
          runId : str
          

        Returns
        -------
        A MJD ascending sorted list with the filenames that match the 
        specified fields, otherwise an empty list []
        """

        res_list = []
        try:
            # NOT USED !!! RUN_ID !!! NOT USED 
            #RUNID: The run id may no be specified
            if runId==None or runId=='*':
                s_run_id = "run_id like '%'"
                runId = ""
            else:
                s_run_id = "run_id=?"
            
            #Type: The type of the image (ANY, DARK, FLAT,SCIENCE, ....)
            if type=='ANY':
                s_type="type>=?"
                type = ""
            elif type=='DOME_FLAT':
                s_type="type like 'DOME_FLAT%' or type=?"
                type = ""
            elif type=='SKY_FLAT':
                s_type="type like 'TW_FLAT%' or type like 'SKY_FLAT%' or type=?"
                type = ""
            else:
                s_type="type=?"
                
            #DetectorID: Detector identifier
            if detectorId=='ANY':
                s_detectorId="detector_id>=?"
                detectorId=""
            else:
                s_detectorId="detector_id=?"
                
            #TEXP: Any 'texp' requirement 
            if  texp==-1:
                s_texp = "texp>=%f" %(texp)
            else:
                ROUND = 0.5  # We do not need an accurate value !
                s_texp = "texp>%f and texp<%f" %(texp-ROUND, texp+ROUND)

            #FILTER: The master dark does not have a 'filter' requirement
            if type=='MASTER_DARK' or filter=="ANY":
                s_filter = "filter>=?"
                filter = ""
            else:
                s_filter = "filter=?"

            #AR
            s_ar = "ra>%s and ra<%s" %(ra-delta_pos, ra+delta_pos)
            #DEC
            s_dec = "dec>%s and dec<%s" %(dec-delta_pos, dec+delta_pos)
            #MJD
            s_mjd = "mjd>%s and mjd<%s" %(mjd-delta_time, mjd+delta_time)
            
            s_select = "select filename from dataset where %s and %s and %s and %s and %s and %s and %s order by mjd" %(s_detectorId, s_type, s_filter, s_texp, s_ar, s_dec, s_mjd)
            print s_select
            cur = self.con.cursor()
            #cur.execute("select filename from dataset where detector_id=? and  type=? and texp>? and texp<?  and filter=? and date=? and run_id=?",
            #                 (detectorId, type, texp-ROUND, texp+ROUND, filter, date, runId))
            cur.execute(s_select,(detectorId, type, filter,))
            
            rows = cur.fetchall()
            res_list = []
            if len(rows)>0:
                res_list = [str(f[0]) for f in rows] # important to apply str() !!
            #print "Total rows selected:  %d" %(len(res_list))
            #print "Files found :\n ", res_list
            return res_list
                             
        except sqlite.DatabaseError,e:
            log.exception("Error in DataSet.GetFile function...")
            raise Exception("Error in DataSet.GetFile: %s",str(e))
            return []

    def GetFilesT( self, type, texp=-1, filter="ANY"):
        """ Get all the files which match with the specified type, texp and filter
            If query does not match, then return a empty list []
        """
              
        return self.GetFiles( "ANY", type, texp, filter, mjd=55000 , ra=0, 
                              dec=0, delta_pos=360*3600/2, delta_time=9999999,
                               runId=None)
          
    def GetOBFiles(self, filter=None):
        """ Get all the files for each Observation Block found. 
            filter:  filter can be specified to restrict the search
            
            Return a list of list, having each list the list of files beloging to.
        """
        
        #print "start GetOBFiles....."
        
        ob_id_list=[]
        ob_file_list=[]

        if filter==None:
            s_filter="filter>=?"
            filter=""
        else:
            s_filter="filter=?"
              
        # First, look for OB_IDs
        #s_select="select ob_id from dataset where %s group by ob_id" %(s_filter)
        s_select="select DISTINCT ob_id from dataset where %s" %(s_filter)
        #print s_select
        cur = self.con.cursor()
        cur.execute(s_select,(filter,))
        rows = cur.fetchall()
        if len(rows)>0:
            ob_id_list = [str(f[0]) for f in rows] # important to apply str() !!
        print "Total rows selected:  %d" %(len(ob_id_list))
        print "OB_IDs found :\n ", ob_id_list
        
        # Finally, look for files of each OB_ID
        for ob_id in ob_id_list:
            s_select = "select filename from dataset where ob_id=? order by mjd"    
            #print s_select
            cur = self.con.cursor()
            cur.execute(s_select,(int(ob_id),))
            #print "done !"
            rows = cur.fetchall()
            if len(rows)>0:
                ob_file_list.append([str(f[0]) for f in rows]) # important to apply str() !!
            #print "%d files found in OB %d" %(len(rows), int(ob_id))
            
        return ob_id_list, ob_file_list
    
    def GetFilterFiles(self, max_mjd_diff=None, max_ra_dec_diff=None, 
                       max_nfiles=None):
        """ 
        Get all SCIENCE and CALIB file groups found for each (Filter,Type) 
        ordered by MJD; no other keyword is looked for (OB_ID, OB_PAT, ...).
        
        In addition, MJD is checked in order to look for time gaps into a (Filter,Type) 
        sequence. It will be quite useful for data grouping when the OT was not used 
        during the observing run.
        
        Parameters
        ----------
        max_mjd_diff: float
            Maximum seconds of temporal distant allowed between two consecutive 
            frames.
        
        max_ra_dec_diff: float
            Maximum seconds of spatial distant allowed between two consecutive 
            frames.
        
        max_nfiles: int
            Maximum number of files allowed into a sequence (only for 'filter' grouping).
         
        Returns
        -------
        A List of types [DARK, DOME_FLAT, TW_FLAT, SCIENCE] and list of list, 
        having each list the list of files beloging to.

        Notes
        -----
        It is mainly useful when the OT is not used for the data acquisition
                    
        """
        
        if max_mjd_diff==None: 
            max_mjd_diff = DataSet.MAX_MJD_DIFF
        
        if max_ra_dec_diff==None:
            max_ra_dec_diff = DataSet.MAX_RA_DEC_DIFF
        
        if max_nfiles==None:
            max_nfiles = DataSet.MAX_NFILES
    
        par_list = [] # parameter tuple list (filter, type)
        filter_file_list = [] # list of file list (one per each filter)
              
        # First, look for Filters on SCIENCE files
        #s_select="select DISTINCT filter,texp from dataset where type='SCIENCE' or type='SKY' order by mjd"
        s_select = "select DISTINCT filter,type from dataset where type<>'DOME_FLAT_LAMP_OFF' and type<>'DOME_FLAT_LAMP_ON' order by mjd"
        

        cur = self.con.cursor()
        cur.execute(s_select,"")
        rows = cur.fetchall()
        par_list = []
        if len(rows)>0:
            par_list = [[str(f[0]), str(f[1])] for f in rows] # important to apply str() ??
        print "Total rows selected:  %d"%(len(par_list))
        print "Filters found :\n ", par_list
        
        # Look for DOME_FLATS_ON/OFF
        s_select="select DISTINCT filter from dataset where type='DOME_FLAT_LAMP_OFF' or type='DOME_FLAT_LAMP_ON' order by mjd"
        cur = self.con.cursor()
        cur.execute(s_select,"")
        rows = cur.fetchall()
        par_list2 = []
        if len(rows)>0:
            par_list2 = [[str(f[0]), "DOME_FLAT"] for f in rows] # important to apply str() ??
        print "(2nd) Total rows selected:  %d" %(len(par_list2))
        print "(2nd) Filters found :\n ", par_list2

        # Concatenate the two list (dome_flat, the rest)
        par_list = par_list + par_list2
        
        print "PAR_LIST=",par_list
        # Finally, look for files of each Filter
        for par in par_list:
            cur = self.con.cursor()
            if par[1]!='DOME_FLAT':
                s_select = "select filename from dataset where filter=? and type=? order by mjd"
                cur.execute(s_select,(par[0],par[1]))    
            else:
                s_select = "select filename from dataset where filter=? and type like 'DOME_FLAT%' order by mjd"
                cur.execute(s_select, (par[0],)) # la coma es imprescindible, no se por que .....
            rows = cur.fetchall()
            if len(rows)>0:
                filter_file_list.append([str(f[0]) for f in rows]) # important to apply str() !!
            #print "====> %d files found for Filter %s" %(len(rows), par[0])

        #
        # Now, look for temporal/spatial/size gap inside the current sequences 
        # found.
        #
        new_seq_list = []
        new_seq_par = []
        k = 0
        
        for seq in filter_file_list:
            group = []
            mjd_0 = self.GetFileInfo(seq[0])[10] # reference for temporal gap
            ra_0 = self.GetFileInfo(seq[0])[7]*3600  # reference for spatial gap
            dec_0 = self.GetFileInfo(seq[0])[8]*3600 # reference for spatial gap
            #print "RA_0=",ra_0
            #print "DEC_0=",dec_0
            for file in seq:
                t = self.GetFileInfo(file)[10]
                ra = self.GetFileInfo(file)[7]*3600 # arcsecs
                dec = self.GetFileInfo(file)[8]*3600 #arcsecs
                #print "DIF_RA=",math.fabs(ra-ra_0)
                #print "MaxRADecDiff=",max_ra_dec_diff
                #print "RA0=",ra_0
                #print "RA=",ra
                #print "DIF_DEC=",math.fabs(dec-dec_0)
                #print "LEN_GROUP=",len(group)
                #print "DIF_MJD=",(t-mjd_0)
                #print "MAX_MJD=",max_mjd_diff
                #print "TYPE=",self.GetFileInfo(seq[0])[2]
                    
                # Note: Currently, this code will not distinguish between next
                # dark sequences:
                #     - 2s 2s 2s  5s 5s 5s 10s 10s 10s 20s 20s 20s ...
                #     - 2s 5s 10s 20s ...
                # both will be classified as dark_model sequences, but
                # if max_nfiles=3, then they will be distinguished, although
                # it will also affect the other sequences. 
                if (math.fabs(t-mjd_0) < max_mjd_diff and
                    math.fabs(ra-ra_0) < max_ra_dec_diff and
                    math.fabs(dec-dec_0) < max_ra_dec_diff and
                    len(group) < max_nfiles):
                    group.append(file)
                    mjd_0 = t
                # Darks do not have coordinates restrictions
                elif (self.GetFileInfo(seq[0])[2]=='DARK' and 
                      math.fabs(t-mjd_0)< max_mjd_diff and len(group) < max_nfiles):
                    group.append(file)
                    mjd_0 = t
                # Flats (dome or sky) do not have coordinates restrictions
                elif (self.GetFileInfo(seq[0])[2].find("DOME_FLAT")>=0 and 
                      math.fabs(t-mjd_0)< max_mjd_diff and len(group) < max_nfiles):
                    group.append(file)
                    mjd_0 = t
                else:
                    log.debug("Sequence split due to temporal gap between sequence frames")
                    new_seq_list.append(group[:]) # very important, lists are mutable !
                    new_seq_par.append(par_list[k][1])
                    mjd_0 = t
                    ra_0 = ra
                    dec_0 = dec
                    group = [file]
            
            new_seq_list.append(group[:]) # very important, lists are mutable !
            new_seq_par.append(par_list[k][1]) # add the type of the group
            k+=1    

        return  new_seq_list, new_seq_par
                 
    def GetFilterFiles_TEST(self, max_mjd_diff=None ):
        """ 
        @summary: Get all SCIENCE and CALIB file groups found for each (Filter,Type) 
        ordered by MJD; no other keyword is looked for (OB_ID, OB_PAT, ...).
        
        In addition, MJD is checked in order to look for time gaps into a (Filter,Type) 
        sequence. It will be quite useful for data grouping when the OT was not used 
        during the observing run.
        
        @param max_mjd_diff: Maximum seconds of temporal distant allowed between 
        two consecutive frames
         
        @return: the list of types [DARK, FLAT, etc] and list of list,
        having each list the list of files beloging to.

        @note: it is mainly useful when the OT is not used for the data acquisition
                    
        """
        
        if max_mjd_diff==None: max_mjd_diff=DataSet.MAX_MJD_DIFF
        
        par_list = [] # parameter tuple list (filter,texp)
        filter_file_list = [] # list of file list (one per each filter)
              
        # First, look for Filters on SCIENCE files
        #s_select="select DISTINCT filter,texp from dataset where type='SCIENCE' or type='SKY' order by mjd"
        s_select="select DISTINCT filter,type from dataset  order by mjd"

        #print s_select
        cur = self.con.cursor()
        cur.execute(s_select,"")
        rows = cur.fetchall()
        if len(rows)>0:
            par_list = [ [str(f[0]), str(f[1])]  for f in rows] # important to apply str() ??
        print "Total rows selected:  %d" %(len(par_list))
        print "Filters found :\n ", par_list
        
        # Finally, look for files of each Filter
        for par in par_list:
            s_select = "select filename from dataset where filter=? and type=? order by mjd"
            #s_select = "select filename from dataset where filter=? and texp=?  order by mjd"    
            cur = self.con.cursor()
            cur.execute(s_select,(par[0],par[1]))
            rows = cur.fetchall()
            if len(rows)>0:
                filter_file_list.append([str(f[0]) for f in rows]) # important to apply str() !!
            #print "%d files found for Filter %s" %(len(rows), par[0])

        # Now, look for temporal gap inside the current sequences were found
        new_seq_list = []
        new_seq_par = []
        k = 0
        for seq in filter_file_list:
            group = []
            mjd_0 = self.GetFileInfo(seq[0])[10]
            for file in seq:
                t = self.GetFileInfo(file)[10]
                if math.fabs(t-mjd_0)<max_mjd_diff:
                    group.append(file)
                    mjd_0 = t
                else:
                    log.debug("Sequence split due to temporal gap between sequence frames")
                    new_seq_list.append(group[:]) # very important, lists are mutable !
                    new_seq_par.append(par_list[k][1])
                    mjd_0 = t
                    group = [file]
            new_seq_list.append(group[:]) # very important, lists are mutable !
            new_seq_par.append(par_list[k][1])
            k+=1    
            
        return  new_seq_list, new_seq_par
    
    def GetSciFilterFiles(self, max_mjd_diff=None ):
        """ 
        @summary: Get all the SCIENCE group files found for each (Filter,TExp) 
        ordered by MJD; no other keyword is looked for (OB_ID, OB_PAT, ...).
        
        In addition, MJD is checked in order to look for time gaps into a (Filter,TExp) 
        sequence. It will be quite useful for data grouping when the OT was not used 
        during the observing run.
        
        @param max_mjd_diff: Maximum seconds of temporal distant allowed between 
        two consecutive frames
         
        @return: the list of parameter-tuples [filter, texp] and list of list,
        having each list the list of files beloging to.
        
        @note: currently NOT USED !!!
            
        """
        
        if max_mjd_diff==None: max_mjd_diff=DataSet.MAX_MJD_DIFF
        
        par_list = [] # parameter tuple list (filter,texp)
        filter_file_list = [] # list of file list (one per each filter)
              
        # First, look for Filters on SCIENCE files
        s_select="select DISTINCT filter,texp from dataset where type='SCIENCE' or type='SKY' order by mjd"
        #s_select="select DISTINCT filter,texp from dataset  order by mjd"

        #print s_select
        cur = self.con.cursor()
        cur.execute(s_select,"")
        rows = cur.fetchall()
        if len(rows)>0:
            par_list = [ [str(f[0]),f[1]] for f in rows] # important to apply str() ??
        print "Total rows selected:  %d" %(len(par_list))
        print "Filters found :\n ", par_list
        
        # Finally, look for files of each Filter
        for par in par_list:
            s_select = "select filename from dataset where filter=? and texp=? and (type='SCIENCE' or type='SKY') order by mjd"
            #s_select = "select filename from dataset where filter=? and texp=?  order by mjd"    
            cur = self.con.cursor()
            cur.execute(s_select,(par[0],par[1]))
            rows = cur.fetchall()
            if len(rows)>0:
                filter_file_list.append([str(f[0]) for f in rows]) # important to apply str() !!
            #print "%d files found for Filter %s" %(len(rows), par[0])

        # Now, look for temporal gap inside the current sequences were found
        new_seq_list = []
        new_seq_par = []
        k = 0
        for seq in filter_file_list:
            group = []
            mjd_0 = self.GetFileInfo(seq[0])[10]
            for file in seq:
                t = self.GetFileInfo(file)[10]
                if math.fabs(t-mjd_0)<max_mjd_diff:
                    group.append(file)
                    mjd_0 = t
                else:
                    log.debug("Sequence split due to temporal gap between sequence frames")
                    new_seq_list.append(group[:]) # very important, lists are mutable !
                    new_seq_par.append(par_list[k])
                    mjd_0 = t
                    group = [file]
            new_seq_list.append(group[:]) # very important, lists are mutable !
            new_seq_par.append(par_list[k])
            k+=1    
            
        return  new_seq_list, new_seq_par

                         
    def GetSequences(self, group_by='ot', max_mjd_diff=None, 
                     max_ra_dec_diff=None, max_nfiles=None ):
        """
        General function to look for Sequences in the current data base of files
        
        Parameters
        ----------
        group_by: str
            parameter to decide what kind of data grouping will be done;
            if 'ot', OT keywords will be used, otherwise ('filter'), Filter and TExp will be
            taken into account for the data grouping.
        
        max_mjd_diff: float
            Maximum seconds of temporal distant allowed between two consecutive 
            frames.
        
        max_ra_dec_diff: float
            Maximum seconds of spatial distant allowed between two consecutive 
            frames.
        
        max_nfiles: int
            Maximum number of files allowed into a sequence.
        
        Returns
        -------
        Two lists:
             - a of list, having each list the list of files beloging 
               to the sequence.
             - a list with the Types for each sequence found (DARK, DOME_FLAT, 
             TW_FLAT, SCIENCE)
        """   
        
        if group_by.lower()=='ot': 
            return self.GetSeqFiles()
        elif group_by.lower()=='filter':
            return self.GetFilterFiles(max_mjd_diff, max_ra_dec_diff, max_nfiles)
        elif group_by.lower()=='none':
            seqs = []
            seq_types = []
            seqs = [self.GetFiles()]
            if len(seqs)>0:
                seq_types = [str(self.GetFileInfo(seqs[0][0])[2])]*len(seqs)
            return seqs, seq_types
        else:
            return [],[]
        
    def GetSeqFiles_TEST(self, filter=None, type=None):
        """
        TO BE DEPRECATED ?
         
        Get all the files for each Observing Sequence (OS) found. 
        By OS we mean a set of files with the next common features:
        
            (OB_ID, OB_PAT, FILTER, TEXP)
            
        INPUTS
        @param filter:  filter can be specified to restrict the search
          
        @param type: it can be SCIENCE, DARKs,FLATs  (None=any type)
          
        @note: this method is thought to work fine for SCIENCE sequences;
        for calibration sequences need improvements
        
        @return: two lists:
             - a list of list of parameters of each list found (OB_ID, OB_PAT, FILTER, TEXP)
             - a of list, having each list the list of files beloging 
               to the sequence.
    
        @see: GetSeqFiles()
        
        @deprecated: currently GetSeqFiles() is used, following PAT_EXPN=PAT_NEXP 
        datagrouping rule 
            
        """
        
        seq_pat_list = [] # list of sequences features (ob_id,ob_pat,filter)
        seq_list = [] # list of list of sequences filenames

        if filter==None:
            s_filter = "filter>=?"
            filter=""
        else:
            s_filter = "filter=?"
              
        if type==None:
            s_type = "type>=''"
        elif type=="SCIENCE":
            s_type = "type='SCIENCE' or type='SKY'"
        elif type=="FLAT":
            s_type = "type='SKY_FLAT' or type='DOME_FLAT' or type='FLAT'"
        else:
            s_type = "type='%s'"%(str(type))
                      
        # First, look for OB_IDs
        #s_select="select ob_id,ob_pat,filter from dataset where %s group by ob_id,ob_pat,filter" %(s_filter)
        s_select="select DISTINCT ob_id,ob_pat,filter,texp from dataset where %s and %s order by mjd" %(s_filter,s_type)
        #print s_select
        cur = self.con.cursor()
        cur.execute(s_select,(filter,))
        rows = cur.fetchall()
        if len(rows)>0:
            seq_pat_list = [[f[0], f[1], f[2], f[3]] for f in rows] 
        print "Total rows selected:  %d" %(len(seq_pat_list))
        print "OS's found :\n ", seq_pat_list
        
        # Finally, look for files of each OB_ID
        for seq in seq_pat_list:
            s_select = "select filename from dataset where ob_id=? and ob_pat=? and filter=? and texp=? order by mjd"    
            #print s_select
            cur = self.con.cursor()
            cur.execute(s_select,(seq[0],seq[1],seq[2],seq[3]))
            #print "done !"
            rows = cur.fetchall()
            if len(rows)>0:
                seq_list.append([str(f[0]) for f in rows]) # important to apply str() !!
            #print "%d files found in OS %s" %(len(rows), str(seq[0])+"_"+str(seq[1])+"_"+str(seq[2])+"_"+str(seq[3]))
            
        return seq_pat_list, seq_list
    
    
    def GetSeqFiles(self, filter=None, type=None):
        """ 
        Get all the files for each Observing Sequence (OS) found. 
        By OS we mean a set of files with the next common features:
        
            - start with PAT_EXPN=1 and end with the PAT_EXPN=PAT_NEXP
            
        INPUTS
        @param filter:  filter can be specified to restrict the search
          
        @param type: it can be SCIENCE, DARKs,FLATs  (None=any type)
          
        @note: this method is thought to work fine for SCIENCE sequences;
        for calibration sequences need improvements
        
        @return: two lists:
             - a list of list, having each list the list of files beloging 
               to the sequence.
             - a list with the Types for each sequence found (DARK, DOME_FLAT, 
             TW_FLAT, SCIENCE)
    
        @see: GetSeqFiles()
            
        """
        
        if filter==None:
            s_filter = "filter>=?"
            filter=""
        else:
            s_filter = "filter=?"
              
        if type==None:
            s_type = "type>=''"
        elif type=="SCIENCE":
            s_type = "type='SCIENCE' or type='SKY'"
        elif type=="FLAT":
            s_type = "type='SKY_FLAT' or type='DOME_FLAT' or type='FLAT'"
        else:
            s_type = "type='%s'"%(str(type))
                      
        # First, look for OB_IDs
        #s_select="select ob_id,ob_pat,filter from dataset where %s group by ob_id,ob_pat,filter" %(s_filter)
        s_select="select filename, ob_id, ob_pat, expn, nexp, filter, texp, type \
                from dataset where %s and %s order by mjd" %(s_filter,s_type)
        #print s_select
        cur = self.con.cursor()
        cur.execute(s_select,(filter,))
        rows = cur.fetchall()
        
        found_first = False
        group = []
        seq_list = [] # list of lists of files from each sequence
        seq_types =[] # list of types for each sequence
        for fits in rows:
            print "%s  %s  %s  %s  %s %s  %s  %s  %s"%(fits[0], fits[1], fits[2], # filename, ob_id, ob_pat 
                                       fits[3], fits[4], fits[5], fits[6], fits[7], # expn, nexp, filter, texp, type
                                       fits[3]==fits[4]) # true/false
            if fits[7].count('MASTER'):
                print "--------> Found a MASTER calibration file; it will not be grouped !!!<----------"
                continue
            # Note: if the beginning (fits[3]==1) or the end of a sequence is 
            # not found (fits[3]==fits[4]), then their files (incomplete sequence) 
            # are added to a unknown_group/sequence.
            if fits[3]==1: #expn == 1 ?
                group = [str(fits[0])] # filename
                found_first = True # update flag
                # special case of only-one-file sequences
                if fits[3]==fits[4]:
                    #detected end of the sequence
                    seq_list.append(group[:]) # very important ==> lists are mutable !
                    # Set the 'nice' type
                    if str(fits[7]).count("DOME_FLAT"): my_type = "DOME_FLAT"
                    elif str(fits[7]).count("TW_FLAT"): my_type = "TW_FLAT"
                    else: my_type = str(fits[7])
                    seq_types.append(my_type)
                    group = []
                    found_first = False  # reset flag
            elif found_first: 
                group.append(str(fits[0]))
                if fits[3]==fits[4]:
                    # Detected end of the sequence
                    seq_list.append(group[:]) # very important ==> lists are mutable !
                    # Set the 'nice' type
                    if str(fits[7]).count("DOME_FLAT"): my_type = "DOME_FLAT"
                    elif str(fits[7]).count("TW_FLAT"): my_type = "TW_FLAT"
                    else: my_type = str(fits[7])
                    seq_types.append(my_type)
                    group = []
                    found_first = False  # reset flag
            else:
                pass
        
        #
        # Look for un-groupped files and build a group/sequence with them
        #
        temp = set([])
        for lista in seq_list:
            temp = temp.union(set(lista))
        un_groupped = set(self.GetFiles()) - temp
        if len(un_groupped)>0:
            seq_list.append(list(un_groupped))
            seq_types.append('UNKNOWN')

        #print "[dataset.GetSeqFiles] OS's found :\n ", seq_list
        
        return seq_list, seq_types

    ############################################################

    def GetSeqFilesC(self, filter=None, type=None):
        """
        !!! NOT USED, TBC !!!
        
        NEW VERSION WITH FRAMENUM for distinguish-FITS
        ==============================================
        Get all the files for each Observing Sequence (OS) found. 
        By OS we mean a set of files with the next common features:
        
            - start with PAT_EXPN=1 and end with the PAT_EXPN=PAT_NEXP
            
        Parameters
        ----------
        filter: str
            filter can be specified to restrict the search
          
        type: str
            it can be SCIENCE, DARKs,FLATs  (None=any type)
          
        Notes
        -----
        This method is thought to work fine for SCIENCE sequences;
        for calibration sequences need improvements
        
        Returns
        -------
        Two lists:
             - a list of list, having each list the list of files beloging 
               to the sequence.
             - a list with the Types for each sequence found (DARK, DOME_FLAT, 
             TW_FLAT, SCIENCE)
    
        See also
        --------
        GetSeqFiles()
            
        """
        
        if filter==None:
            s_filter = "filter>=?"
            filter = ""
        else:
            s_filter = "filter=?"
              
        if type==None:
            s_type = "type>=''"
        elif type=="SCIENCE":
            s_type = "type='SCIENCE' or type='SKY'"
        elif type=="FLAT":
            s_type = "type='SKY_FLAT' or type='DOME_FLAT' or type='FLAT'"
        else:
            s_type = "type='%s'"%(str(type))
                      
        # First, look for OB_IDs
        #s_select="select ob_id,ob_pat,filter from dataset where %s group by ob_id,ob_pat,filter" %(s_filter)
        s_select="select filename, ob_id, ob_pat, expn, nexp, filter, texp, type \
                from dataset where %s and %s order by mjd" %(s_filter,s_type)
        #print s_select
        cur = self.con.cursor()
        cur.execute(s_select,(filter,))
        rows = cur.fetchall()
        
        found_first = False
        group = []
        seq_list = [] # list of lists of files from each sequence
        seq_types =[] # list of types for each sequence
        for fits in rows:
            print "%s  %s  %s  %s  %s %s  %s  %s  %s"%(fits[0], fits[1], fits[2], # filename, ob_id, ob_pat 
                                       fits[3], fits[4], fits[5], fits[6], fits[7], # expn, nexp, filter, texp, type
                                       fits[3]==fits[4]) # true/false
            if fits[7].count('MASTER'): 
                print "--------> Found a MASTER calibration file; it will not be grouped !!!<----------"
                continue
            if fits[3]==1 : #and framenum==1: #expn == 1 ?
                group = [str(fits[0])] #filename
                found_first = True # update flag
                # special case of only-one-file sequences
                if fits[3]==fits[4]:
                    #detected end of the sequence
                    seq_list.append(group[:]) # very important ==> lists are mutable !
                    # Set the 'nice' type
                    if str(fits[7]).count("DOME_FLAT"): my_type = "DOME_FLAT"
                    elif str(fits[7]).count("TW_FLAT"): my_type = "TW_FLAT"
                    else: my_type = str(fits[7])
                    seq_types.append(my_type)
                    group = []
                    found_first = False  # reset flag
            elif found_first:
                group.append(str(fits[0]))
                if fits[3]==fits[4]:
                    #detected end of the sequence
                    seq_list.append(group[:]) # very important ==> lists are mutable !
                    # Set the 'nice' type
                    if str(fits[7]).count("DOME_FLAT"): my_type = "DOME_FLAT"
                    elif str(fits[7]).count("TW_FLAT"): my_type = "TW_FLAT"
                    else: my_type = str(fits[7])
                    seq_types.append(my_type)
                    group = []
                    found_first = False  # reset flag
            else:
                pass
        
        print "[dataset.GetSeqFilesC] OS's found :\n ", seq_list
        
        return seq_list, seq_types
                       
                       
    ############################################################    
    def GetFileInfo( self, filename ):
        """
        Query the database fields of a specified filaname.

        Parameteres
        -----------
        filename: str
            filename to query

        Returns
        -------
        A list with some database fields, i.e.:
                      
            date(0), ut_time(1), type(2), filter(3), texp(4), detector_id(5),
            run_id(6), ra(7), dec(8), object(9), mjd(10)
        """

        try:
            # CAUTION: if we modify the query values in the SELECT, it will 
            # affect a lot of code in PAPI !!!! need to be re-writed !!
            s_select = "select date, ut_time, type, filter, texp, detector_id,\
             run_id, ra, dec, object, mjd from dataset where filename=?"
             
            cur = self.con.cursor()
            cur.execute(s_select, (filename,))
            
            rows = cur.fetchall()
            if len(rows)==0:
                # Any match
                return None
            elif len(rows)>1: # it should not happen !!
                for row in rows:
                    print "Two rows were found !!: ", row #only for debug
                return rows[0] # return only the first
            else:
                return rows[0]
                
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetFile function...")
            raise

    ############################################################    
    def GetMasterDark( self, detectorId, texp, date, create=False, runId=None):
        """
        Return the filename (with path) of a master dark with the specified features.

        Parameters
        ----------
        detectorId : str
        texp : float
        date : str
        create : bool 
            If the required master filename does not exist, it will be created whether create=True,
                therwise None will be returned
        runId : str

        Returns
        -------
        A list with the filenames that match the specified fields, otherwise None
        """

        ROUND = 0.5
        type = 'MASTER_DARK'
        try:
            cur=self.con.cursor()
            #The run id may no be specified
            if runId==None or runId=='*':
                cur.execute("select filename from dataset where detector_id=? and  type=? and texp>? and texp<?  and date=? ",
                            (detectorId, type, texp-ROUND, texp+ROUND, date))
            else:
                cur.execute("select filename from dataset where detector_id=? and  type=? and texp>? and texp<?  and date=? and run_id=?",
                            (detectorId, type, texp-ROUND, texp+ROUND, date, runId))
                
            rows=cur.fetchall()
            if len(rows)>0:
                res_list = [str(f[0]) for f in rows] # important to apply str() !!
            #print "Total rows selected:  %d" %(len(res_list))
            #print "Files found :\n ", res_list
            return res_list
                
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetMasterDark function...")
            raise
            return []

    ############################################################    
    def GetMasterFlat( self, detectorId, type, filter, date, create=False, 
                       runId=None):
        """
          \brief TBC
          
          \brief Return the filename (with path) of a master FLAT with the specified features.

          \param detectorId
          \param type  (MASTER_DOME_FLAT, MASTER_SKY_FLAT)
          \param filter
          \param date
          \param runId
          \param create : If the required master filename does not exist, it will be created whether create=True,
          otherwise None will be returned
          
          \return A list with the filenames that match the specified fields, otherwise None
        """

        if (type!='MASTER_DOME_FLAT' and type!='MASTER_SKY_FLAT'):
            log.error("Wrong Master Flat type specified in GetMasterFlat")
            return None

        try:
            cur=self.con.cursor()
            #The run id may no be specified
            if runId==None or runId=='*':
                cur.execute("select filename from dataset where detector_id=? and  type=? and filter=? and date=?",
                            (detectorId, type, filter, date))
            else:
                cur.execute("select filename from dataset where detector_id=? and  type=? and filter=? and date=? and run_id=?",
                            (detectorId, type, filter, date, runId))
                
            rows=cur.fetchall()
            if len(rows)==0 and create==True:
                #Then, we try to compute it
                #TODO
                pass
            elif len(rows)>1:
                for row in rows:
                    print row #only for debug
                return rows
            else:
                return rows
                
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetMasterFlat function...")
            raise
        
            
    ############################################################    
    def ListDataSet( self ):
        """
        List all entries in the dataset

        Return 0 if all was successful, otherwise <0
        """

        FIELD_MAX_WIDTH = 20
        cur = self.con.cursor()
        cur.execute("select * from dataset where type LIKE '%' order by ut_time","")
        
        # Print a header.
        for fieldDesc in cur.description:
            print fieldDesc[0].ljust(FIELD_MAX_WIDTH) ,
        print # Finish the header with a newline.
        print '-' * 120
            
        # For each row, print the value of each field left-justified within
        # the maximum possible width of that field.
        fieldIndices = range(len(cur.description))
        for row in cur:
            for fieldIndex in fieldIndices:
                fieldValue = str(row[fieldIndex])
                print fieldValue.ljust(FIELD_MAX_WIDTH) ,
                    
            print # Finish the row with a newline.
     
    def ListDataSetNames( self ):
        """
        List all entries in the dataset.

        Return 0 if all was successful, otherwise <0
        """

        FIELD_MAX_WIDTH = 20
        cur = self.con.cursor()
        cur.execute("select filename from dataset where type LIKE '%' order by ut_time","")
        
        # Print a header.
        for fieldDesc in cur.description:
            print fieldDesc[0].ljust(FIELD_MAX_WIDTH) ,
        print # Finish the header with a newline.
        print '-' * 120
            
        # For each row, print the value of each field left-justified within
        # the maximum possible width of that field.
        fieldIndices = range(len(cur.description))
        for row in cur:
            for fieldIndex in fieldIndices:
                fieldValue = str(row[fieldIndex])
                print fieldValue.ljust(FIELD_MAX_WIDTH) ,
                    
            print # Finish the row with a newline.
  

################################################################################
#This is only for testing purposes
################################################################################
if __name__ == "__main__": 

    # Get and check command-line options
    usage = "usage: %prog [options]"
    desc = """
This module load into a SQLite databse the files found in the input source. 
"""
    parser = OptionParser(usage, description=desc)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source",
                  help="Source of files to insert into DB (file or directory)")
    
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1 or len(args)!=0:
       parser.print_help()
       sys.exit(0)

    if not options.source:
        parser.print_help()
        parser.error("Incorrect number of arguments " )
    
    if os.path.isfile(options.source):
        try:
            ds = DataSet("file")
            ds.createDB()
            ds.load(options.source)
        except Exception,e:
            log.erro("Error running task %s"%str(e))
            sys.exit(0)
    elif os.path.isdir(options.source):
        try:
            ds = DataSet("directory")
            ds.createDB()
            ds.load(options.source)
        except Exception,e:
            log.erro("Error running task %s"%str(e))
            sys.exit(0)

    ds.ListDataSet()

    sys.exit(0)
