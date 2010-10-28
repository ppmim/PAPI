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
from pysqlite2 import dbapi2 as sqlite
import datahandler
import os
import sys

###### Enable logging
from misc.paLog import log


############################################################
class DataSet:  

    """
    \brief
    Class used to define a data set of frames load from a local directory or a Data Base  
    """
    TABLE_COLUMNS="(id, filename, date, ut_time, mjd, type, filter, texp, detector_id, run_id, ra, dec, object)"

    ############################################################
    def __init__( self , _source):
        """
        \brief The constructor
        
        \param source : can be a 'directory' name, a 'filename' containing the list file or a 'db_address'
        """
        self.con=None #connection
        self.source=_source
        self.id = 0

    ############################################################    
    def createDB (self):
        """
        \brief Create the dataset table

        """
        
        self.con = sqlite.connect(":memory:")
        cur = self.con.cursor()
        cur.execute("create table dataset " + DataSet.TABLE_COLUMNS  )
        self.id = 0

    ############################################################    
    def  load( self , source ):

        """
        \brief Load the source for files and insert them into the dataset DB

        """

        # 1. Load the source
        if type(source)==type(list()): contents = source 
        elif os.path.isdir(source):
            log.debug("Loadding Source Directory %s" %source)
            contents = [os.path.join(source, file) for file in os.listdir(source)]
        elif os.path.isfile(source):
            log.debug("Loadding Source File %s" %source)
            contents = [line.replace( "\n", "") for line in fileinput.input(source)]
        else:
            #TODO
            log.error("Error, source not supported !!")
            pass

        # 2. Insert loaded data into in memory DB
        #    -Load and check the FITS file
        #    -Insert into 'dataset' table a new row with data from FITS file
        for file in contents:
            try:
                self.insert( file )
            except:
                log.error("Error while inserting file %s " %file)
                #raise
                continue
                        
    ############################################################
    def insert ( self, filename ):
        """
          \brief Insert new FITS file into dateset

          \param filename : input filename to insert into the dataset

          \return 0 if all was successful, otherwise <0
        """

        log.debug("Inserting file %s into dataset" % filename)
        try:
            fitsf=datahandler.ClFits ( filename )
        except:
            log.exception( "Unexpected error reading FITS file %s" %filename )
            raise
            return -1
        
        data = (self.id, filename, fitsf.date_obs, fitsf.time_obs, fitsf.mjd, fitsf.type, fitsf.filter, fitsf.exptime, fitsf.detectorID, fitsf.runID, fitsf.ra, fitsf.dec, fitsf.object) 
        cur = self.con.cursor()
        try:
            cur.execute("insert into dataset" + DataSet.TABLE_COLUMNS +"values (?,?,?,?,?,?,?,?,?,?,?,?,?)", data)
            self.con.commit()

        except sqlite.DatabaseError:
            self.con.rollback()
            log.exception("error inserting into DB")
            raise

        self.id+=1
        return 0

    
    ############################################################
    def delete ( self, filename, date=None ):

        """
          \brief Delete a row file from the dataset table

          \paran filename file to be deleted from dataset
          \return 0 if all was successful, otherwise <0
          
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
          \brief Return the filename (with path) of a data frame with specified fields. We can ask for
          any type of file (calib, science, master calib, reduced, ...)

          \param detectorId
          \param type
          \param texp
          \param filter
          \param date
          \param runId

          \return A list with the filenames that match the specified fields, otherwise None
        """

        try:
            #The run id may no be specified
            if runId==None or runId=='*':
                s_run_id="run_id like '%'"
                runId=""
            else:
                s_run_id="run_id=?"
            
            #The master flat does not have a 'texp' requirement 
            if type=='MASTER_DOME_FLAT' or type=='MASTER_SKY_FLAT':
                s_texp="texp>? and texp<?"
                ROUND=10000
            else:
                s_texp="texp>? and texp<?"
                ROUND=0.5  # We do not need an accurate value !

            #The master dark does not have a 'filter' requirement
            if type=='MASTER_DARK' or filter=="ANY":
                s_filter="filter like '%'"
                filter=""
            else:
                s_filter="filter=?"

            s_select="select filename from dataset where detector_id=? and  type=? and %s and %s and date=? and %s" %(s_texp,s_filter, s_run_id)
            print s_select
            cur=self.con.cursor()
            #cur.execute("select filename from dataset where detector_id=? and  type=? and texp>? and texp<?  and filter=? and date=? and run_id=?",
            #                 (detectorId, type, texp-ROUND, texp+ROUND, filter, date, runId))
            cur.execute(s_select,(detectorId, type, texp-ROUND, texp+ROUND, filter, date, runId))
            
            rows=cur.fetchall()
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
    def GetFiles( self, detectorId, type, texp, filter, mjd, ra, dec, delta_pos, delta_time, runId=None):
        """
          \brief Return the filenames (with path) of a data set with specified fields. We can ask for
          any type of file (calib, science, master calib, reduced, ...)

          \param detectorId
          \param type
          \param texp
          \param filter
          \param mjd           Modified julian date of observation (days)
          \param ra
          \param dec
          \param delta_pos     arcsec for ra/dec search box
          \param delta_time    secs for data-time obs search time window
          \param runId

          \return A list with the filenames that match the specified fields, otherwise None
        """

        try:
            #RUNID: The run id may no be specified
            if runId==None or runId=='*':
                s_run_id="run_id like '%'"
                runId=""
            else:
                s_run_id="run_id=?"
            
            #DetectorID: Detector identifier
            if detectorId=='ANY':
                s_detectorId="detector_id>=?"
                detectorId=""
            else:
                s_detectorId="detector_id=?"
                
            #TEXP: Any 'texp' requirement 
            if  texp==-1:
                s_texp="texp>%f" %(texp)
            else:
                ROUND=10  # We do not need an accurate value !
                s_texp="texp>%f and texp<%f" %(texp-ROUND, texp+ROUND)

            #FILTER: The master dark does not have a 'filter' requirement
            if type=='MASTER_DARK' or filter=="ANY":
                s_filter="filter>=?"
                filter=""
            else:
                s_filter="filter=?"

            #AR
            s_ar="ra>%s and ra<%s" %(ra-delta_pos, ra+delta_pos)
            #DEC
            s_dec="dec>%s and dec<%s" %(dec-delta_pos, dec+delta_pos)
            #MJD
            s_mjd="mjd>%s and mjd<%s" %(mjd-delta_time, mjd+delta_time)
            
            s_select="select filename from dataset where %s and type=? and %s and %s and %s and %s and %s order by mjd" %(s_detectorId, s_filter, s_texp, s_ar, s_dec, s_mjd)
            print s_select
            cur=self.con.cursor()
            #cur.execute("select filename from dataset where detector_id=? and  type=? and texp>? and texp<?  and filter=? and date=? and run_id=?",
            #                 (detectorId, type, texp-ROUND, texp+ROUND, filter, date, runId))
            cur.execute(s_select,(detectorId, type, filter))
            
            rows=cur.fetchall()
            if len(rows)==0:
                # Any match
                print "Empty response"
                return None
            elif len(rows)>1:
                for row in rows:
                    print row #only for debug
                print "Total rows selected:  %d" %(len(rows))   
                return rows
            else:
                return rows
                         
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetFile function...")
            raise
            return None

    
    ############################################################    
    def GetFileInfo( self, filename ):
        """
          \brief Return the database fields of a specified filaname.

          \param filename

          \return A list with database fields (date, ut_time, type, filter, texp, detector_id, run_id, object)
        """

        try:
            s_select="select date, ut_time, type, filter, texp, detector_id, run_id, ra, dec, object from dataset where filename=?"
            cur=self.con.cursor()
            cur.execute(s_select, (filename,))
            
            rows=cur.fetchall()
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
            return None

     ############################################################    
    def GetMasterDark( self, detectorId, texp, date, create=False, runId=None):
        """
          \brief TBC
          \brief Return the filename (with path) of a master dark with the specified features.

          \param detectorId
          \param texp
          \param date
          \param create : If the required master filename does not exist, it will be created whether create=True,
          otherwise None will be returned
          \param runId

          \return A list with the filenames that match the specified fields, otherwise None
        """

        ROUND=0.5
        type='MASTER_DARK'
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
            if len(rows)==0:
                # No match
                # Then, we try to compute it
                # TODO
                pass
                return None
            elif len(rows)>1:
                for row in rows:
                    print row #only for debug
                return rows
            else:
                return rows
                
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetMasterDark function...")
            raise
            return None

     ############################################################    
    def GetMasterFlat( self, detectorId, type, filter, date, create=False, runId=None):
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
                return None
            elif len(rows)>1:
                for row in rows:
                    print row #only for debug
                return rows
            else:
                return rows
                
        except sqlite.DatabaseError:
            log.exception("Error in DataSet.GetMasterFlat function...")
            raise
            return None
        
     ############################################################    
    def ListDataSet( self ):
        """
           \brief List all entries in the dataset

           \return 0 if all was successful, otherwise <0
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
            \brief List all entries in the dataset

            \return 0 if all was successful, otherwise <0
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


#######################################################################################
# Init Database instance (the same instance for all modules !!)
filesDB=DataSet("directory")
def initDB(source_directory=None):
    filesDB.createDB()
    if source_directory != None:
        filesDB.load(source_directory)
        
    

#######################################################################################
#This is only for testing purposes
#######################################################################################
if __name__ == "__main__": 

    #log.initLog("/tmp/panic.log",logging.DEBUG)
   
    t1=DataSet("directory")
    t1.createDB()
    #t1.load("/disk-a/caha/panic/DEVELOP/PIPELINE/PAPI/simu_panic14F/p1/")
    t1.load("/disk-a/caha/panic/DATA/SIMU_PANIC_3/")
    t1.ListDataSet()
    print 'GET FILES !!!!!!!!!!!!!!!!!!!!!!!! 1/36=%f' %(1/36)
    t1.GetFiles('ANY', 'SCIENCE', -1, 'ANY', 54846.85034991 , 83.9975, -5.240611, 0.0166, 0.00556, runId=0)
    #(detectorId, type, texp, filter, mjd, ar, dec, delta_pos, delta_time, runId=None):
    #print "\n--------> Now, we delete file orion0021_x4.fits " 
    
    sys.exit()
    
    t1.delete("/disk-a/caha/panic/DATA/data_mat/QL1/orion0021_x4.fits","2007-10-18")
    t1.ListDataSet()

    print "\n--------> Now, insert the file orion0021_x4.fits "
    t1.insert("/disk-a/caha/panic/DATA/data_mat/QL1/orion0021_x4.fits")
    t1.ListDataSet()

    print "\n--------> Test for GetFile........"
    r=t1.GetFile("O2k", "SCIENCE", 48, 'KS', '2007-10-18','0')
    print "Result from GetFile : %s" %r 
    
    
    

