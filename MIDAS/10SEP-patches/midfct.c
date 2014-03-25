/*===========================================================================
  Copyright (C) 1995-2011 European Southern Observatory (ESO)
 
  This program is free software; you can redistribute it and/or 
  modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public 
  License along with this program; if not, write to the Free 
  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge, 
  MA 02139, USA.
 
  Correspondence concerning ESO-MIDAS should be addressed as follows:
	Internet e-mail: midas@eso.org
	Postal address: European Southern Observatory
			Data Management Division 
			Karl-Schwarzschild-Strasse 2
			D 85748 Garching bei Muenchen 
			GERMANY
===========================================================================*/

/*+++++++++++++++++++++ Module MIDFCT +++++++++++++++++++++++++++++++++++++++
.LANGUAGE   C
.IDENTIFICATION  Module midfct.c
.AUTHOR   Klaus Banse
.KEYWORDS Midas utility routines.
.ENVIRONMENT VMS and UNIX
.COMMENTS
 holds MID_ACCFRM, MID_ACCFITS, 
       MID_CREFRE, MID_FINDFR, MID_INITFR, MID_RETNAM
       MID_FCTCTRL
.VERSION  [1.50] 871110: new repartition KB

 110310		last modif
------------------------------------------------------------------------*/
 
 
#include <fileexts.h>		/* includes all others: e.g. <midas_def.h> */
#include <computer.h>
#include <fsydef.h>
#include <osparms.h>
#include <osyparms.h>
 
#include  <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <libgen.h>
 
 
#define   BIT_0   0x1
 
#define  READ         0
#define  WRITE        1
#define  READ_WRITE   2
#define  APPEND       3
 

static int swapshort = SWAPSHORT;
static int swapint = SWAPINT;
static int floatfmt = FLOATFMT;

/*
from 900606 on: VERS_005	not supported after 04SEP
from 921006 on: VERS_006		"
from 950712 on: VERS_007		"
from 961030 on: VERS_010	not supported after 07SEP
from 010308 on: VERS_100
from 020313 on: VERS_101
from 030120 on: VERS_105
from 060601 on: VERS_110
from 090401 on: VERS_120
*/

static char Midvers[] = "VERS_";

static struct FCT_STRUCT  *fctpntr;

static int compress_flag;

/*

*/
 
int MID_ACCFITS(name,flag,extensio,entrx)

/*++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
frame control table, FITS file access
.PURPOSE
access a FITS file and convert the FITS header to Midas descr. format
.ALGORITHM
Build internal name FITZ`orig-name' and use MID_fitsin to convert 
FITS header to Midas descr. format
.RETURNS
status:    int  return status
--------------------------------------------------*/

char  *name;		/* IN : complete frame name */
int   flag;		/* IN : control flag: 
			   0 - open an existing FITS file,
			   1 - open an existing FITS file again, 
			   2 - create new entry + open (existing) FITS file 
			    for 1. time      */
int   extensio;		/* IN : extension no. */
int   *entrx;		/* OUT: FCT entry number */
 
{
int  kstat, last_compress, first_time, mm;
int  lennam;

char  work[84], bname[160];
char  *myname;

struct FCB_STRUCT *fcbp;




/*
	FITS files are handled directly in Midas as follows:
	The FITS file abc.xyz is read in as any FITS file in Midas, thereby 
	converted to the Midas format, and stored in virtual memory, only.
	This file is internally renamed to FITZabc.xyz if it's a single 
	image/table FITS file.
	FITS extensions are counted with the prime header regarded as 
	extension 0. So, a FITS file with two extensions would have a primary 
	header ("empty", i.e. NAXIS=0, if ext. 1 is a table), and ext. 1 and 
	ext. 2.
	The extension no. is appended to the name, e.g. extension 3 of abc.xyz 
	would have the internal name FITZabc.xyz003. 
	Therefore, FITS files with max. 999 extensions can be processed 
	directly in Midas, currently.
*/
	


/* FITS files must be in current dir. */

mm = CGN_JNDEXC(name,FSY_DIREND);
if (mm > -1)
   {				/* we have a FITS file with a dir. */
   char  curdir[1024];

   (void) getcwd(curdir,(size_t)1024);
   if (curdir[0] != '\0')
      {
      char *dirc, *dirnam;

      dirc = strdup((const char*) name);
      dirnam = dirname(dirc);		/* extract dir. path only */
      kstat = strcmp(curdir,dirnam);
      free(dirc);

      if (kstat == 0)
         {				/* o.k. - file path = current dir. */
         myname = name + mm + 1;	/* take off path */
         goto move_on;
         }
      }

   (void) strcpy(curdir,
     "direct access to FITS files only possible in current working dir...");
   MID_LOG('G',curdir,(int)strlen(curdir));
   printf("%s\n",curdir);			/* avoid using STTPUT() ... */
   return (ERR_INPINV);
   }
else
   myname = name;

move_on:
(void) strcpy(bname,"FITZ");
(void) strcpy(&bname[4],myname);
lennam = (int)strlen(bname);

if (extensio > -1)
   {
   (void) sprintf(work,"%3.3d",extensio);
   (void) strcat(bname,work);
   }

last_compress = compress_flag;	/* save flag from previous ACCFRM call */

kstat = MID_ACCFRM(bname,flag,entrx,&mm);  /* maybe already there... */
if (kstat != ERR_NORMAL)
   {				/* NO - we have to create new frame */
   MID_POPERR();		/* pop errors message from ACCFRM */
   kstat = 				/* convert FITS to Midas */
   MID_fitsin(1,myname,extensio,bname,work,entrx,&mm);
					/* this goes via SCFCRE ... */
   if (kstat != 0) 
      {
      if (kstat == -9)
         return (ERR_FITEXT);
      else
         return (ERR_VERNOR);
      }
   first_time = 1;
   }
else
   first_time = 0;


/* indicate FITS format */

fctpntr = FCT.ENTRIES + (*entrx);
fctpntr->COMPRESS = last_compress;	/* in case ACCFRM compressed... */
fcbp = fctpntr->FZP;
fctpntr->SIZE = fcbp->FITSINF1;	/* use that as nopix ... */
fctpntr->O_NAMLEN = lennam;

if (first_time == 1) 
   fctpntr->FILTYP = mm;		/* FITS format > 0 */

return (ERR_NORMAL);
}
/*

*/
 
int MID_ACCFRM(name,flag,entrx,aux)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
frame control table, access
.PURPOSE
Access the frame control table to retrieve all necessary
information concerning data frames
.ALGORITHM
Search the FCT for an existing frame entry of the same name. 
If a FCT entry does not exist, it is created and a bulk
data file is physically opened. This operation returns
the I/O channel assigned to the file and (for VMS) other useful info.
Finally read in the FCB.
.RETURNS
status:	   int  return status
--------------------------------------------------*/
 
char  *name;	/* IN : complete frame name */
int   flag;	/* IN : control flag: 
			0 - open an existing disk file,
			1 - open an existing disk file again, 
			2 - create new entry + open (existing) disk file 
			    for 1. time      */
int   *entrx;	/* OUT: FCT entry number */
int   *aux;	/* OUT: additional info:
		        [0] = 0 - file was not already in FCT
		        [0] = 1 - yes, file was already in FCT  */
 
{
int   entrx1, n, mm, status;
int   ioch, fid, fidcount;

#if vms
int vmsfida, vmsfidb;
static char  c1='_', c2='G', c3='Z';

#else
static char  c1='.', c2='g', c3='z';
#endif

char  strbuf[8], work[168], *nampntr;

struct FCB_STRUCT  *fcbp;

struct FCT_STRUCT  *fctpntro;

struct LDB_STRUCT  *ldbp;

extern char DATA_PATH[328];




status = ERR_NORMAL;
strbuf[0] = 'M';			/* make sure, it's initialized */
strbuf[1] = '\0';

compress_flag = 0;
mm = *entrx;				/* maybe value should be reused */

if (flag == 2) goto new_entry;


/* flag = 0 or 1 -> first, search for existing entry in FCT */

entrx1 = MID_FINDFR(name);

if (entrx1 < 0)			         /* if not found, test more */
   {
   ioch = (int)strlen(name) - 2;
   if ((name[ioch] == c1) && (name[ioch+1] == 'Z'))
      compress_flag = 1;
   else if ((name[ioch-1] == c1) && (name[ioch] == c2) 
                                 && (name[ioch+1] == c3))
      {
      ioch --;			/* keep pointed to '.'  */
      compress_flag = 2;
      }
   else
      goto new_entry;


   /* we have a compressed frame here */

   name[ioch] = '\0';                     /* name.typ.Z => name.typ */
   entrx1 = MID_FINDFR(name);		/* try again */
   if (entrx1 < 0)	         /* still not found, so do it ... */
      {
      char  *Nullptr = (char *) 0;

#if vms
      if (compress_flag == 1)
         (void) snprintf(work,(size_t) 160,"compress -d %s_Z",name);
      else
         (void) snprintf(work,(size_t) 160,"gzip -d %s_GZ",name);

#else
      if (compress_flag == 1)
         (void) snprintf(work,(size_t) 160,"uncompress %s.Z",name);
      else
         (void) snprintf(work,(size_t) 160,"gzip -d %s.gz",name);
#endif

      oshcmd(work,Nullptr,Nullptr,Nullptr);
      goto new_entry;
      }
   }

fctpntr = FCT.ENTRIES + entrx1;

if (flag == 1)				/* open frame again */
   {
   register int mr;

   fctpntro = fctpntr;			/* save entry of that match */
   entrx1 = MID_CREFRE(name,-1);	/* create one more entry */
   fctpntr = FCT.ENTRIES + entrx1;	/* copy entry */
   fctpntr->IOCHAN = fctpntro->IOCHAN;

#if vms
   fctpntr->FILEIDA = fctpntro->FILEIDA;
   fctpntr->FILEIDB = fctpntro->FILEIDB;
   (void) strcpy(fctpntr->DEVICE,fctpntro->DEVICE);
#else
   fctpntr->FILEID = fctpntro->FILEID;
#endif

   for (mr=0; mr<4; mr++)
      fctpntr->KAUX[mr] = fctpntro->KAUX[mr];
   fctpntr->SIZE = fctpntro->SIZE;
   fctpntr->PROT = fctpntro->PROT;
   fctpntr->COMPRESS = fctpntro->COMPRESS;
   fctpntr->NOBYTE = fctpntro->NOBYTE;
   fctpntr->FORMAT = fctpntro->FORMAT;
   fctpntr->DATTYP = fctpntro->DATTYP;
   fctpntr->PIXPBL = fctpntro->PIXPBL;
   fctpntr->STBLOK = fctpntro->STBLOK;
   fctpntr->FILTYP = fctpntro->FILTYP;
   fctpntr->FITSEXT = fctpntro->FITSEXT;
   fctpntr->LINK[0] = fctpntro->LINK[0];
   fctpntr->LINK[1] = fctpntro->LINK[1];
   fctpntr->CR_FLAG = fctpntro->CR_FLAG;
   fctpntr->O_NAMLEN = fctpntro->O_NAMLEN;
   for (mr=0; mr<3; mr++)
      fctpntr->FITSADDR[mr] = fctpntro->FITSADDR[mr];
   fctpntr->FITSOUT = fctpntro->FITSOUT;
   fctpntr->CATALOG[0] = fctpntro->CATALOG[0];
   fctpntr->CATALOG[1] = fctpntro->CATALOG[1];
   fctpntr->FZP = fctpntro->FZP;
   }	

*aux = 1;			/* indicate that frame already in FCT */
*entrx = entrx1;
return ERR_NORMAL;



/* build up new entry in FCT */

new_entry:				/* create new entry in FCT  */
*aux = 0;			/* indicate that frame is new in FCT */
nampntr = name;
fidcount = 0;
entrx1 = MID_CREFRE(name,mm); 	/* try to create one more entry */

if (unlikely(entrx1 < 0))			/* filename too long */
   {
   MID_ERROR("MIDAS","MID_ACCFRM:",status,0);
   return ERR_FILNAM;
   }

fctpntr = FCT.ENTRIES + entrx1;		/* o.k.  */
fctpntr->COMPRESS = compress_flag;


open_file:				/* open bulk data file */

#if vms
FSY_OPNBDF(nampntr,fctpntr->NAMLEN,&ioch,&vmsfida,&vmsfidb,
           fctpntr->DEVICE,&status);
if ((status & BIT_0) == 0 ) 
   fid = -1;
else
   {
   fid = 0;
   fctpntr->FILEIDA = vmsfida;
   fctpntr->FILEIDB = vmsfidb;
   fctpntr->DEVICE[16] = '\0';
   fctpntr->IOCHAN = ioch;
   status = ERR_NORMAL;
   }

#else
fid = open(nampntr,O_RDWR);
if (fid == -1)		/* if not o.k. for READ+WRITE, try to open READ only */
   {
   fid = open(nampntr,O_RDONLY);
   if (fid > -1) fctpntr->PROT = 2;	/* set read only flag */
   }

fctpntr->FILEID = fid;
fctpntr->IOCHAN = fid;			/* the same in Unix */
#endif


if (fid < 0)
   {
   if (fidcount < 4)
      {
      n = fidcount*80;
      (void) strncpy(work,&DATA_PATH[n],80);
      if (work[0] != '^')
         {
         work[80] = ' ';
         n = CGN_INDEXC(work,' ');
         (void) strcpy(&work[n],name);		/* use new path */
         nampntr = work;
         fidcount ++;
         goto open_file;
         }
      }

   fctpntr->NAME[0] = ' ';                      /* clear entrx1 again */
   free(fctpntr->FZP);
   MID_ERROR("FSY","MID_ACCFRM:",ERR_FRMNAC,0);
   return ERR_FRMNAC;
   }
   

fcbp = fctpntr->FZP;			/* pointer to individual FCB */
if (flag != 2) 
   {			/* get FCB +  check compatibility with host */
   status = OSY_RVB(fctpntr->IOCHAN,(char *)fcbp,OUR_BLOCK_SIZE,1);
   if (unlikely(status != ERR_NORMAL)) 
      {
      status = ERR_FRMNAC;
      (void) strcpy(strbuf,"OSY");
      goto error_return;
      }
    
   if (strncmp(fcbp->VERSION,Midvers,5) == 0)   /* compare only "VERS_" */
      {
      if (swapshort == 12)	
         {
         if (unlikely(fcbp->SWPSHORT != '=')) goto error_return_0;
         }
      else 				/* if (swapshort == 21) */
         {
         if (unlikely(fcbp->SWPSHORT != 's')) goto error_return_0;
         }

      if (swapint == 1234)
         {
         if (unlikely(fcbp->SWPINT != '=')) goto error_return_0;
         }
      else if (swapint == 4321)
         {
      if (unlikely(fcbp->SWPINT != 's')) goto error_return_0;
         }
      else if (swapint == 2143)
         {
         if (unlikely(fcbp->SWPINT != 'h')) goto error_return_0;
         }
      else 				/* if (swapint == 3412) */
         {
         if (unlikely(fcbp->SWPINT != 'w')) goto error_return_0;
         }
 
      if (floatfmt == IEEEFLOAT)
         {
         if (unlikely(fcbp->FLOTFMT != '=')) goto error_return_0;
         }

#if vms
      else if (floatfmt == VAXFLOAT)
         {
         if (fcbp->FLOTFMT != 'V')
            {
            if (fcbp->FLOTFMT != 'G') goto error_return_0;
            SCTPUT("Warning: data is VAX Gfloat instead of Float...");
            }
         }
      else if (floatfmt == VAXGFLOAT)
         {
         if (fcbp->FLOTFMT != 'G')
            {
            if (fcbp->FLOTFMT != 'V') goto error_return_0;
            SCTPUT("Warning: data is VAX Float instead of Gfloat...");
            }
         }
#endif

      else if (floatfmt == HPFLOAT)
         {
         if (unlikely(fcbp->FLOTFMT != 'H')) goto error_return_0;
         }
      }

   else
      {
      status = ERR_VERNOR;		/* version not recognizable... */
      goto error_return;		/* tell caller, to try FITS access */
      }

   fctpntr->SIZE = fcbp->NDVAL;
   fctpntr->NOBYTE = fcbp->NOBYT;
   fctpntr->FORMAT = fcbp->DFORMAT;
   fctpntr->PIXPBL = fcbp->PIXPBL;
   fctpntr->STBLOK = fcbp->D1BLOCK;
   fctpntr->CATALOG[0] = fcbp->BDTYPE[0];

					/*  also read first LDB from disk */
   status = cacheLDB(1,fctpntr->IOCHAN,fcbp->PTRLDB,&ldbp);
   if (unlikely(status != ERR_NORMAL))
      {
      (void) strcpy(strbuf,"MIDAS");
      goto error_return;
      }
   }

fctpntr->CATALOG[1] = 'N';
*entrx = entrx1;
return status;


 
/* here for errors */
 
error_return_0:
 status = ERR_FMTBAD;
 (void) strcpy(strbuf,"MIDAS");
 
error_return:
 if ((char *)fctpntr->FZP != (char *) 0) free ((char *)fctpntr->FZP);
 OSY_DASSGN(entrx1,mm);

 fctpntr->NAME[0] = ' ';		/* clear FCT entry again */
 fctpntr->NAME[1] = '\0';
 MID_ERROR(strbuf,"MID_ACCFRM:",status,0);
 return status;
}	
/*

*/
 
int MID_CREFRE(name,wanted)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  frame control table
.PURPOSE
  create an entry in the frame control table (FCT) for
  specified name + access flag and return the entry no.
  if no room, return status = 0, else status = 1
.ALGORITHM
  straight forward
.RETURNS
  entrx:	I*4		FCT entry no. (>= 0)
  or -2    if filename too long
--------------------------------------------------*/
 
char   *name;		/* IN : complete name of bulk data frame */
int wanted;		/* IN : if in [0,FCT.MAXENT[, it is interpreted \
                        	as wanted index and we try our best at it */
{
int  lna, mm, nr;



#if vms
char upname[FCT_NAME_LEN];
int nll;


nll = FCT_NAME_LEN-1;
CGN_UPCOPY(upname,name,nll);		/* for VMS use upper case names */
#endif




for (nr=0; nr<FCT_NAME_LEN; nr++)		/*  check name length */
   {
   if (name[nr] == '\0') 
      {
      lna = nr;
      goto start_ok;
      }
   }
return (-2); 					/*  name too long!  */


start_ok:
if ((wanted >= 0) && (wanted < FCT.MAXENT))	/* if wanted in valid range */
   {						/* try to fill that entry */
   nr = wanted;
   fctpntr = FCT.ENTRIES + nr;
   if (fctpntr->NAME[0] == ' ') goto do_it;
   }

 
/* look for empty slot in FCT table */

look:
fctpntr = FCT.ENTRIES;
for (nr=0; nr<FCT.MAXENT; nr++)
   {
   if (fctpntr->NAME[0] == ' ') goto do_it;
   fctpntr ++;
   }


/*  no more space in FCT - we have to extend it  */

mm = FCT.MAXENT + 8;			/* increase in chunks of 8 */ 
MID_FCTIN(mm);				/* update also FCT.MAXENT there */
goto look;				/* and search again */

do_it:
#if vms
(void) strcpy(fctpntr->NAME,upname);	/* copy 'upname' + init entry */
#else
(void) strcpy(fctpntr->NAME,name);	/* copy 'name' + init entry */
#endif

fctpntr->NAMLEN = lna;
fctpntr->BDADDR[0] = (char *) 0;
fctpntr->BDADDR[1] = (char *) 0;
fctpntr->PNTR = (char *) 0;
fctpntr->PROT = 3;			/* default to read_write */
fctpntr->COMPRESS = 0;	
fctpntr->CATALOG[1] = 'N';
fctpntr->KAUX[0] = 0;
fctpntr->KAUX[1] = 0;
fctpntr->KAUX[2] = 0;
fctpntr->LINK[0] = 0;
fctpntr->LINK[1] = 0;
fctpntr->ACCESS = 'I';
fctpntr->CR_FLAG = 0;
fctpntr->O_NAMLEN = 0;
fctpntr->DATTYP = D_OLD_FORMAT;
fctpntr->FILTYP = 0;			/* default to Midas file type */
fctpntr->FITSEXT = 0;	
fctpntr->FZP = 
 (struct FCB_STRUCT *) malloc((unsigned int) sizeof(struct FCB_STRUCT));
fctpntr->FITSADDR[0] = (char *) 0;
fctpntr->FITSADDR[1] = (char *) 0;
fctpntr->FITSADDR[2] = (char *) 0;
fctpntr->FITSOUT = ' ';

return (nr);				/* return entry into FCT */
}
/*

*/
 
int MID_FINDFR(name)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  Frame Control Table, bulk data frames
.PURPOSE
  get the entry in the FCT of a bulk data frame
.ALGORITHM
  search the FCT for given name + access flag and return its entry
  if not found return -1
.RETURNS
  entry_no in FCT
-------------------------------------------------------------*/	
 
char  *name;	/* IN : complete name of bulk data frame */
 
{
int  nr;
 

#if vms
char upname[FCT_NAME_LEN];
int nll;

nll = FCT_NAME_LEN-1;
CGN_UPCOPY(upname,name,nll);            /* for VMS use upper case names */
#endif


fctpntr = FCT.ENTRIES;
for (nr=0; nr<FCT.MAXENT; nr++)		/*loop through FCT and compare name */
   {
#if vms
   if ( strcmp(fctpntr->NAME,upname) == 0 ) return nr;
#else
   if ( strcmp(fctpntr->NAME,name) == 0 ) return nr;
#endif

   fctpntr ++;
   }
	
return (-1);
}
/*

*/
 
int MID_FCTCTRL(flag,id,values)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  Frame Control Table, bulk data frames
.PURPOSE
  access/control the FCT contents of a frame from high level code
.ALGORITHM
  function depends on flag
.RETURNS
  if o.k. return 0, else != 0
-------------------------------------------------------------*/	
 
int  flag;	/* IN: action flag */
int  id;	/* IN: FCT entry no.  */
int  *values;	/* IN/OUT: data array, depends on flag */


{
if ((id < 0) || (id >= FCT.MAXENT)) return (-1);
fctpntr = FCT.ENTRIES + id;

if (flag == 1)
   {
   fctpntr->KAUX[2] = *values;
   }

return 0;
}
/*

*/
 
int MID_RETNAM(id,name,lnam)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  Frame Control Table, bulk data frames
.PURPOSE
  get the FCT name of a file id (imno in SCF modules)
.ALGORITHM
  return the name in given FCT_entry
  if internal FITS file, remove leading FITZ and cut off extens. no's
.RETURNS
  if not found return -1, if name buffer not large enough return -2
  else 0
-------------------------------------------------------------*/	
 
int  id;	/* IN: FCT entry no.  */
char  *name;	/* OUT: complete name of bulk data frame */
int  lnam;	/* IN: max. length of name  */

{
int  ilen;



if ((id < 0) || (id >= FCT.MAXENT)) return (-1);


fctpntr = FCT.ENTRIES + id;
if (fctpntr->NAME[0] == ' ') return (-1);	/* file already closed */

if (fctpntr->O_NAMLEN > 0)		/* internal FITS file */
   {					/* FITZorigname... */
   ilen = fctpntr->O_NAMLEN - 4;
   if (lnam <= ilen) return (-2);		/* name too long */

   (void) strcpy(name,fctpntr->NAME+4);
   name[ilen] = '\0';				/* cut off ext no's ... */
   }
else
   {
   if (lnam <= fctpntr->NAMLEN) return (-2);	/* name too long */
   (void) strcpy(name,fctpntr->NAME);
   }

return (0);
}

