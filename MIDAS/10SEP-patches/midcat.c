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

/*+++++++++++++++++++++ Module MIDCAT +++++++++++++++++++++++++++++++++++++++
.LANGUAGE   C
.IDENTIFICATION  MID* interfaces related to Catalogs
.COMMENTS
holds MID_CCRE, MID_COPN, MID_CCLO + ascii_tst, readcat, rewicat, mergecat
.AUTHOR   Klaus Banse
.KEYWORDS Midas utility routines.
.ENVIRONMENT independent

.VERSION  [1.30]   861110: (K. Banse)

 110303		last modif
------------------------------------------------------------------------*/
 
#include <stdio.h>
#include <string.h>

#include <fileexts.h>
 


/*

*/
 
int MID_CCLO(catno)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  catalog file
.PURPOSE
  close catalog 
.ALGORITHM
  straight forward
.RETURNS stat
-------------------------------------------------------------------------*/
 
int   catno;	/* IN: file id of opened catalog 
		       if = -1, close all catalogs  */

{
int   catid, n;


if (catno < 0) 
   {				/*  CLEAR the complete catalog structure   */
   for (n=0; n<CAT_MAXNO; n++)
      {
      if (CATAL[n].NAME[0] != ' ') 
         {
         CATAL[n].NAME[0] = ' ';
         catid = CATAL[n].FID;
         (void) osaclose(catid);
         }
      }
   }
 
else
   {				/*  CLOSE only specified catalog  */
   if (catno >= CAT_MAXNO)  return (ERR_INPINV);

   if (CATAL[catno].NAME[0] != ' ')
      {
      CATAL[catno].NAME[0] = ' ';
      catid = CATAL[catno].FID;
      if (osaclose(catid) != 0) return (ERR_CATBAD);
      }
   } 
 
return (ERR_NORMAL);
}
 
/*

*/
 
int MID_CKLO(catnam)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  catalog file
.PURPOSE
  close catalog 
.ALGORITHM
  straight forward
.RETURNS stat
-------------------------------------------------------------------------*/
 
char  *catnam;	/* IN: name of catalog to close */

{
int n, catid, status;



status =  ERR_INPINV;

for (n=0; n<CAT_MAXNO; n++)
   {
   if (strcmp(CATAL[n].NAME,catnam) == 0)
      {
      CATAL[n].NAME[0] = ' ';
      catid = CATAL[n].FID;
      if (osaclose(catid) != 0) 
         status = ERR_CATBAD;
      else
         status = ERR_NORMAL;
      }
   } 
 
return (status);
}
 
/*

*/
 
int MID_COPN(catfile,type,cimno)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  catalog file
.PURPOSE
  open catalog file
.ALGORITHM
  straight forward
.RETURNS stat
-------------------------------------------------------------------------*/
 
char	*catfile;	/* IN: catalog file (with \0 )   */
int   *type;	/* OUT: filetype of catalog       */
int   *cimno;	/* OUT: catfile index in CATAL structure  */

{
int  mm, status, n, namtyp, part_no;
int  catid, found;
register int nr, mr;

char	  work[100]; 
register char  cc;
 


*cimno = -1;			/*  get cleaned frame name  */
mm = FNAME_LEN;
status = CGN_CLEANF(catfile,7,work,mm,&namtyp,&part_no);
if (status != 0) return (ERR_FILNAM);                /* invalid syntax */


/*  first search for matching name, then for free slot   */
 
for (n=0; n<CAT_MAXNO; n++)
   {
   if (CATAL[n].NAME[0] != ' ')
      {
      if (strcmp(work,CATAL[n].NAME) == 0)	/* we found the catalog  */
         {
         *type = CATAL[n].TYPE;
         *cimno = n;
         return (ERR_NORMAL);
         }
      }
   }
 
catid = CGN_OPEN(work,READ_WRITE);
if (catid <= 0) 
   {
   catid = CGN_OPEN(work,READ);
   if (catid <= 0) return (ERR_FILNAM);
   }

 
for (n=0; n<CAT_MAXNO; n++)
   {
   if (CATAL[n].NAME[0] == ' ')
      {
      (void) strcpy(CATAL[n].NAME,work);
      found = n;
      goto read_catal;
      }
   }
 
return (ERR_CATOVF);
 
 
read_catal:
status = osaread(catid,work,80);
if (status < 1) goto bad_file;

CATAL[found].TYPE_SET = 1;
n = CGN_INDEXC(work,'=');
if (n < 0) 
   {					/* use image catalog as default */
   CATAL[found].TYPE_SET = 0;
   mm = F_IMA_TYPE;
   (void) strcpy(CATAL[found].DESCR,"IDENT");
   goto next_step;
   }


cc = work[++n];
if ( (cc == 'I') || (cc == 'i') )
   mm = F_IMA_TYPE;
else if ( (cc == 'T') || (cc == 't') )
   mm = F_TBL_TYPE;
else if ( (cc == 'F') || (cc == 'f') )
   mm = F_FIT_TYPE;
else if ( (cc == 'A') || (cc == 'a') )
   mm = F_ASC_TYPE;
else
   goto bad_file;


/* look for ", descr.name */

for (nr=(n+1); ;nr++)		/* start at work[n+1] */
   {
   cc = work[nr];
   if (cc == '\0') break;

   else if (cc == ',')
      {
      mr = nr;
      while(1)
         {
         cc = work[++mr];
         if (cc != ' ')		/* skip any leading blanks */
            {
            if (cc == '\0')
               goto default_dsc;
            else
               {
               (void) strcpy(CATAL[found].DESCR,&work[mr]);
               goto next_step;
               }
            }
         }
      }
   }
default_dsc:
(void) strcpy(CATAL[found].DESCR,"IDENT");
      
next_step:
CATAL[found].FID = catid;
CATAL[found].TYPE = mm;		/*  store type of catalog  */
CATAL[found].RECNO = 1;		/*  next record will be # 1  */
 
*type = mm;
*cimno = found;
return (ERR_NORMAL);
 

bad_file:
osaclose(catid);
CATAL[found].NAME[0] = ' ';
return (ERR_CATBAD);
}

/*

*/

int MID_CCRE(catfile,type,descr,cimno)
 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.KEYWORDS
  catalog file
.PURPOSE
  create up to 3 catalog files in one program
.ALGORITHM
  straight forward
.RETURNS stat
-------------------------------------------------------------------------*/
 
char  *catfile;		/* IN: catalog file (with \0 )   */
int   type;		/* IN: same as filetype in SCFOPN */
char  *descr;		/* IN: char. descr used for Ident field */
int   *cimno;		/* OUT: catfile index in CATAL structure  */

{
int   catid;
int   n, mm, found, namtyp, part_no, status;

char	 work[100]; 

 
/*  get cleaned frame name  */

*cimno = -1;

mm = FNAME_LEN;
CGN_CLEANF(catfile,7,work,mm,&namtyp,&part_no);


/*  first search for matching name, then for free slot   */
 
for (n=0; n<CAT_MAXNO; n++)
   {
   if (CATAL[n].NAME[0] != ' ')
      {
      if (strcmp(work,CATAL[n].NAME) == 0)
         {
         found = n;
         status = osaclose(CATAL[found].FID);
         if (status != 0) return (ERR_CATBAD);
         goto sect_1000;
         }
      }
   }
 
for (n=0; n<CAT_MAXNO; n++)
   {
   if (CATAL[n].NAME[0] == ' ')
      {
      (void) strcpy(CATAL[n].NAME,work);
      found = n;
      goto sect_1000;
      }
   }
 
return (ERR_CATOVF);
 
 
sect_1000:
catid = CGN_OPEN(work,WRITE);
CATAL[found].FID = catid;

if (catid > 0)
   {
   CATAL[found].TYPE =  type;		/*  store type of catalog  */
   mm = (int)strlen(descr);
   if (mm > 47) 
      {
      CATAL[found].NAME[0] = ' ';
      osaclose(catid);
      return (ERR_INPINV);
      }

   (void) strcpy(CATAL[found].DESCR,descr);

   if (type == F_IMA_TYPE) 
      (void) snprintf(work,(size_t)100," =Image catalog, %s",descr);
   else if (type == F_TBL_TYPE) 
      (void) snprintf(work,(size_t)100," =Table catalog, %s",descr);
   else if (type == F_FIT_TYPE) 
      (void) snprintf(work,(size_t)100," =Fit file catalog, %s",descr);
   else
      (void) strcpy(work," =ASCII file catalog");
 
   status = osawrite(catid,work,(int)strlen(work));
   CATAL[found].RECNO = 1;		/*  next record will be # 1  */
   *cimno = found;
   return (ERR_NORMAL);
   }
else
   {
   CATAL[found].NAME[0] = ' ';
   return (ERR_CATBAD);
   }
}

/*

*/
 
int ascii_tst(inpntr,outpntr)
char  *inpntr, *outpntr;

{
char  cbuf[88], ftype[12];

int   reclen, totlen, istat, ascid;
register int  nr;




reclen = CGN_INDEXC(inpntr,'.');
totlen = (int) strlen(inpntr) - 1;
if (*(inpntr+totlen) == ':') return (-9);	/* now come the subdirs */

if ( (reclen <= 0) || ((totlen-reclen) > 8) )
   goto asc_open;


/* isolate file type  */

   (void) strcpy(ftype,inpntr+reclen);
   for (nr=1; nr<10; nr++)
      {					 /* get rid of trailing blanks */
      if (ftype[nr] == ' ')
         {
         ftype[nr] = '\0';
         break;
         }
      }


/* first look for file types which are rejected */

if (strcmp(ftype,".obj") == 0)                  /* VAX/VMS object code */
   return (1);
else if (strcmp(ftype,".o") == 0)               /* Unix object code */
   return (1);
else if (strcmp(ftype,".exe") == 0)             /* VMS + MIDAS executable */
   return (1);
else if (strcmp(ftype,".dvi") == 0)             /* Tex/Dvi file */
   return (1);
else if (strcmp(ftype,".jpg") == 0)             /*  JPEG file */
   return (1);
else if (strcmp(ftype,".bdf") == 0)             /* MIDAS image */
   return (1);
else if (strcmp(ftype,".tbl") == 0)             /* MIDAS table */
   return (1);
else if (strcmp(ftype,".fit") == 0)             /* MIDAS fit file */
   return (1);
else if ( (strcmp(ftype,".fits") == 0) ||
          (strcmp(ftype,".tfits") == 0) ||
          (strcmp(ftype,".mt") == 0) )		/* FITS file */
   return (1);


/* then look for known file types  */

if (strcmp(ftype,".cat") == 0)
   {
   (void) strcpy(outpntr,"MIDAS catalog ^");
   return (0);
   }
else if (strcmp(ftype,".ctx") == 0)
   {
   (void) strcpy(outpntr,"MIDAS context file ^");
   return (0);
   }
else if (strcmp(ftype,".prg") == 0)
   {
   (void) strcpy(outpntr,"MIDAS procedure ^");
   return (0);
   }
else if (strcmp(ftype,".tex") == 0)
   {
   (void) strcpy(outpntr,"TEX file ^");
   return (0);
   }
else if (strcmp(ftype,".inc") == 0)
   {
   (void) strcpy(outpntr,"FORTRAN include file ^");
   return (0);
   }
else if (strcmp(ftype,".h") == 0)
   {
   (void) strcpy(outpntr,"C include file ^");
   return (0);
   }
else if (strcmp(ftype,".hh") == 0)
   {
   (void) strcpy(outpntr,"C++ include file ^");
   return (0);
   }
else if ( (strcmp(ftype,".f") == 0) ||
          (strcmp(ftype,".for") == 0) )
   {
   (void) strcpy(outpntr,"FORTRAN source code ^");
   return (0);
   }
else if (strcmp(ftype,".java") == 0)
   {
   (void) strcpy(outpntr,"JAVA source code ^");
   return (0);
   }
else if (strcmp(ftype,".c") == 0)
   {
   (void) strcpy(outpntr,"C source code ^");
   return (0);
   }
else if ( (strcmp(ftype,".cc") == 0) ||
          (strcmp(ftype,".C") == 0) )
   {
   (void) strcpy(outpntr,"C++ source code ^");
   return (0);
   }
else if (strcmp(ftype,".ps") == 0)
   {
   (void) strcpy(outpntr,"Postscript file ^");
   return (0);
   }
else if (strcmp(ftype,".sh") == 0)
   {
   (void) strcpy(outpntr,"shell script ^");
   return (0);
   }

/*  for all other files we really have to read the first record  */

asc_open:

ascid = CGN_OPEN(inpntr,READ);
if (ascid == -1) return (-1);

istat = 1;

asc_read:
reclen= osaread(ascid,cbuf,20);
if (reclen < 0) goto asc_close;

if (reclen == 0) goto asc_read;

if (0 < cbuf[0])              /* Yes, it's ASCII */
   {
   CGN_UPSTR(cbuf);
   if (strncmp(cbuf,"SIMPLE  =",9) == 0)	/* FITS file */
      istat = 1;				/* don't use */
   else if (strncmp(cbuf,"#!",2) == 0)
      {
      (void) strcpy(outpntr,"shell script ^");
      istat = 0;
      }
   else
      {
      (void) strcpy(outpntr,"ASCII file ^");
      istat = 0;
      }
   }

asc_close:
(void) osaclose(ascid);

return istat;
}

/*

*/

int readcat(id,kno,buf,kk)
int  id;		/* IN: catalog file id */
int  kno;		/* IN: catalog id no. */
int  *kk;		/* OUT: rec_flag = 0, deleted record,
					 = 1, valid record */
char *buf;		/* OUT: buffer with catalog record (if kk = 1) */

/* return actual length of catalog record */


{
int ilen;


*kk = 0;				/* default to non-valid record  */

ilen = osaread(id,buf,CATREC_LEN);
if (ilen > 0) 
   {
   CATAL[kno].RECNO ++;
   if (*buf != '!') *kk = 1;		/* valid record */
   }

return (ilen);	
}

   
   
   

int rewicat(catid,kno)
int  catid;
int  kno;

{
int istat;
char  cbuf[84];


istat = (int) osaseek(catid,0L,FILE_START);
if (istat >= 0) 
   {
   if (CATAL[kno].TYPE_SET == 1)
      istat = osaread(catid,cbuf,CATREC_LEN);	/* skip over type record  */
   CATAL[kno].RECNO = 1;
   }

return (istat);
}

/*

*/

void mergecat(flag,buf1,off,buf2,lbuf2)
int  flag;		/* IN: 1 = Ident buffer,
			       0 = Num buffer */
char *buf1;		/* IN/OUT: buffer with catalog record */
int  *off;		/* IN/OUT: offset of next char. in buf1 */
char *buf2;		/* IN: buffer to be merged into catalog record */
int  lbuf2;		/* IN: length of buf2 */


{
int kk, actlen;
register int nr;



kk = *off;

if (flag == 1)			/* handle Ident field */
   {
   if (lbuf2 >= CATIDENT_LEN)
      {
      (void) strncpy(&buf1[kk],buf2,CATIDENT_LEN);
      kk += CATIDENT_LEN;
      }
   else
      {
      (void) strncpy(&buf1[kk],buf2,(size_t)lbuf2);
      kk += lbuf2;
      nr = CATIDENT_LEN - lbuf2;
      memset(&buf1[kk],32,(size_t)nr);
      kk += nr;
      }
   buf1[kk++] = '^';
   }
else
   {
   actlen = 0;

   for (nr=lbuf2-1; nr>(-1); nr--)
      {
      if (buf2[nr] != ' ') 
         {
         actlen = nr + 1;
         break;
         }
      }

   nr = CATREC_LEN - kk;		/* available space */
   if (actlen > nr) actlen = nr;
   if (actlen > 0)
      {
      (void) strncpy(&buf1[kk],buf2,(size_t)actlen);
      kk += actlen;
      }
   }

*off = kk;
buf1[kk] = '\0';

}
