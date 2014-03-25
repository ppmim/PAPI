/*===========================================================================
  Copyright (C) 1985-2011 European Southern Observatory (ESO)
 
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

/*+++++++++++++++++++++   Module SCC   ++++++++++++++++++++++++++++++++++++
.LANGUAGE   C
.IDENTIFICATION  Module SCC.C
.COMMENTS
Holds catalog interfaces:
SCCAD, SCCSUB, SCCCRE, SCCCRA, create_cat
.AUTHOR   Klaus Banse	ESO - Garching
.KEYWORDS MIDAS catalogues
.ENVIRONMENT independant

.VERSION  [1.00]  870911: creation from FORTRAN version 3.55 as of 850919

 110303		last modif
------------------------------------------------------------------------*/
 
#include <fileexts.h>
#include <fsydef.h>

#include <stdio.h>
#include <string.h>

/*

*/
 
int SCCSUB(catfile,name)
 
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PURPOSE
  remove given entry from catalog
.ALGORITHM
  search for given entry + change first 2 chars to //
.RETURNS
  int return status
------------------------------------------------------------------------- */
 
char    *catfile;	/* IN: catalog file  */
char	*name;		/* IN: name of frame to be removed  */
 
{
char output[84];
char catrec[CATREC_LEN+4], catrc[CATREC_LEN+4];
 
int   catid, catlen, loop_cnt;
int   n, ll, stat, cattyp;
int   saveoff, cimno, valid;

 
stat = MID_COPN(catfile,&cattyp,&cimno);
if (stat != ERR_NORMAL) 
   {
   if (stat == ERR_FILNAM)
      SCTPUT("(ERR) SCCSUB: - FILNAM");
   else
      SCTPUT("(ERR) SCCSUB: - CATOVF");
   return stat;
   }

catid = CATAL[cimno].FID;


/*  position record pointer to beginning of catalog file  */

if (CATAL[cimno].RECNO > 1)
   {
   stat = rewicat(catid,cimno);		/* that also skips over type record  */
   if (stat < 0) goto osa_error;
   }

(void)strcpy(output,name);

#if vms
CGN_LOWSTR(output);
#endif

n = CGN_INDEXC(output,' ');		/* cutoff trailing blanks */
if (n > 0) output[n] = '\0';

n = CGN_JNDEXC(output,FSY_DIREND);		/* search backwards */
ll = CGN_JNDEXC(output,FSY_TYPMARK);		/* search backwards for type */

if (ll <= n) 
   (void)strcat(output,FSY_DEFPNTR[cattyp-1]);	/* append default file type  */

loop_cnt = 0;


/*  loop through catalog + search for given name */

read_loop:
catlen = readcat(catid,cimno,catrec,&valid);
if (catlen < 0) goto not_found;			/*  EOF reached  */

loop_cnt ++;
if (valid == 0) goto read_loop;
 

/* compare name with name_field of catal record  */

ll = CGN_INDEXC(catrec,' ');		/* find end of name */
if (ll < 1)
   {
   (void) printf("SCCSUB: no file delimiter...\n");
   ll = 1;
   }
(void) strncpy(catrc,catrec,ll);
catrc[ll] = '\0';

if (strcmp(catrc,output) == 0) 
   {
   ll = CGN_COPY(catrc,catrec);			/* save record */

#if vms
   stat = (int) osaseek(catid,0L,FILE_START);		/* rewind file */
   if (stat < 0) goto osa_error;

   for (n=0; n<loop_cnt; n++)		/* read previous records again */
      (void) osaread(catid,catrec,CATREC_LEN);	/* including catal_type_rec */

   stat = (int) osaseek(catid,0L,FILE_CURRENT);

#else
   saveoff = (int) osaseek(catid,0L,FILE_CURRENT);   /* get current position */
   if (saveoff < 0) goto osa_error;
   saveoff -= (ll + 1);				/* move to beginning of line */

   stat = (int) osaseek(catid,0L,FILE_START);		/* rewind file */
   if (stat < 0) goto osa_error;

   stat = (int) osaseek(catid,(long)saveoff,FILE_START); 
#endif

   if (stat < 0) goto osa_error;

   catrc[0] = '!';			/* mark as deleted  */
   catrc[1] = ' ';
   stat = osawrite(catid,catrc,catlen);
   if (stat < catlen) 
      goto osa_error;
   else
      goto all_ok;
   }
goto read_loop;

all_ok:
stat = MID_CCLO(cimno);
return stat;

not_found:
stat = ERR_INPINV;
SCTPUT("(ERR) SCCSUB: - INPINV");
return stat;

osa_error:
stat = ERR_CATBAD;
SCTPUT("(ERR) SCCSUB: - CATBAD");
return stat;
}

/*

*/
 
int SCCCRE(catfile,type,flag)
 
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PURPOSE
  create a catalog file
.ALGORITHM
  1. read each filename from directory output file `directory.dat'
  2. read relevant descriptors from that file and store in catalog file
.RETURNS
  int return status
------------------------------------------------------------------------- */
 
char	*catfile;	/* IN: catalog file  (with \0)  */
int    type;	/* IN: type of files to use, same as in SCFCRE  */
int    flag;	/* IN: = 1, if file 'dirfile.ascii' exists */

{
  
/*  use SCCCRA with descr IDENT */

return (SCCCRA(catfile,type,flag,"IDENT"));

}

/*

*/
 
int SCCCRA(catfile,type,flag,catdsc)
 
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PURPOSE
  create a catalog file
.ALGORITHM
  1. read each filename from directory output file `directory.dat'
  2. read relevant descriptors from that file and store in catalog file
.RETURNS
  int return status
------------------------------------------------------------------------- */
 
char   *catfile;	/* IN: catalog file  (with \0)  */
int    type;		/* IN: type of files to use, same as in SCFCRE  */
int    flag;		/* IN: = 1, if file 'dirfile.ascii' exists */
char   *catdsc;		/* IN: char. descr. used for Ident field */

{
int   stat, mycno;

  
/*  get catalog file access  */
 
stat = MID_CCRE(catfile,type,catdsc,&mycno);
if (stat != ERR_NORMAL)
   {
   if (stat = ERR_INPINV)
      SCTPUT("(ERR) SCCCRE: - INPINV");
   else if (stat = ERR_CATBAD)
      SCTPUT("(ERR) SCCCRE: - CATBAD");
   else 
      SCTPUT("(ERR) SCCCRE: - CATOVF");
   return stat;
   }
else 
   return (create_cat(catfile,type,flag,mycno));

}

/*

*/
 
int create_cat(catfile,type,flag,cno)
 
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PURPOSE
  actual routine to create a catalog file
.ALGORITHM
  1. read each filename from directory output file `directory.dat'
  2. read relevant descriptors from that file and store in catalog file
.RETURNS
  int return status
------------------------------------------------------------------------- */
 
char	*catfile;	/* IN: catalog file  (with \0)  */
int    type;	/* IN: type of files to use, same as in SCFCRE  */
int    flag;	/* IN: = 1, if file 'dirfile.ascii' exists */
int    cno;	/* IN: index into catalog structure */

{
char output[84], *catdsc, cbuf[CATIDENT_LEN+4];
char catrec[CATREC_LEN+4];

int   e_c, e_l, e_d, deja, typflg;
int   kcount, iav, nullo, ibuf[5];
int   ival[3], nax, imno=-1;
int   catid, uni;
int   did, kk, n, stat;
int   catoff, flag_FITS;
 
void  mergecat();


  
catid = CATAL[cno].FID;
catdsc = CATAL[cno].DESCR;			/* point to catalog descr. */
kcount = 0;
if (flag != 1)	goto catalog_close;		/* nothing to be done...*/
 
	
/*  open directory file 'dirfile.ascii'  */

did = CGN_OPEN("dirfile.ascii",READ);
if (did == -1)
   {
   SCTPUT("No file `dirfile.ascii' found...");
   flag = 0;
   goto catalog_close;
   }

e_c = ERRO_CONT;			/* save error continuation flag  */
e_l = ERRO_LOG;				/* and error log flag  */
e_d = ERRO_DISP;			/* and display flag  */
ERRO_CONT = 1;				/* continue on errors  */
ERRO_LOG = 0;				/* and log no errors  */
ERRO_DISP = 0;				/* and display no errors...  */
	
	
read_loop:
stat = osaread(did,output,60);		/* max. 60 chars filenames... */
if (stat < 0) 
   goto end_of_file;		/*  EOF reached  */
else if (stat == 0) 
   goto read_loop;

 
#if vms
n = CGN_JNDEXC(output,FSY_TESTA);		/* remove version no. */
if (n > 0) output[n] = '\0';
CGN_LOWSTR(output);    
#endif
	
n = CGN_JNDEXC(output,FSY_DIREND);		/* search backwards */
kk = CGN_JNDEXC(output,FSY_TYPMARK);		/* search backwards */
if ((type != F_ASC_TYPE) && (kk <= n))
   {
   char  newrec[160];

   (void) sprintf(newrec,
          "\"%s\" not added to catalog - missing type...",output);
   SCTPUT(newrec);
   goto read_loop;
   }

if (strncmp(output,"middumm",7) == 0)
   goto read_loop;			/* avoid middummX frames ... */

 
/*  fill name field of catalog record  */

memset((void *)catrec,32,(size_t)CATREC_LEN);
catrec[CATREC_LEN] = '\0';
catoff = CGN_COPY(catrec,output);
catrec[catoff++] = ' ';			/* remove the `\0'  again  */
 
 
/* handle ASCII files specially...  */
 
if (type == F_ASC_TYPE)
   {
   iav = ascii_tst(output,&catrec[catoff]);	/* check, what ASCII file */
   if (iav == 0)
      {
      kcount ++;
      stat = osawrite(catid,catrec,(int)strlen(catrec));
      CATAL[cno].RECNO ++;
      }
   else
      {
      char  newrec[160];

      if (iav == -9)			/* we reached a subdirectory... */
         goto end_of_file;

      (void) sprintf(newrec,
             "\"%s\"  no text file, omitted ...",output);
      SCTPUT(newrec);
      }
   goto read_loop;
   }
 
 
/*  look for relevant descriptors   */
	
n = SCFINF(output,0,ibuf);
if (n == ERR_NORMAL)			/* see, if already opened before */
   deja = 0;
else
   deja = -1;

flag_FITS = 0;				/* default to no FITS format */
ibuf[1] = -999;
stat = SCFINF(output,9,ibuf);			/* get type of frame */
if (stat != ERR_NORMAL) 
   {
   char  newrec[160];

   (void) sprintf(newrec,"Warning: Could not open file %s ...",output);
   SCTPUT(newrec);
   goto read_loop;
   }

kcount ++;

if (type != ibuf[1])
   { 
   char  newrec[160];

   (void) sprintf(newrec,
          "Warning: File %s not of same type as catalog ...",output);
   SCTPUT(newrec);

   if (ibuf[1] == F_IMA_TYPE)
      {
      (void) SCFOPN(output,D_OLD_FORMAT,0,F_IMA_TYPE,&imno);
      typflg = -1;
      }
   else if (ibuf[1] == F_TBL_TYPE)
      {
      (void) SCFOPN(output,D_OLD_FORMAT,0,F_TBL_TYPE,&imno);
      typflg = -2;
      flag_FITS = ibuf[0];	/* save TBL (0) or FITS (>0) format */
      }
   else                                      /* only FIT files left */
      {
      (void) SCFOPN(output,D_OLD_FORMAT,0,F_FIT_TYPE,&imno);
      typflg = -3;
      }
   }

else
   {
   typflg = 0;
   (void) SCFOPN(output,D_OLD_FORMAT,0,type,&imno);	/* SCFINF was o.k. */
   if (type == F_TBL_TYPE) 
      flag_FITS = ibuf[0];	/* save TBL (0) or FITS (>0) format */
   }

stat = SCDGETC(imno,catdsc,1,CATIDENT_LEN,&iav,cbuf);
if (stat != ERR_NORMAL)  	
   {
   if (typflg == 0)                         /* matching frame type */
      iav = CGN_COPY(cbuf,"   ");
   else if (typflg == -1)
      iav = CGN_COPY(cbuf,"is image");
   else if (typflg == -2)
      iav = CGN_COPY(cbuf,"is table");
   else                                      /* only typflg = -3 possible */
      iav = CGN_COPY(cbuf,"is Fit_file");
   }


/* merge CATAL.DESCR buffer into catalog record */
 
mergecat(1,catrec,&catoff,cbuf,iav);

if (typflg == 0)			/* for matching file types */
   {					/* get auxiliary info */
   if (type == F_IMA_TYPE)			/*  handle images */
      {
      stat = SCDRDI(imno,"NAXIS",1,1,&iav,&nax,&uni,&nullo);
      if (stat == ERR_NORMAL)  	
         {
         kk = nax;
         if (kk > 3) kk = 3;			/* only space for three Npix */
         stat = SCDRDI(imno,"NPIX",1,kk,&iav,ival,&uni,&nullo);
         if (stat == ERR_NORMAL)  	
            {
            if (nax == 1)
               (void)sprintf(cbuf,"%d %d",nax,ival[0]);
            else if (nax == 2)
               (void)sprintf(cbuf,"%d %d,%d",nax,ival[0],ival[1]);
            else
               (void)sprintf(cbuf,"%d %d,%d,%d",nax,ival[0],ival[1],ival[2]);

            iav = (int) strlen(cbuf);
            mergecat(0,catrec,&catoff,cbuf,iav);	/* merge Num field */
            }
         else
            (void) sprintf(cbuf,"NPIX");
         }
      else
         (void) sprintf(cbuf,"NAXIS");
      }

   else if (type == F_TBL_TYPE)		/*  and tables  */
      {
      stat = SCDRDI(imno,"TBLCONTR",3,2,&iav,ival,&uni,&nullo);
      if (stat == ERR_NORMAL)  	
         {
         (void)sprintf(cbuf," %5d %5d",ival[0],ival[1]);
         iav = (int) strlen(cbuf);
         mergecat(0,catrec,&catoff,cbuf,iav);	/* merge Num field */
         }
      else
         (void) sprintf(cbuf,"TBLCONTR");
      }

   if (stat != ERR_NORMAL)
      {
      char  newrec[160];

      (void) sprintf(newrec,
             "Warning: descr %s of %s is corrupted...",cbuf,output);
      SCTPUT(newrec);
      }
   } 


/*  clear up, but release frame only if it was not opened before... */

if (deja == -1) 
   {
   if (flag_FITS > 0)
      stat = TCTCLO(imno);
   else
      stat = SCFCLO(imno);
   if (stat != ERR_NORMAL)
      {
      char  newrec[160];

      (void) sprintf(newrec,
             "Warning: could not close correctly file: %s ...",output);
      SCTPUT(newrec);
      }
   }
	
stat = osawrite(catid,catrec,catoff);	/*  write new record	*/
CATAL[cno].RECNO ++;

goto read_loop;				/* look for next file */
 
	
/*  come here at EOF - and reset keywords  */
 
end_of_file:
osaclose(did);
ERRO_CONT = e_c;
ERRO_LOG = e_l;
ERRO_DISP = e_d;

catalog_close:
stat = MID_CCLO(cno);
(void) SCKWRI("OUTPUTI",&kcount,10,1,&uni);
return stat;
}
 
/*

*/

int SCCADD(catfile,name,ident)
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PURPOSE
  add catalog entry
.ALGORITHM
  add an entry with given name + identification field to the catalog
.RETURNS
  int return status
------------------------------------------------------------------------- */
 
char    *catfile;	/* IN: catalog file  */
char	*name;		/* IN: name of frame to be added  */
char	*ident;		/* IN: Ident info of frame to be added  */
 
{
char  output[200], *catdsc, cbuf[80];
char  catrec[CATREC_LEN+4], catrc[CATREC_LEN+4];
char  *osmsg();

int   ibuf[5], iav, nullo, nax, ival[3], e_c, e_l;
int   uni, deja;
int   imno=-1;
int   catid, catoff, catlen, n, stat, kk;
int   cattyp, typflg, saverec, cimno, valid, noblank;

void mergecat();



stat = MID_COPN(catfile,&cattyp,&cimno);
if (stat != ERR_NORMAL) 
   {
   if (stat == ERR_FILNAM)
      SCTPUT("(ERR) SCCADD: - FILNAM");
   else
      SCTPUT("(ERR) SCCADD: - CATOVF");
   return stat;
   }

catid = CATAL[cimno].FID;
catdsc = CATAL[cimno].DESCR;


/*  position record pointer to begin of catalog file  */

if (CATAL[cimno].RECNO > 1)
   {
   stat = rewicat(catid,cimno);		/* that also skips over type record  */
   if (stat < 0) goto osa_error;
   }

(void)strcpy(output,name);

#if vms
CGN_LOWSTR(output);
#endif

n = CGN_INDEXC(output,' ');             /* cutoff trailing blanks */
if (n > 0) output[n] = '\0';

n = CGN_JNDEXC(output,FSY_DIREND);              /* search backwards */
kk = CGN_JNDEXC(output,FSY_TYPMARK);            /* search backwards for type */

if (kk <= n)
   (void)strcat(output,FSY_DEFPNTR[cattyp-1]);  /* append default file type  */

if (strncmp(output,"middumm",7) == 0)          /*  avoid middummX frames ... */
   {
   char  newrec[160];

   (void) sprintf(newrec,
          "Warning: dummy file %s not stored in catalog...",output);
   SCTPUT(newrec);
   return (ERR_NORMAL);
   }


/* check that file is of correct type  */

typflg = 0;				/* init to success */
if (cattyp != F_ASC_TYPE)
   {
   n = SCFINF(output,0,ibuf);
   if (n == ERR_NORMAL)			/* see, if already opened before */
      deja = 0;
   else
      deja = -1;

   ibuf[1] = -999;
   stat = SCFINF(output,1,ibuf);		/* get type of frame */
   if (stat != ERR_NORMAL) 
      {
      char  newrec[160];

      (void) sprintf(newrec,
             "Could not open file %s ",output);
      SCTPUT(newrec);
      SCTPUT("(ERR) SCCADD: - INPINV");
      return stat;
      }

   else if (cattyp != ibuf[1])
      {
      char  newrec[160];

      (void) sprintf(newrec,
             "Warning: File %s not of same type as catalog ...",output);
      SCTPUT(newrec);

      if (ibuf[1] == F_IMA_TYPE)
         {
         (void) SCFOPN(output,D_OLD_FORMAT,0,F_IMA_TYPE,&imno);
         typflg = -1;
         }
      else if (ibuf[1] == F_TBL_TYPE)
         {
         (void) SCFOPN(output,D_OLD_FORMAT,0,F_TBL_TYPE,&imno);
         typflg = -2;
         }
      else					/* only FIT files left */
         {
         (void) SCFOPN(output,D_OLD_FORMAT,0,F_FIT_TYPE,&imno);
         typflg = -3;
         }
      }
   }

else
   {
   deja = 0;				/* avoid SCFCLO later on */
   kk = ascii_tst(output,catrec);	/* catrec is overwritten later on */
   if (kk != 0)
      {
      char  newrec[160];

      (void) sprintf(newrec,
             "Warning: File %s not an ASCII file ...",output);
      SCTPUT(newrec);
      iav = CGN_COPY(cbuf,"no ASCII file");
      typflg = 1;
      }
   }


/*  copy name into 1. field of catalog record  */

memset((void *)catrec,32,(size_t)CATREC_LEN);
catrec[CATREC_LEN] = '\0';
catoff = CGN_COPY(catrec,output);
catrec[catoff++] = ' ';                 /* remove the \0  again  */

if (typflg == 1)
   {					/* we don't have Ident + other info */
   mergecat(1,catrec,&catoff,cbuf,iav);		/* merge Ident field */
   mergecat(0,catrec,&catoff," ",1);		/* merge Num field */
   goto find_loop;				/* now store it in catalog */
   }


/* work on catalog descriptor */

e_c = ERRO_CONT;			/* save error continuation flag  */
e_l = ERRO_LOG;			/* and error log flag...  */
ERRO_CONT = 1;			/* continue on errors  */
ERRO_LOG = 0;			/* and log no errors...  */

iav = (int)strlen(ident);		/* is it just blank? */
noblank = 0;
for (n=0; n<iav; n++)
   {
   if (ident[n] != ' ') 
      {
      noblank = 1;
      break;
      }
   }

if (typflg == 0)				/* SCFINF was o.k. */
   (void) SCFOPN(output,D_OLD_FORMAT,0,CATAL[cimno].TYPE,&imno);

if (noblank != 0)
   {
   if (iav > CATIDENT_LEN) iav = CATIDENT_LEN;
   (void) strncpy(cbuf,ident,iav);
   }

else
   {
   stat = SCDGETC(imno,catdsc,1,CATIDENT_LEN,&iav,cbuf);
   if (stat != ERR_NORMAL)
      {
      if (typflg == 0)				/* matching frame type */
         iav = CGN_COPY(cbuf,"   ");
      else if (typflg == -1)
         iav = CGN_COPY(cbuf,"is image");
      else if (typflg == -2)
         iav = CGN_COPY(cbuf,"is table");
      else 	
         iav = CGN_COPY(cbuf,"is Fit_file");
      }
   }
                                /* merge Ident buffer into catalog record */
mergecat(1,catrec,&catoff,cbuf,iav);

if (typflg == 0)			/* for matching file types */
   {					/* get auxiliary info */ 
   if (CATAL[cimno].TYPE == F_IMA_TYPE)		/* handle images */
      {
      nax = -1;
      (void) SCDRDI(imno,"NAXIS",1,1,&iav,&nax,&uni,&nullo);
      if (nax < 1)		/* empty header file */
         iav = sprintf(cbuf,"%d",nax);
      else
         {
         ival[0] = ival[1] = ival[2] = -1;
         if (nax > 3) 
            kk = 3;          /* only space for three Npix */
         else
            kk = nax;
         stat = SCDRDI(imno,"NPIX",1,kk,&iav,ival,&uni,&nullo);
         if (nax == 1)
            iav = sprintf(cbuf,"%d %d",nax,ival[0]);
         else if (nax == 2)
            iav = sprintf(cbuf,"%d %d,%d",nax,ival[0],ival[1]);
         else
            iav = sprintf
                  (cbuf,"%d %d,%d,%d",nax,ival[0],ival[1],ival[2]);
         }
      mergecat(0,catrec,&catoff,cbuf,iav);   /* merge Num field */
      }


   else if (CATAL[cimno].TYPE == F_TBL_TYPE)		/* and tables  */
      {
      stat = SCDRDI(imno,"TBLCONTR",3,2,&iav,ival,&uni,&nullo);
      if (stat == ERR_NORMAL)
         {
         (void)sprintf(cbuf," %5d %5d",ival[0],ival[1]);
         iav = (int) strlen(cbuf);
         mergecat(0,catrec,&catoff,cbuf,iav);      /* merge Num field */
         }
      else
         (void) sprintf(cbuf,"TBLCONTR");
      }

   if (stat != ERR_NORMAL)
      {
      char  newrec[160];

      (void) sprintf(newrec,
             "Warning: descr %s of %s is corrupted...",cbuf,output);
      SCTPUT(newrec);
      }
   }

/*  clear up, but release frame only if it was not opened before... */
 
if (deja == -1) (void) SCFCLO(imno);

ERRO_CONT = e_c;
ERRO_LOG = e_l;


/* search through complete file, if entry already there ...   */

find_loop:
catlen = readcat(catid,cimno,catrc,&valid);
if (catlen < 0)
   {                                            /* EOF reached */
   stat = (int) osaseek(catid,0L,FILE_END);     /* position write pointer */
   if (stat < 0)
      goto osa_error;
   else
      goto not_found;
   }
else if (valid == 0)			/* was a deleted record */
   goto find_loop;


/* compare name with name_field of catal record  */

n = CGN_INDEXC(catrc,' ');            /* find end of name */
if (n < 1)
   {
   (void) printf("SCCADD: no file delimiter...\n");
   n = 1;
   }
(void) strncpy(cbuf,catrc,n);
cbuf[n] = '\0';

if (strcmp(cbuf,output) == 0)
   {
   saverec = CATAL[cimno].RECNO - 1;            /* save record no. */
   if ((stat = rewicat(catid,cimno)) < 0) goto osa_error;

   while (saverec != CATAL[cimno].RECNO)
      {
      stat = readcat(catid,cimno,catrc,&valid);
      if (stat < 0) goto osa_error;
      }
   stat = (int) osaseek(catid,0L,FILE_CURRENT);		/* set write pointer */
   if (stat < 0) goto osa_error;

   CATAL[cimno].RECNO = saverec;

   if (catlen < catoff)			/* new length is larger */
      {
      catrc[0] = '!';
      catrc[1] = ' ';
      stat = osawrite(catid,catrc,catlen);		/* delete this entry */
      stat = (int) osaseek(catid,0L,FILE_END);     /* position write pointer */
      if (stat < 0) goto osa_error;		/* to EOF */

      (void) 
      sprintf(output,"updated entry #%-4.4d moved to end of catalog",saverec);
      SCTPUT(output);
      CATAL[cimno].RECNO = 9999999;		/* force rewind next time */
      }

   else if (catlen > catoff)			/* new length is smaller */
      {
      memset((void *)catrc,32,(size_t)CATREC_LEN);
      (void) strncpy(catrc,catrec,catoff);	/* new record -> old record */
      catrc[catlen] = '\0';
      (void) strcpy(catrec,catrc);
      }
   }
else
   goto find_loop;

	
/*  write new record	*/

not_found:
stat = osawrite(catid,catrec,catoff);
if (stat < catoff)
   {
   puts(osmsg());
   stat = ERR_INPINV;
   SCTPUT("(ERR) SCCADD: - INPINV");
   return stat;
   }
else
   {
   CATAL[cimno].RECNO ++;
   return ERR_NORMAL;
   }
 
osa_error:
stat = ERR_CATBAD;
SCTPUT("(ERR) SCCADD: - CATBAD");
return stat;
}

