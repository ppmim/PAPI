/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		test_function.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS: #include<midas_def.h>
.KEYWORDS
.PURPOSE		test function modules
.IN/OUTPUT		varies
.RETURNS
.COMMENTS
.VERSION
--------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Parameter and Symbolic Constant Files */
#include "pipe_para.h"		/* pipeline related parameters/constants */

/* detector related constants: OPRIME_para.h for OMEGA PRIME */
                                   /* O2000_para.h  for OMEGA 2000  */
/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif




/* Header Files of Secundary Modules */
#include "flat_correction.h"
#include "subframe_statistics.h"
#include "normalize_frame.h"
#include "extract_name.h"
#include "load_stack.h"
#include "copy_frame.h"


/* variable declerations */
char inframe1[80], inframe2[80], outframe[80];		/* names of input and output images */
char ident_1[80], ident_2[80];				/* ascii identifier of image */
char ident_out[80] = "test image";
char cunit[40];						/* unit of each axis */

float *pntr_1, *pntr_2, *pntr_out;			/* pointers to mapped data of input */
							/* images 1 and 2 and output frame  */
double start[2], step[2];				 /* start coordinate, step : ARRAYS! */ 

int stat;  	 					/* holds current return status */
int no_of_char;						/* number of characters of read in names */
int image_no1, image_no2;				/* image numbers of inframe 1 / 2 */
int image_no_out;					/* iamge number of output frame */
int naxis, npix[2];					/* # dimens of image, size of each dim */
int loop;						/* loop variable */
int unit = 0;
int objectnumber = 20;



void main()

{

  SCSPRO("test_function.c");


  /* set up program for error handling*/  
  //  SCECNT("PUT",&cont,&log_flag,&disp);

  /* copy image names from PRG paramters to character arrays */
  stat = SCKGETC("IN_A",1,75,&no_of_char,inframe1);
  //  stat = SCKGETC("IN_B",1,75,&no_of_char,inframe2);
  //  stat = SCKGETC("OUT_A",1,75,&no_of_char,outframe);
 

  //  stat = SCTPUT(inframe2);
  //  stat = SCTPUT(outframe);


  /* open the two input image files from disk, using high level interface */
  /*    stat = SCIGET(inframe1,D_R4_FORMAT,F_I_MODE,F_FIMA_TYPE,2,&naxis,npix,start,step,id ent_1,cunit,(char *)&pntr_1,&image_no1); */
  			        /* read_mode  */

  /* open an .FITS input frame in read-write-mode */
/*   stat = SCIGET(inframe1,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,ident_1,cunit,(char *)&pntr_1,&image_no1); */


  printf("Status SCIGET 1 = %d\n", stat);  

  if (stat != 0)
    SCETER(stat,"Grund fuer Abbruch 1");


/*   stat = SCIGET(inframe2,D_R4_FORMAT,F_I_MODE,F_FIMA_TYPE,2,&naxis,npix,start,step,ident_2,cunit,(char *)&pntr_2,&image_no2); */

/*   printf("Status SCIGET 2 = %d\n", stat);   */
/*   printf("Image_no 2 = %d\n", image_no2);  */

/*   if (stat != 0) */
/*     SCETER(stat,"Grund fuer Abbruch 2"); */


  /* create a frame for the output, using high level Interface */
  /*  stat = SCIPUT(outframe,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,ident_out,cunit,(char *)&pntr_out,&image_no_out); 	*/
  //  printf("Status SCIPUT = %d\n", stat);     /* .fits Format */

  //  if (stat != 0)
  //    SCETER(stat,"Grund fuer Abbruch 3");


  /* change continuation flag*/
  //  cont = 0;
  //  SCECNT("PUT",&cont,&log_flag,&disp);



  /* call test function */
  printf("Status vor Test = %d\n", stat); 
  
  stat = 100;				/* set status to 100 to check the function */

  //  stat = Flat(pntr_1, pntr_2, pntr_out, TOT_PIX, image_no_out, inframe1 );
  //  stat = Statistics(pntr_1, image_no1);
  //  stat = Normalize(pntr_1,image_no1,TOT_PIX);
  //  stat = Name(inframe1, outframe);


  //  stat = Copy(pntr_1,pntr_out, TOT_PIX, image_no_out, inframe1);





  /* Test the findobj routine */
  /* write necessary keywords */
  stat = SCKWRC("in_ima",1,inframe1,1,60,&unit);
  stat = SCKWRI("n_obj",&objectnumber,1,1,&unit);	/* find 20 elements */

  stat = SCTPUT(inframe1);
  printf("\nJetzt gehts zu findobj\n");

  system("/disk-v/fassbend/pipeline/OP_data/findobj.exe");









  printf("Status nach Test= %d\n\n", stat); 

  /* unmap all frames */
  //  stat = SCFUNM(image_no1);   
  //  stat = SCFUNM(image_no2); 

  //  printf("Status vor UNMAP des output frames: %d\n", stat);
  //  stat = SCFUNM(image_no_out); 
  //  printf("Status nach UNMAP des output frames: %d\n", stat);

  /* write descriptors of new frame: all but standard descriptors */
  //  stat = SCDCOP(image_no1,image_no_out,3,inframe1);


  /* close all frames */
  //  stat = SCFCLO(image_no1);
  //  stat = SCFCLO(image_no2);
  //  stat = SCFCLO(image_no_out);



  stat = SCTPUT("Meine erste Testfunktion ist erfolgreich verlaufen!");

  /* close down everything */


  SCSEPI();
 
}
