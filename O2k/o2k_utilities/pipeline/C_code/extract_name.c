/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		extract_name.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 4 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		image name
.PURPOSE	     	extract part of input image name to create a related name for the
			output image 
.IN/OUTPUT		IN: pointer to name of input image
			IN: pointer to name of output image
.RETURNS		Status: 0 = ok
.COMMENTS		Depending on the DETECTOR used, output image names are created from
			the input image names.
			The output names will have the following form:
			For OPRIME: red_???????????.fits , where ? are the last characters
			of the input image name. 
			For OMEGA2000: the original names will be taken with the prefixes
					  red_name or sky_name for reduced and sky images
.CALL			int Name(char *p_in_name, char *p_out_name)
.VERSION		1.00   		02.09.02
			1.1			23.10.02	name for .fits files
			2.0			22.12.02	names for OMEGA2000
			2.1			06.03.03	changed all output to SCTPUT
			3.0			08.04.03	added absolute path handling
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<string.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "extract_name.h"

/* Pipeline Parameters, including DETECTOR */
#include "pipe_para.h"



int Name(char *p_in_name, char *p_out_name)
{
  /* Local Variables */
  int loop = 0;							/* loop variable */

  char *p_slash_pos;						/* pointer to slash-position */




  /* check which DETECTOR is used */

  if(DETECTOR == 1)						/* 1 = OMEGA_PRIME */
    {
      /* first 4 characters of out_name have been initialized @ declaration with: red_ */

      /* set pointers to start positions */
      p_out_name = p_out_name + 4 ;				/* start with 5th character */
      p_in_name = p_in_name + 20;				/* start with 20th character */


      for(loop=0 ; loop<=11 ; loop++)
	{
	  *p_out_name = *p_in_name;		  		/* take part of inframe name */

	  p_out_name++; 
	  p_in_name++;	
	}

      *p_out_name = 'f';
      p_out_name++;  

      *p_out_name = 'i';
      p_out_name++;
  
      *p_out_name = 't';
      p_out_name++;

      *p_out_name = 's';
      p_out_name++;

      *p_out_name = '\0';				   	/* result is: red_???????????.fits */

      return 0;
    }



  else if(DETECTOR == 0)					/* 0 = OMEGA2000 */
    {
 
      /* check in-name for last slash : / = ASCII character #47 */
      /* returns pointer to slash or NULL if not found */
      p_slash_pos = strrchr(p_in_name,47);


      /* extract name according to last slash finding */
      /* no slash found */
      if(p_slash_pos == NULL)
	{
	  /* add full name to the prefixes red_ or sky_ with strncpy()*/
	  /* start with the 5th position */
	  strncpy(p_out_name+4,p_in_name,80);
	}

      /* slash found */
      else
	{
	  /* start with first position after last slash */
	  strncpy(p_out_name+4,p_slash_pos+1,80);
	}
      


      return 0;
    }







  else
    {
      SCTPUT("Error: Wrong Detector Parameter");

      return 7;							/* 7 = input invalid */
    }

}
