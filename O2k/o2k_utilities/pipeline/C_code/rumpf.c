/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS: #include<midas_def.h>
.KEYWORDS
.PURPOSE	       
.IN/OUTPUT		
.RETURNS
.COMMENTS
.VERSION
--------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Parameter and Symbolic Constant Files */
#include "pipe_para.h"		/* pipeline related parameters/constants */
#include "OPRIME_para.h"	/* detector related constants: OPRIME_para.h for OMEGA PRIME */
								   /* O2000_para.h  for OMEGA 2000  */

/* Header Files of Secundary Modules */



main()

{

  SCSPRO("");



  SCSEPI();
 
}
