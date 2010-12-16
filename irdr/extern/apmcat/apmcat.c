/*
 * Get data from the APM catalogues server. 
 *
 * This is a standalone program which makes a request of the
 * APM catalogues system based on command line arguments and returns
 * either a list or a postscript finding chart on standard output.
 *
 * VMS Notes:
 *    This program assumes Multinet TCP/IP support.
 *    This program expects to run as a foreign command.
 *    If standard output (SYS$OUTPUT) is redirected to a file (using
 *      ASSIGN or DEFINE) the resulting file will have control
 *      characters before each record.  
 *   
 * This version should compile under Ultrix and VMS/Multinet.
 *
 * Created: 1-February-1995 by T. McGlynn
 *          Goddard Space Flight Center
 *          Universities Space Research Association
 *          Code 668.1
 * Modified by Geraint Lewis and Mike Irwin to support APM online catalogues
 *          1-April-1996
 *
 */
/*#define APMCAT "131.111.68.247" - alternative form if name resolver is crap*/
#define APMCAT "www.ast.cam.ac.uk"
#define PORT	80
#define ERROR  -1

#ifdef VMS  /* Assumes Multinet support */
#include "multinet_root:[multinet.include.sys]types.h"
#include "multinet_root:[multinet.include.sys]socket.h"
#include "multinet_root:[multinet.include.netinet]in.h"
#include <stdio.h>
#include "multinet_root:[multinet.include]netdb.h"

   #include "multinet_root:[multinet.include]netdb.h"

/* To use with the UCX TCPIP stack, replace the directories in 
multinet_root: in ALL FOUR CASES above with 

   DISK$SYSTEM:[VMS$COMMON.DECC$LIB.REFERENCE.DECC$RTLDEF]

where DISK$SYSTEM is the (logical) name of the disk where the files 
reside and compile with CC/STANDARD=VAXC.  (Tested with UCX V4.1-12G and
DEC C V5.5-002 on OpenVMS V7.1 Alpha) - Phillip Helbig */

#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stdio.h>
#include <netdb.h>

#endif

/* Define return types of functions */
char	*mo_code();
char	*form_request();
char	*strchr();
void	lowcase();


/* Request buffer */
char	req_str[2048];

char	*out_file = 0;

/* Is this an informational request */
int	info_req = 0;


int main(argc, argv)

int	argc;
char	*argv[];
{
    int		name_given=0;		/* Name or lists? */

    /* Request fields */
    char	survey[64], equinox[64], email[64];
    char	box[64], img[64], list[64];
    char	coords[1024], numbers[64];

    int		i;
    char	*p;
    
    /* Initialize all fields to zero */
    survey[0]  = 0;
    equinox[0] = 0;
    email[0] = 0;
    box[0] = 0;
    img[0] = 0;
    list[0] = 0;
    coords[0] = 0;
    numbers[0] = 0;
    
    /* Check that there are enough arguments  */
    if (argc < 2) 
    {
	printf("Usage: apmcat ra dec [optional arguments]\n");
	printf("\n");
	printf("   where the optional arguments are of the form\n");
	printf("   keyword=value.  Valid arguments include:\n");
	printf("              survey=poss1      [ukst]\n");
	printf("              box=5             [box size in arcmins]\n");
        printf("              equinox=b1950     [j2000]\n");
	printf("	      ps=image.ps       [name of postscript image]\n");
	printf("	      list=image.lis    [name of list file]\n");
	printf("              numbers=n         [numbers on plot ? n/y]\n");
	printf("              email=null        [email address]\n");
	printf("    The above values are the defaults.\n");
	printf("    The default is to return a postscript \n");
	printf("    image with the name image.ps \n\n");
	printf("    Please specify either a list or a postscript image\n");
	printf("\n");
	printf("Examples:\n");
	printf("    apmcat \"00 40 00.00\" \"11 41 00\" survey=poss1 ps=chart.ps\n");
	printf("    apmcat \"10 31 21.3 -5 10 23\" survey=ukst list=chart.list\n");
	printf("\n");
	exit();
    } 

    /* Convert all arguments to lower case */
    for (i=0; i<argc; i++)
    	lowcase(argv[i]);
    
    /* Check if we have an RA and DEC fields or a single field. */
    
    if (argc == 2)
    	name_given=1;
    else if (strchr(argv[2], '='))
        name_given=1;

    /* Default output file */
    
    out_file = "image.ps" ;

    if (name_given) {
	i = 2;
	sprintf(coords, "RADEC=%s", argv[1]);
    }
    else
    {
        i=3;
	sprintf(coords, "RADEC=%s %s", argv[1], argv[2]);
    }

    /* Loop over other arguments */	     
    for (;i<argc; i++) 
    {
	
        p = argv[i];    
	
	if (!strncmp(p, "email", 5) ) 
        {
	    strcpy(email,"EMAIL");
	    strcat(email,p+5);
	    
	}

	else if (!strncmp(p, "survey=",6) )
	{
	    strcpy(survey, "CAT=");
	    strcat(survey, p+7);
	}
	
	else if (!strncmp(p, "equinox",7) )
	{
	    strcpy(equinox, "EQUINOX=");
	    strcat(equinox, p+8);
	}
	else if (!strncmp(p, "box",3) )
	{
	    strcpy(box, "BOX=");
	    strcat(box, p+4);
	}
	else if (!strncmp(p, "numbers",7) )
	{
	    strcpy(numbers, "NUMBERS=");
	    strcat(numbers, p+8);
	}
	else if (!strncmp(p, "ps=",3) )
	{
	    strcpy(img, "IMG=");
	    strcat(img, p+3);
	    out_file = p+3 ;
	}
	else if (!strncmp(p, "list=",5) )
        {
	    strcpy(list, "LIST=");
	    strcat(list, p+5);
	    out_file = p+5 ;
	}

    }

    p= form_request(coords, survey, equinox, list, box, numbers,
		    img, email);
    
    strcpy(req_str, p);
    
    /*    fprintf(stderr, "%s\n", req_str); */
    getimage();

    
    return 0;
}

char *form_request(coords, survey, equinox, list, box, numbers,
		   img, email)

char	*coords, *survey, *equinox, *list, *box, *numbers;
char	*img, *email;

{
    static char buf[2048];
    char	lbuf[128];
    char	*t = buf;
    
    
    /* Concatenate argument fields together.  Add in defaults
     * for arguments not specified.
     */
    strcpy(buf, mo_code(coords));
    
    if ( *list && *img ) {
      printf("\n\n Please specify either an image OR a list\n\n");
      exit();
    } 
     
    strcat(buf, "&");
    if (*survey)
    {
	
	strcat(buf, mo_code(survey));
    } else 
    {
	strcat(buf, "SURVEY=poss1");
    }

    strcat(buf, "&");
    if (*list)
    {
	strcat(buf, mo_code(list));
    } else 
    {
	strcat(buf, "LIST=off");
    }
    
    strcat(buf, "&");
    if (*equinox)
    {
	strcat(buf, mo_code(equinox));
    } else 
    {
	strcat(buf, "EQUINOX=b1950");
    }
	
    strcat(buf, "&");
    if (*box)
    {
	strcat(buf, mo_code(box));
    } else 
    {
	strcat(buf, "BOX=5");
    }

    strcat(buf, "&");
    if (*numbers)
    {
	strcat(buf, mo_code(numbers));
    } else 
    {
	strcat(buf, "NUMBERS=n");
    }

    strcat(buf, "&");
    if (*img)
    {
	strcat(buf, mo_code(img));
    } else 
    {
	strcat(buf, "PS=file.ps");
    }

    strcat(buf, "&");
    if (*email)
    {
	strcat(buf, mo_code(email));
    } else 
    {
	strcat(buf, "EMAIL=null");
    }
    
    strcat(buf, "\013\012\012");

    return buf;
}


/* Send request and receive data.
 */
int getimage()

{
	int s, n;
	char buf[16384];
        char nbuf[16384];
   	char nbuf2[16384];

	struct sockaddr_in sin;
	struct hostent *hp;
	struct servent *sp;
        int val;
        int i;

	hp = gethostbyname(APMCAT);
	if (hp == NULL) {
		return ERROR;
	}
	/*
	 *  Create an IP-family socket on which to make the connection
	 */

	s = socket(hp->h_addrtype, SOCK_STREAM, 0);
	if (s < 0) {
		return ERROR;
	}


	sin.sin_family = hp->h_addrtype;
	memcpy(&sin.sin_addr, hp->h_addr, hp->h_length);
	sin.sin_port = htons(PORT);

	/*
	 *  Connect to that address...
	 */

	if ((val=connect(s, &sin, sizeof(sin))) < 0) {
		return ERROR;
	}

	/* Send first part of request */
	strcpy(buf, "POST /apmcatbin/post-query HTTP/1.0\012");
   
        if (send(s, buf, strlen(buf),0) < 0)
    		return ERROR;

        strcpy(buf, "User_Agent: SkyView Image Selector\012");
        strcat(buf, "Content-type: application/x-www-form-urlencoded\012");
        sprintf(nbuf2, "Content-length: %d\012\012", strlen(req_str));

        strcat(buf, nbuf2);
        strcat(buf, req_str);

        /* Send request contents */
        if (send(s, buf, strlen(buf), 0) < 0)
   		return ERROR;
         /* printf("%s",buf); */

	/* Now get back the data */
        if (ptrans(s) < 0)
    		return ERROR;
}

char	*mo_code(string)

char	*string;

{
        /* Perform Mosaic encoding.  Only alphanumerics are
	 * unchanged.  Spaces are replaced by +.  All others by
	 * %xx where xx is the Hex code for the value.
	 * Note this routine uses a static pointer so the value
	 * should be copied before this routine is called again.
	 */
	static char buf[2048];
    	char	*s,*t;
        char	c;
    
    	s = string;
        t = buf;
        do {
	
	    c = *s;
	    *t++ = *s++;
	}
    
    	while (c != '=');
    
    	c = *s;
    	while (c) 
    	{
		if (c == ' ') *t++ = '+';
	    
	        /* Assume ASCII sequencing */
	        else if ( !  ((c >= 'a' && c <= 'z') ||
			 (c >= 'A' && c <= 'Z') ||
			 (c >= '0' && c <= '9')
			 )) {
			sprintf(t, "%%%+2x", c);
		        t += 3;
		} else 
	        {
		    *t = c;
		    t++;
		}
	    
		
	        c = *++s;
	}
    
        *t = 0;
	return buf;
}

/* Read from input until we can fill the buffer */
int getbuf(s, rbuf, qlen)

int s;
char	*rbuf;
int qlen;
{
    

	static int left=0;
        static char bigbuf[32768];
        int	n;
	int	len = qlen;
        int     i;

        while (len > left) {
	    n = recv(s, bigbuf+left, sizeof(bigbuf)-left, 0);
	    if (n <= 0) {
		if (left <= 0)
		    	return( -1);
		else {
			len = left;
			break;
		}
	    }
	    else 
	    {
		left = left + n;
	    }
	}
        memcpy(rbuf, bigbuf, len);
        left = left - len;
        if (left > 0)
	  for (i = 0; i < left; i++)
	    bigbuf[i] = bigbuf[i+len];
/*	        memmove(bigbuf, bigbuf+len, left);     */
    	return len;
}

/* Transfer data from socket to standard output
 */
int ptrans(s)

int	s;

{
	char	line[16384];
    	int 	lcnt = 0;
    	int	scaling = 0;
        char	lbuf[20];
	char	bigbuf[2880];
    	int	nlset=0;
    	char	c;
    	int	fp=1;

        /* Look for a double newline to signal the beginning of the data */
        while (1) {
	    if (getbuf(s, &c, 1) < 0) 
		return ERROR;
	    if (c == '\012') 
	    {
		if (nlset)
			break;
		else
			nlset = 1;
	    } else nlset = 0;
	}
    
        if (out_file) if (*out_file) 
        {
	    if (info_req)
		fp = creat(out_file, 0644);
	    else
            {
	        fp = creat(out_file,
#ifdef VMS  /* Use default file mode for VMS */
		       0,"bls=2880","rfm=fix", "mrs=2880"
#else
		       0644
#endif
		       );
	    }
	    if (fp <=0) 
	    {
		printf("ERROR: Unable to open output file %s\n", out_file);
		perror("  creat");
	    }
	}

	lcnt=getbuf(s, line, 1);
        while(1) {
	    lcnt=getbuf(s, line, 2880);
	    if (lcnt <= 0) break;
	    write(fp, line, lcnt);
	}
	return 0;
}

void	lowcase(s)

char	*s;

{
    
    /* This procedure converts a string, s, to lower case */
    
    char	*p;
    p = s;
    while (*p) 
    {
	if (isalpha(*p)) *p = tolower(*p);
	p++;
    }
}
