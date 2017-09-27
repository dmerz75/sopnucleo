
                     /* this is the file identify.c */ 

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


void identify(d1,d2,d3,x,y,z)
double d1,d2,d3,*x,*y,*z;
{

	  *x = d1;
	  *y = d2;
	  *z = d3;

}

void dist(x1,x2,y1,y2,z1,z2,d)
double x1,x2,y1,y2,z1,z2, *d;
/* this fuction calculates the distance between 2 points in 3D space */

{

  *d = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}



int string_to_integer(char string[])
/* this function converts a string into an integer */

{
  int i, integer_value, result = 0;

  i = 0;
  while(string[i] != '\0'){

    if(string[i] >= '0' && string[i] <= '9'){
      integer_value = string[i] - '0';
      result = result*10 + integer_value;
    }

    i++;

  }

  return(result);
  
}


/* determine the number (between 1 and 20) corresponding to the type of the 
   amino acid when the 1-letter code is used */

int identity_aa_one(char aa){

  int typ;

  if(     (aa == 'C') || (aa == 'c')) typ = 1;
  else if((aa == 'F') || (aa == 'f')) typ = 2;
  else if((aa == 'L') || (aa == 'l')) typ = 3;
  else if((aa == 'W') || (aa == 'w')) typ = 4;
  else if((aa == 'V') || (aa == 'v')) typ = 5;
  else if((aa == 'I') || (aa == 'i')) typ = 6;
  else if((aa == 'M') || (aa == 'm')) typ = 7;
  else if((aa == 'H') || (aa == 'h')) typ = 8;
  else if((aa == 'Y') || (aa == 'y')) typ = 9;
  else if((aa == 'A') || (aa == 'a')) typ = 10;
  else if((aa == 'G') || (aa == 'g')) typ = 11;
  else if((aa == 'P') || (aa == 'p')) typ = 12;
  else if((aa == 'N') || (aa == 'n')) typ = 13;
  else if((aa == 'T') || (aa == 't')) typ = 14;
  else if((aa == 'S') || (aa == 's')) typ = 15;
  else if((aa == 'R') || (aa == 'r')) typ = 16;
  else if((aa == 'Q') || (aa == 'q')) typ = 17;
  else if((aa == 'D') || (aa == 'd')) typ = 18;
  else if((aa == 'K') || (aa == 'k')) typ = 19;
  else if((aa == 'E') || (aa == 'e')) typ = 20;
  else typ = 0;

  return(typ);

}

