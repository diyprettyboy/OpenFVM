/* file calcu.c */
/* (c) Donald Axel GPL - license */
/* ANSI - C program demonstration, command line calculator */
/* Recursive descent parser */
/* Improve: Make a HELP command. Add more variables.       */

// modified original file calcu.c by OpenFVM team in 12/12/2005

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parser.h"

/* Description: */
/* Test the application with a given string.            */

/*
int main(int argc, char *argv[])
{
    
    int i, rv, jj;
    
    x = 2.0;
    y = 1.0;
    z = 0.0;

    jj = 0;
    while (++jj < argc) {
        strcat(gs, argv[jj]);
    }
    
    if (argc == 1)
        return calcu();
    
    strcat(gs, "\n");
    
    for (i = 0; i < 21; i++)
    {
    	
    	y = (double) i / 20.0 * 0.4;
	  
    	strcat(gs, "6*1.0/0.4^2*(2*y*0.2-y^2)\n");
    
    	rv = evaluate(gs, &oldval);
  
    	if (!rv)
        	printf("%.3f\n", oldval);
    	else
        	printf("Calc:%s\n%*s\n", gs, rv, "^Error");
    }
    
    return rv;
    
}*/

/* Description: */
/* calcu() sets up a string which is then evaluated as an expression  */
/* If (argc>1) main sets up string for evaluate() and prints result.  */
/* stricmp does not stop at '\n' - so we have to compare with "xx\n"  */
/* gettok() could solve that problem. TRY to use gettok().            */

int nextchar()
{
    ++cp;
    while (*cp == ' ')
        ++cp;
    return *cp;
}


int eatspace()
{
    while (*cp == ' ')
        ++cp;
    return *cp;
}

int calcu()
{
    FILE *ifil;
    char line[MAXL];
    int rpos;
    double r;

    ifil = stdin;
    while (1) {
        errorp = NULL;
        printf("Calc:");
        
	if (!fgets(line, MAXL, ifil))
            break;
        
	#ifdef WIN32
	if (strlen(line) && strnicmp(line,"QUIT",4) && stricmp(line,"Q\n"))
            rpos = evaluate(line, &r);
	else
            break;
	#endif
	
	#ifndef WIN32
        if (strlen(line) && strncasecmp(line,"QUIT",4) && strcasecmp(line,"Q\n"))
            rpos = evaluate(line, &r);
	else
            break;
        #endif
	
	if (!rpos) {
            printf("%-18g\n", r);
            oldval = r;
        } else {                /* prints Error in field min. 12 wide */
            printf("%*s\n", rpos, "^Error");
        }
    }
    return rpos;                /* if interactive rpos should always be 0 */
}

int evaluate(char *s, double *r)
{
    cp = s;
    eatspace();
    *r = expression();
    eatspace();
    if (*cp == '\n' && !errorp)
        return (0);
    else
        return (cp - s) + 11;
}


double expression()
{
    double e;
    int opera2;

    /* printf("test arg:%s\n",cp); */

    e = product();
    while ((opera2 = *cp) == '+' || opera2 == '-') {
        nextchar();
        if (opera2 == '+')
            e += product();
        else
            e -= product();
    }
    eatspace();
    return e;
}


double product()
{
    double dp;
    int ope;

    dp = potens();
    while ((ope = *cp) == '*' || ope == '/') {
        nextchar();
        if (ope == '*')
            dp *= potens();
        else
            dp /= potens();
    }
    eatspace();
    return dp;
}


double potens()
{
    double dpo;

    dpo = signedfactor();
    while (*cp == '^') {
        nextchar();
        dpo = exp(log(dpo) * signedfactor());
    }
    eatspace();
    return dpo;
}


double signedfactor()
{
    double ds;
    if (*cp == '-') {
        nextchar();
        ds = -factor();
    } else
        ds = factor();
    eatspace();
    return ds;
}


double factor()
{
    double df;

    /* while (*cp!='\n') {
       putchar(*cp++);
       } 
     */

    switch (*cp) {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
        df = strtod(cp, &cp);
        break;
    case '(':
        nextchar();
        df = expression();
        if (*cp == ')')
            nextchar();
        else
            errorp = cp;
        break;
    case 'M':
        nextchar();
        df = oldval;
        break;

    case 'x':
        nextchar();
        df = x;
        break;

    case 'y':
        nextchar();
        df = y;
        break;

    case 'z':
        nextchar();
        df = z;
        break;
			
    default:
        df = stdfunc();
    }
    /* printf("ddt: df = %lf, *cp = %c\n",df,*cp); */

    eatspace();
    return df;
}


char *functionname[] =
{
    "abs", "sqrt", "sin", "cos", "atan", "log", "exp", "\0"
};

double stdfunc()
{
    double dsf;
    char **fnptr;
    int jj;

    eatspace();
    jj = 0;
    fnptr = functionname;
    while (**fnptr) {
        /* printf("%s\n",*fnptr); */
        if (strncmp(*fnptr, cp, strlen(*fnptr)) == 0)
            break;

        ++fnptr;
        ++jj;
    }
    if (!**fnptr) {
        errorp = cp;
        return 1;
    }
    cp += (strlen(*fnptr) - 1);
    nextchar();
    dsf = factor();
    switch (jj) {
    case 0: dsf = fabs(dsf);  break;
    case 1: dsf = sqrt(dsf); break;
    case 2: dsf = sin(dsf);  break;
    case 3: dsf = cos(dsf);  break;
    case 4: dsf = atan(dsf); break;
    case 5: dsf = log(dsf);  break;
    case 6: dsf = exp(dsf);  break;
    default:{
            errorp = cp;
            return 4;
        }
    }
    eatspace();
    return dsf;
}


/* end calcu.c */

