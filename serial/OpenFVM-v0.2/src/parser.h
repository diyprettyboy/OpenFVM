
/* (c) Donald Axel GPL - license */
/* ANSI - C program demonstration, command line calculator */
/* Recursive descent parser */
/* Improve: Make a HELP command. Add more variables.       */

// modified original file calcu.c by OpenFVM team in 12/12/2005

#define MAXL 8196

char gs[MAXL];
char *cp;
char *errorp;
double oldval;
double x, y, z;

/* local prototypes: */
int calcu();
int evaluate(char *line, double *prev_result);

/* More local prototypes. This could, of course, be a separate file. */
double expression();
double product();
double potens();
double signedfactor();
double factor();
double stdfunc();
