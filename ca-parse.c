/* file: ca-parse.c							*/

#define INT_EXT	extern 	/* causes all variables in ca.h to be external	*/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "ca.h"
#include <values.h>

void run_evolve();
void iterate_ca(struct CA *source_ca);
int disp_plot();
int disp_window();


/************************************************************************
 *									*
 * check_for_string() ensures that the next input parsed is a string	*
 * or a string variable.  If it is, a pointer to the string (encoded	*
 * as a double) is returned.  If the next input is not a string -1 is	*
 * returned.  *next_char is updated to point past the input or at the	*
 * end of the input string if -1 is returned.				*
 *									*
 ************************************************************************/
double check_for_string(s, next_char)
char *s, **next_char;
{
    double rv;

    if ((rv = set_option(s, next_char)) > 0)	/* make sure positive	*
						 * (pointer) is returned*/
        if (return_type == STRING)	/* make sure it is a string	*/
	    return(rv);					/* return it	*/

    if (return_type == ERROR) /* error has been displayed already	*/
	return(-1);

				/* either NULL returned or not a string	*/
    *next_char = s + strlen(s);
    parse_err = 1;
    printf("parse error: string argument expected=> %s \n", s);
    return_type = ERROR;
    return(-1);
} /* end check_for_string() */


/************************************************************************
 *									*
 * parse_ca() sets values associated with cellular automata returns	*
 * -1 if value passed is out of range or command is not valid		*
 *									*
 ************************************************************************/
double parse_ca(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;
    double lambda, max_rule, max_double;
    char *rule_string,c;
    int i= 0, j = 0, ts;
    long temp;


    switch(*(s + 1)) {				/* command letter	*/

      case 'n':					/* neighborhood size	*/
	rv = set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, or oops	*/
	    if (return_type != ERROR) {
	   	string_error(s);
	    }
	    break;					/* return -1	*/
	}
	if (rv < 0) {
	    root.neighborhood = 1;
	    parse_error(s);
	    break;
	}
	else if (rv > 0) {
	    if (pow((double) root.states, (double)(rv*2 + 1)) > MAX_RULES) {
		printf(
"neighborhood exceeds rule table size, reduce number of states first\n");
		parse_error(s);
		*next_char = s + strlen(s);
		return(-1);
	    }
	    root.neighborhood = (int) rv;
	}
	else			/* rv == 0, return current value	*/
	    rv = root.neighborhood;

	if (verbose)
	    printf("root CA neighborhood size = %d\n", root.neighborhood);
	break;


      case 's':					/* number of states	*/
	rv = set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		    string_error(s);
	    break;					/* return -1	*/
	}
	if (rv < 0) {
	    root.states = 2;
	    parse_error(s);
	}
	else if (rv > 0) {
	    if (pow(rv, (double) (root.neighborhood*2 + 1)) > MAX_RULES) {
		printf(
"number of states exceeds rule table size, reduce neighborhood first\n");
		parse_error(s);
		*next_char = s + strlen(s);
		return(-1);
	    }
	    root.states = (int) rv;
	}
	else				/* rv==0, return current value	*/
	    rv = root.states;

	if (verbose)
	    printf("root CA number of states = %d\n", root.states);
	break;


      case 'l': 	/* lambda parameter (% of rules that are not 0)	*/
	lambda = set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		    string_error(s);
	    break;					/* return -1	*/
	}
	if ((lambda < 0) || (lambda > 1)) {	/* lambda out of range	*/
	    printf("Error: lambda value must be between 0 and 1\n");
	    parse_error(s);
	    rv = -1;
	    break;					/* return -1	*/
	}
	if (lambda != 0)			/* if lambda==0, don't	*
						 * change root.lambda	*/
	    root.lambda = lambda;

				/* calculate current lambda value	*/
	ts = (int) pow((double)root.states, (double)(2*root.neighborhood + 1));
	for (i=0, j=0; i<ts; i++)
	    if (root.rule_table[i])
		j++;

	rv = ((double) j) / ((double) ts);

	if (verbose)
	    printf(
	"lambda of current rule= %g. lambda of next random rule = %g\n", rv,
								root.lambda);
	break;


      case 'r': 		/* build rule table from integer (only	*
				 * works for small rule tables		*/
	ts = (int) pow((double) root.states, (double)(root.neighborhood*2 + 1));

	max_rule = pow((double) root.states, (double) ts);
	max_double = DBL_MAX;

	if (max_rule > max_double) {
	    printf("max rule = %g, unable to parse rule numbers beyond %g\n",
							max_rule, max_double);
	    rv = -1;
	    *next_char = s + strlen(s);
	    break;
	}

	rv = set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		    string_error(s);
	    break;
	}
	if (rv < 0) {
	    parse_error(s);
	    break;
	}
	if (get_space(&rule_string, ts+1, sizeof(char))) {
	    rv = -1;
	    break;
	}
	temp = (long) rv;
	for (i=0; i<ts; i++) {
	    root.rule_table[i] = temp % root.states;
	    temp /= root.states;
	    rule_string[ts - i - 1] = root.rule_table[i] + '0';
	}
	rule_string[ts] = '\0';
	if (verbose)
	    printf("rule table = \"%s\"\n", rule_string);
	return_type = STRING;
	rv = (double)(long long) rule_string;
	break;


      case 't':				/* max time of iterations	*/
	rv = (int) set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		    string_error(s);
	    break;					/* return -1	*/
	}
	if ((rv < 0) || (rv > MAX_TIME)) {	/* out of bounds	*/
	    printf("Error iterations must be between 0 and %d\n", MAX_TIME);
	    parse_error(s);
	    break;
	}
	else if (rv > 0)
	    max_time = (int) rv;
	else				/* rv==0, return current value	*/
	    rv = (double) max_time;

	if (verbose)
	    printf("maximum iterations of CA = %d\n", max_time);
	break;


      case 'w':						/* world size	*/
	rv = set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		    string_error(s);
	    break;					/* return -1	*/
	}
	if ((rv < 0) || (rv > MAX_WORLD)) {	/* out of bounds	*/
	    printf("world size out of range, max = %d\n", MAX_WORLD);
	    parse_error(s);
	    rv = -1;
	    break;
	}
	else if (rv > 0) {
	    world_size = (int) rv;
				/* generate enough space for running 	*
				 * complete iterations			*/
	    free(ca_states);
	    free(root.world);
	    if (get_space(&ca_states, MAX_TIME*world_size, sizeof(char)) ||
	        get_space(&root.world, world_size, sizeof(char))) {
	        break;					/* return -1	*/
	    }
	}
	else				/* rv==0, return current value	*/
	    rv = (double) world_size;

	if (verbose)
	    printf("CA world size = %d\n", world_size);
	break;

      default:							/* oops	*/
	parse_error(s);
	*next_char = s+strlen(s);
	break;						/* return -1	*/
    } /* end switch *(s+1) */

    return (rv);
} /* end parse_ca() */


/************************************************************************
 *									*
 * parse_display() displays graphic representations of CA activity	*
 *									*
 ************************************************************************/
double parse_display(s, next_char)
char *s, **next_char;
{
    int i,j;
    char buf[10];
    double rv = 0;

    return_type = NRV;				/* no return value	*/
    switch (*(s + 1)) {
      case 'w':					/* display X window	*/
	if (getenv("DISPLAY") == NULL) {
	    printf("Not in X window environment\n");
	    *next_char = s+strlen(s);
	    break;
	}
				/* start a new process to display	*/
	if (fork() == 0) {				/* child	*/
#ifndef GPROF				/* X11 lib can't support gprof	*/
	    disp_window();		/* displays window on screen	*/
#endif
	    exit(0);				/* exit when finished	*/
	}
	window_count++;

        *next_char = s+2;/* continue parsing input after Dw command	*/
	break;


      case 'a':				/* display iterations fully	*/
      case 's':			/* display automata one line at a time	*/
	*next_char = s+2;	/* continue parsing after Da or Ds	*/

	if (*(s+1) == 's')	/* Ds, display one line at a time	*/
	  printf("Enter for new line, x to quit\n");
	else 		/* Da, display entire contents of ca_states	*/
	  printf("Type ^c to exit display\n");

	int_hit = 0;				/* clear ^C indication	*/

	for(i=0; i<time_step; i++) {		/* for each time step	*/
	  for(j=0; j<world_size; j++){		/* for each position	*/
	    if (ca_states[i*world_size + j] == 0)
	      printf(" ");		/* space for quiescent state	*/
	    else if (root.states == 2)	/* use X for two state rules	*/
	      printf("X");
	    else			/* use state # for > 2 states	*/
	      printf("%1d",ca_states[i*world_size + j]);


	  }
	  if (*(s+1) == 's') {	/* Ds, look for input after each line	*/
	    if (fgets(buf, sizeof(buf), stdin) == NULL) {
			break;
		}
	    if ((*buf == 'x') || (*buf == 'X'))
	      break;
	  }
	  else {		/* Da, look for ^C after each line	*/
	    printf("\n");
	    if (int_hit) {		/* ^C struck, clear indicator	*/
	      int_hit = 0;
	      break;					/* exit loop	*/
	    }
	  }
	}
	break;					/* end case 'a' and 's'	*/

	/* TEK added 4/18/91 -- UNIX(tm) plot support...
	 * specific terminal depends on compilation library.
	 */
	case 'p':	/* UNIX(tm) plot */
        *next_char = s+2;/* continue parsing input after Dw command	*/

	disp_plot();	/* Display on graphics terminal */

	break;

      default:							/* oops	*/
	  *next_char = s + strlen(s);
	  parse_error(s);
	  rv = -1;
	  break;
    } /* end switch *(s+1) */

    return(rv);
} /* end parse_display() */


/************************************************************************
 *									*
 * parse_edit() displays the user with the current state of either the	*
 * rules or the world and allows the user to enter new values for them	*
 *									*
 ************************************************************************/
double parse_edit(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;
    char *line, *input, *rule_string, c;
    int i, rule_next_char, ts;

    switch (*(s + 1)) {
      case 'w':					/* edit world state	*/

				/* make sure next input is a string	*/
	rv = check_for_string(s+2, next_char);
	if (rv <= 0)
	    return(rv);

	input = (char*) (long long) rv;


	for (i=0; i<world_size; i++) {	/* check input for errors	*/
	  if (input[i] == '\0')
	    break;

	  if ((input[i] < '0') || (input[i] >= root.states + '0')) {
	    parse_error(&input[i]);		/* input out of range	*/
	    free(input);
	    return(-1);
	  }
	}
				/* allocate space for return line	*/
	if (get_space(&line, world_size + 1, sizeof(char)))
	     return(-1);

	for (i=0; i<world_size; i++) {	/* copy input back into world	*/
	  if (input[i] == '\0')		/* input is short of world_size	*/
	    break;
	  line[i] = input[i];
	  root.world[i] = input[i] - '0';
	}

	free(input);				/* release space	*/

	for ( ; i<world_size; i++)	/* copy world state into 	*
					 * remaining space in line[]	*/
	    line[i] = root.world[i] + '0';
	line[i] = '\0';					/* terminate	*/

	rv = (double)(long long) line;
	return_type = STRING;		/* indicate string returned	*/

	if (verbose)
	    printf("world = \"%s\"\n", line);

	break;						/* end 'w'	*/


      /******************************************************************
       * 'r' rules, a string of rules is expected.  The string can take	*
       * one of many forms: "+" "-" "z" or "####rr###r...#r##" where	*
       * # is a number >= 0 and < root.states and r is the letter r.	*
       * If the rule is "+" or "-" then the next sequential rule is 	*
       * generated, note that all rules map the 0 state surrounded by	*
       * 0 states back to state 0, i.e. rule "01101110" becomes 	*
       * "01110000" under "+" and "01101100" under "-".  "z" causes	*
       * the next state to be the set of zero rules: "000...0000"	*
       * A string of numbers interspersed with r's causes the rule to 	*
       * become those states that are specified by numbers and random 	*
       * for those rules that have been specified by 'r'.		*
       * e.g. "0110r010" can generate the rule "01100010" or "01101010" *
       * A zero length input string ("") means return a string 		*
       * representing the current rules.				*
       ******************************************************************/
      case 'r': /* set new CA rules, return string pointing to rules	*/
	if ((rv = check_for_string(s+2, next_char)) < 0)
	    return(rv);					/* not a string	*/

	rule_string = (char*) (long long) rv;
					/* calculate number of rules	*/
	ts = (int) pow((double) root.states, (double)(root.neighborhood*2 + 1));

	if (strlen(rule_string) != 0) {	/* non-zero string length	*/
	  switch(rule_string[0]) {
	  case '+':				/* next sequential rule	*/
	      for (i=0; i<ts; i++) {		/* start at rule 1	*/
		  root.rule_table[i]++;
		  if (root.rule_table[i] == root.states)/* carry	*/
		      root.rule_table[i] = 0;
		  else					/* no carry	*/
		      break;
	      }
	      break;

	  case '-':				/* decrement rule table	*/
	      for (i=0; i<ts; i++) {
		  if (root.rule_table[i] == 0)		/* borrow	*/
		      root.rule_table[i] = root.states - 1;
		  else {				/* no borrow	*/
		      root.rule_table[i]--;
		      break;
		  }
	      }
	      break;

	  case 'z':					/* zero rule	*/
	      for (i=0; i<ts; i++) {
		  root.rule_table[i] = 0;
	      }
	      break;

	  default:      			/* interpret string	*/
	    for (i=0; (i<ts) && (c = rule_string[i]); i++) {
	      if (c=='r') {				/* random state */
		  root.rule_table[ts - i -1] +=
			1 + (int)((double) (root.states - 1) * my_random());
		  root.rule_table[ts - i -1] %= root.states;
	      }
	      else {			/* copy state from string	*/
	    	  root.rule_table[ts - i - 1] = c - '0';
		  root.rule_table[ts - i -1] %= root.states;
	      }
	    }
	  } /* end switch(rule_string[0]) */
	} /* end strlen(rv) != 0 */

	free(rule_string);
			/* make space for whole rule table string	*/
	if (get_space(&rule_string, ts + 1, sizeof(char)))
	    break;
	for (i=ts-1; i>=0; i--)		/* copy rule table to string	*/
	    rule_string[ts-1-i] = root.rule_table[i] + '0';

	rv = (double)(long long) rule_string;/* return pointer to string	*/
	return_type = STRING;	/* indicate string being returned	*/

	if (verbose)
	    printf("rule table = \"%s\"\n", rule_string);
	break;						/* end 'r'	*/


      default:						/* oops		*/
	*next_char = s + strlen(s);
	parse_error(s);
	break;
    } /* end switch */

    return(rv);
} /* end parse_edit() */


/************************************************************************
 *									*
 * parse_execute() runs iterations of CA or of GA populations		*
 *									*
 ************************************************************************/
double parse_execute(s, next_char)
char *s, **next_char;
{
    double rv;

    *next_char = s+2;
    rv = 0;
    return_type = NRV;				/* no return value	*/
    switch (*(s + 1)) {
	case 'c':		/* iterate CA for max_time steps	*/
	  iterate_ca(&root);
	  break;

	case 'p':			/* run until number of 		*
			 	 	 * generations specified are	*
					 * exhausted			*/
	  run_evolve();
	  break;

	default:					/* oops		*/
	  parse_error(s);
	  rv = -1;
    }

    return(rv);
} /* end parse_execute() */


/************************************************************************
 *									*
 * parse_ga() sets values associated with genetic algorithm		*
 *									*
 ************************************************************************/
double parse_ga(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;

    switch (*(s + 1)) {
	case 'm':				/* mutation rate	*/
	  mutation_rate = (float)set_option(s+2, next_char);
	  if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
	      string_error(s);
	    break;					/* return -1	*/
	  }
	  if ((mutation_rate<0) || (mutation_rate>1.0)) {
	      printf("mutation rate must be between 0 and 1\n");
	      parse_error(s);			/* out of range		*/
	      break;
	  }
	  if (verbose)
	      printf("GA mutation rate = %g\n", mutation_rate);
	  rv = mutation_rate;
	  break;

	case 'a':				/* augmentation rate	*/
	  augmentation_rate = (float)set_option(s+2, next_char);
	  if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
	      string_error(s);
	    break;					/* return -1	*/
	  }
	  if ((augmentation_rate<0) || (augmentation_rate>1.0)) {
	      parse_error(s);
	      break;
	  }
	  if (verbose)
	      printf("GA augmentation rate = %g\n", augmentation_rate);
	  rv = augmentation_rate;
	  break;

	case 'g':			/* generations to run 		*/
	  generations = (int)set_option(s+2, next_char);
	  if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
	      string_error(s);
	    break;					/* return -1	*/
	  }
	  if (generations<=0) {
	      parse_error(s);
	      break;
	  }
	  if (verbose)
	      printf("GA number of generations = %d\n", generations);
	  rv = (double) generations;
	  break;

	case 's':				/* population size	*/
	  rv = (double)(int)set_option(s+2, next_char);
	  if (return_type != NUMBER) {	/* string returned, oops	*/
	      string_error(s);
	      break;					/* return -1	*/
	  }
	  if (rv<0) {				/* out of range		*/
	      parse_error(s);
	      break;
	  }
	  else if (rv > 0)			/* copy to pop_size	*/
	      pop_size = rv;
	  else				/* rv==0, return current value	*/
	      rv = pop_size;

	  if (verbose)
	      printf("population size = %d\n",pop_size);
	  break;

	default:					/* oops		*/
	  *next_char = s+strlen(s);
	  parse_error(s);
	  break;
    } /* end switch */

    return(rv);
} /* end parse_ga() */


/************************************************************************
 *									*
 * parse_jump() compares its first two arguments based on its secondary	*
 * command.  If the comparison is unsuccesful 0 is returned.		*
 * If the comparison is succesful the third argument is used		*
 * as a target address. 						*
 * The current option file is scanned looking for the			*
 * target label.  If the label is found 1 is returned and the		*
 * pointer into the file is left at the line immediately following the	*
 * label.  If the target label is not found -1 is returned		*
 *									*
 ************************************************************************/
double parse_jump(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;
    double first, second;
    char target[60];
    char input_line[200], label[200];
    int offset, first_flag, second_flag;

    first = set_option(s+2, next_char);
    first_flag = return_type;		/* was first argument string?	*/
    if (first_flag == ERROR)
	return(-1);

    second = set_option(*next_char, next_char);
    second_flag = return_type;		/* was second argument string?	*/
    if (second_flag == ERROR)
        return(-1);


    if (first_flag != second_flag) {			/* oops		*/
	printf("unable to compare string and non-string terms\n");
	parse_error(s);
	*next_char = s + strlen(s);
	return(-1);
    }

    target[0] = '\0';			/* terminate so sscanf error 	*
					 * will be noticed		*/

    switch(*(s+1)) {
      case 'e':				/* jump if equal		*/
	if (first_flag == STRING) {	/* comparing strings		*/
	  if (strcmp((char*) (long long) first, (char*) (long long) second)) {
	    rv = 0; break;
	  }
	}
	else {				/* comparing numbers		*/
	  if (first != second)  {
	      rv = 0.0; break;		/* unsuccesful			*/
	  }
	}
	rv = 1; break;			/* succesful			*/

      case 'n':				/* jump if not equal		*/
	if (first_flag == STRING) {	/* strings			*/
	  if (strcmp((char*) (long long) first, (char*) (long long) second) == 0) {
	    rv = 0; break;		/* unsuccesful			*/
	  }
	}
	else {				/* numbers			*/
	  if (first == second)  {
	      rv = 0.0; break;		/* unsuccesful			*/
	  }
	}
	rv = 1; break;			/* succesful 			*/

      case 'l':				/* jump if less than		*/
	if (first_flag == STRING) {	/* comparing strings		*/
	  if (strcmp((char*) (long long) first, (char*) (long long) second) >= 0) {
	    rv = 0; break;		/* unsuccesful 			*/
	  }
	}
	else {				/* numbers			*/
	  if (first >= second)  {
	      rv = 0.0; break;		/* unsucc			*/
	  }
	}
	rv = 1; break;			/* succcess			*/

      case 'g':				/* jump if greater than		*/
	if (first_flag == STRING) {	/* strings			*/
	  if (strcmp((char*) (long long) first, (char*) (long long) second) <= 0) {
	    rv = 0; break;		/* unsuccess			*/
	  }
	}
	else {				/* numbers			*/
	  if (first <= second)  {
	      rv = 0.0; break;		/* unsuccess			*/
	  }
	}
	rv = 1; break;			/* success			*/

      default:						/* oops		*/
	*next_char = s+strlen(s);
	parse_error(s);
	return(-1);;
    }

					/* scan past close parens	*/
    while ((**next_char == ')') || isspace(**next_char))
	(*next_char)++;
					/* scan in target label	*/
    sscanf(*next_char, "%s %n", target, &offset);
    *next_char += offset;

    return_type = NRV;
    if ((rv == -1) || (rv == 0.0))
      return(rv);

    rewind(option_file_ndx->fp);	/* go to beginning of file	*/
    while(1) {				/* read until target label or	*
					 * end of file			*/
      if (fgets(input_line, 200, option_file_ndx->fp) == NULL) {
	printf("target label :%s not found in file %s\n", target,
							option_file_ndx->name);
	*next_char = s + strlen(s);
	return_type = ERROR;
	return(-1);
	break;
      }
      if (input_line[0] == ':') {
	   sscanf(input_line+1, "%s ", label);
	   if (strcmp(target, label)==0)
		return(1);
      }
    } /* end while(1) */
} /* end parse_jump() */



/************************************************************************
 *									*
 * get_log_file() looks up in the logging_files[] array the name 	*
 * passed.  If the logging file is open then the file pointer is 	*
 * returned.  Otherwise the logging file is opened and the new file	*
 * pointer is returned.							*
 *									*
 ************************************************************************/
 FILE *get_log_file(log_file_name)
 char *log_file_name;
 {
     struct file_struct *fsp, *first_opening;

     first_opening = NULL;

    for (fsp=&logging_files[0]; fsp!=logging_files_ndx; fsp++) {
	if (strcmp(fsp->name, log_file_name) == 0) {
			/* name matches, return fp, may open first	*/
	    if (fsp->fp == NULL)
		fsp->fp = open_file(log_file_name, "a");
	    return (fsp->fp);
	}
	if (fsp->fp == NULL)			/* save an empty slot	*/
	    first_opening = fsp;
    }
					/* did not find open file	*/
    if (first_opening) {		/* use existing slot		*/
	strcpy(first_opening->name, log_file_name);
	first_opening->fp = open_file(log_file_name, "a");
	return(first_opening->fp);
    }
					/* use last slot		*/
    if (logging_files_ndx == &logging_files[MAX_FILES]) {
	printf("Error: too many logging files open\n");
	return(NULL);
    }
    logging_files_ndx++;
    strcpy(fsp->name, log_file_name);
    fsp->fp = open_file(log_file_name, "a");
    return(fsp->fp);

} /* end get_log_file() */


/************************************************************************
 *									*
 * parse_logging() opens, closes, and maintains the current logging	*
 * file.								*
 *									*
 ************************************************************************/
double parse_logging(s, next_char)
char *s, **next_char;
{
    struct file_struct *lfp;
    char *log_file_name;
    double rv = 0;

    return_type = NRV;
    switch (*(s + 1)) {
      case 'n':					/* open named file	*/
				/* make sure next argument is a string	*/
	rv = check_for_string(s+2, next_char);

	if (rv <= 0)				/* not a string		*/
	    return(rv);

	log_file_name = (char*) (long long) rv;
	log_file = get_log_file(log_file_name);

	if (log_file != NULL) {		/* log file opened ok		*/
	  if (verbose)
	    printf("log file: %s\n", log_file_name);
	}
	else {				/* problem with log file	*/
	    rv = -1;
	    return_type = ERROR;
	}

	free(log_file_name);			/* free the space	*/
	break;


      case 'f':					/* flush log file	*/
	if (log_file == NULL) {		/* oops, no log file to flush	*/
	    printf("no log file selected\n");
	    *next_char = s + strlen(s);
	    parse_err = 1;
	    return_type = ERROR;
	    break;
	}
	fflush(log_file);			/* flush it		*/
	*next_char = s+2;			/* update input pointer	*/
	break;

      case 's':				/* print remaining values in	*
					 * line in log file		*/
	if (log_file == NULL) {
	    printf("no log file open\n");
	    *next_char = s + strlen(s);
	    parse_err = 1;
	    return_type = ERROR;
	    break;
	}
	*next_char = s + 2;		/* start at next position	*/
	while (**next_char) {	/* as long as there is input to read	*/

	    rv = set_option(*next_char, next_char);
	    if (return_type == STRING) {	/* print the string	*/
		fprintf(log_file, "%s", (char*) (long long) rv);
		free((char*) (long long) rv);
	    }
	    else if (return_type == NUMBER)	/* print the value	*/
		fprintf(log_file, "%g ", rv);
	    else if (return_type == ERROR)
		break;
	}
	fprintf(log_file, "\n");
	if (return_type != ERROR)
	    return_type = NRV;
	break;


      case 'p':				/* plot a point in the log file	*/
	if (log_file == NULL) {
	    printf("no log file open\n");
	    *next_char = s + strlen(s);
	    parse_err = 1;
	    return_type = ERROR;
	    break;
	}

	rv = set_option(s+2, next_char);	/* plot the first point	*/
	if (return_type != NUMBER) {
	    if (return_type != ERROR)
		string_error(s);
	    break;
	}
	fprintf(log_file, "%g ", rv);

	rv = set_option(*next_char, next_char);	/* plot the second point*/
	if (return_type != NUMBER) {
	    if (return_type != ERROR)
		string_error(s);
	    break;
	}
	fprintf(log_file, "%g\n", rv);
	return_type = NRV;
	break;

      case 'c':					/* close logging file	*/
	if (log_file == NULL) {
	    printf("no log file open\n");
	    *next_char = s + strlen(s);
	    parse_err = 1;
	    return_type = ERROR;
	    break;
	}

			/* search through and remove table entry	*/
	for (lfp = &logging_files[0]; lfp < logging_files_ndx; lfp++) {
	    if (lfp->fp == log_file) {
		lfp->fp = NULL;
		if (lfp == (logging_files_ndx - 1))
		    logging_files_ndx--;
		break;
	    }
	}
	rv = (double) fclose(log_file);
	log_file = NULL;
	*next_char = s+2;
	break;

      default:					/* oops			*/
	parse_error(s);
	break;
    } /* end switch */

    return(rv);
} /* end parse_logging() */


/************************************************************************
 *									*
 * parse_measure() measures the entropy and mutual information content	*
 * of CA iterations stored in ca_states.  The calling routine specifies	*
 * where to start, how much to look at, which dimension to group the	*
 * blocks and the size of the blocks					*
 *									*
 ************************************************************************/
double parse_measure(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;
    int i_time, separation, block_size, n;
    char temporal;
    char *msg, *next_item, *limit;
    int line_no;

    switch(*(s + 1)) {
      case 'M':				/* mutual information		*/
	temporal = *(s + 2);		/* third character specifies	*
					 * 's'=space, 't'=time		*/
	if ((temporal != 'l') && (temporal != 's') && (temporal != 't')) {
	    parse_error(s);					/* oops	*/
	    *next_char = s + strlen(s);
	    break;
	}

				/* first parameter is starting position	*/
	i_time = (int) set_option(s+3, &next_item);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		string_error(s);
	    break;					/* return -1	*/
	}
			/* second parameter is distance to block to be 	*
			 * compared with for mutual information		*/
	separation = (int) set_option(next_item, &next_item);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		string_error(s);
	    break;					/* return -1	*/
	}
				/* third parameter is size of blocks 	*
				 * used to generate state configurations*
				 * to measure entropy			*/
	block_size = (int) set_option(next_item, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
	        string_error(s);
	    break;					/* return -1	*/
	}

	if ((next_item == *next_char) || (i_time < 0) ||
					(block_size <= 0) || (separation < 0)) {
	    parse_error(s);	/* error in one or more parameters	*/
	    break;
	}

	if (temporal=='t')
	    limit = &ca_states[(time_step-block_size) * world_size];
	else if (temporal=='s')
	    limit = &ca_states[time_step * world_size];
	else
	    limit = &ca_states[(i_time + 1) * world_size];

				/* calculate mutual information based	*
				 * on parameters read in		*/
	rv = mut_inf_space_time(i_time, separation, limit, block_size,
								temporal=='t');

	if (verbose) {			/* print results		*/
	    if (temporal=='t')
		msg = "temporal";
	    else
		msg = "spatial";
	    printf (
	"%s MI = %g: block size = %d, start time = %d, separation = %d\n",
				msg, rv, block_size, i_time, separation);
	}
	return_type = NUMBER;
	break;					/* end mutual info	*/

      case 'E':						/* Entropy	*/
	switch (*(s + 2)) {
	  case 's':		/* spatial, measure entropy of adjacent	*
				 * blocks from line_no to end of 	*
				 * ca_states				*/
	    line_no = (int) set_option(s+3, next_char);
	    if (return_type != NUMBER) {/* string returned, oops	*/
	      if (return_type != ERROR)
	        string_error(s);
	      break;					/* return -1	*/
	    }

	    block_size = (int) set_option(*next_char, next_char);
	    if (return_type != NUMBER) {/* string returned, oops	*/
	      if (return_type != ERROR)
	        string_error(s);
	      break;					/* return -1	*/
	    }

	    if (block_size <= 0) {
		parse_error(s);
		break;
	    }

	    rv = space_entropy(&ca_states[line_no * world_size],
				(time_step-line_no)*world_size, block_size, 0);
	    if (verbose)
		printf(
		"spatial entropy = %g: block size = %d, start time = %d\n",
						rv, block_size, line_no);
	    return_type = NUMBER;
	    break;

	  case 'l':		/* calculate entropy of adjacent blocks	*
				 * in a single time step		*/
	    line_no = (int) set_option(s+3, next_char);
	    if (return_type != NUMBER) {/* string returned, oops	*/
	      if (return_type != ERROR)
	        string_error(s);
	      break;					/* return -1	*/
	    }

	    block_size = (int) set_option(*next_char, next_char);
	    if (return_type != NUMBER) {/* string returned, oops	*/
	      if (return_type != ERROR)
	        string_error(s);
	      break;					/* return -1	*/
	    }

	    if (block_size <= 0) {
		parse_error(s);
		break;
	    }

	    rv = space_entropy(&ca_states[line_no*world_size], world_size,
								block_size, 0);

	    if (verbose)
		printf(
		"spatial entropy = %g: block size = %d, at time = %d\n",
						rv, block_size, line_no);
	    return_type = NUMBER;
	    break;

	  case 't':		/* calculate entropy of blocks made of	*
				 * sequential entries in time.  (single	*
				 * location)				*/
	    line_no = (int) set_option(s+3, next_char);
	    if (return_type != NUMBER) {/* string returned, oops	*/
	      if (return_type != ERROR)
	        string_error(s);
	      break;					/* return -1	*/
	    }

	    block_size = (int) set_option(*next_char, next_char);
	    if (return_type != NUMBER) {/* string returned, oops	*/
	      if (return_type != ERROR)
	    	string_error(s);
	      break;					/* return -1	*/
	    }

	    if (block_size <= 0) {
		parse_error(s);
		break;
	    }

	    rv=space_entropy(&ca_states[line_no*world_size],
		(time_step - line_no - (block_size - 1)) * world_size,
								block_size, 1);

	    if (verbose)
		printf(
		    "temporal entropy = %g: block size = %d, start time = %d\n",
						rv, block_size, line_no);
	    return_type = NUMBER;
	    break;


	  default:					/* oops		*/
	    *next_char = s + strlen(s);
	    parse_error(s);
	    break;
	} /* end switch *(s+2) */
	break;

      default:						/* oops		*/
	*next_char = s + strlen(s);
	parse_error(s);
	break;

    } /* end switch *(s+1) */

    return(rv);
}


/************************************************************************
 *									*
 * parse_random() initializes random seed, returns random numbers,	*
 * generates random worlds and random rule sets				*
 *									*
 ************************************************************************/
double parse_random(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;
    int offset;
    unsigned init_seed;
    struct timeval tv;


    switch(*(s + 1)) {
      case 'i': 			/* generate random seed		*/
	init_seed = (unsigned) set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		string_error(s);
	    break;					/* return -1	*/
	}

	if (verbose)
	    printf("initialization seed = %d\n", init_seed);
	srandom(init_seed);		/* initialize random numbers	*/
	rv = (double) init_seed;
	return_type = NUMBER;
	break;

      case 't': 			/* return usec time of day	*/
	gettimeofday(&tv, NULL);
	*next_char = s+2;
	rv = (double) tv.tv_usec;
	return_type = NUMBER;
	break;

      case 'r':				/* random rule generation	*/
	*next_char = s+2;
	rv = random_rule_table(&root);	/* returns a string pointer	*/
	break;

      case 'w':				/* random world generation	*/
	*next_char = s+2;
	rv = random_world(&root);	/* returns a string pointer	*/
	break;

      case 'n':				/* return a random number >0 <1	*/
	*next_char = s+2;
	rv = my_random();
	return_type = NUMBER;
	break;

      default:						/* oops		*/
	*next_char = s + strlen(s);
	parse_error(s);
	break;
    }
    return(rv);
} /* end parse_random() */


/************************************************************************
 *									*
 * parse_string() supports various operations on string variables:	*
 * copy substring onto another string, append strings, return an 	*
 * element, measure the length.						*
 *									*
 ************************************************************************/
double parse_string(s, next_char)
char *s, **next_char;
{
    double rv = -1.0;
    int index, length;
    char *string, *insert_string, *sub_string;

    switch (*(s + 1)) {
      case 'i':			/* insert one string into another	*/
					/* position of insert		*/
	index = (int) set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		string_error(s);
	    break;					/* return -1	*/
	}

					/* read source string		*/
	string = (char*) (long long) check_for_string(*next_char, next_char);
	if ((long long) string <= 0)
	    break;

					/* read insert string		*/
	insert_string = (char*) (long long) check_for_string(*next_char, next_char);
	if ((long long) insert_string <= 0) {
	    free(string);
	    break;
	}

	rv = (double)(long long) string;	/* return will be source string	*/

				/* copy until end of either string	*/
	for (string = &string[index]; *string && *insert_string;
						string++, insert_string++) {
	    *string = *insert_string;
	}
	free(insert_string);		/* release space		*/
	return_type = STRING;		/* indicate string return	*/
	break;

      case 'l':				/* measure length of string	*/
	string = (char*) (long long) check_for_string(s+2, next_char);
	if ((long long)string <= 0)
	    break;

	rv = (double) strlen(string);
	free(string);
	return_type = NUMBER;
	break;


      case 's':			/* return a substring of string of 	*
				 * specified length starting at index	*/
	index = (int) set_option(s+2, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		string_error(s);
	    break;					/* return -1	*/
	}

	length = (int) set_option(*next_char, next_char);
	if (return_type != NUMBER) {	/* string returned, oops	*/
	    if (return_type != ERROR)
		string_error(s);
	    break;					/* return -1	*/
	}

	string = (char*) (long long) check_for_string(*next_char, next_char);

	if ((index < 0) || (length <= 0) || ((long long) string <= 0)) {
	    *next_char = s + strlen(s);
	    break;
	}

	sub_string = string + index;		/* pointer to start	*/
	sub_string[length] = '\0';		/* terminate		*/
	rv = get_string(sub_string);		/* copy exact length	*/
	free(string);				/* release space	*/
	return_type = STRING;			/* string returned	*/
	break;


      case 'a':					/* append two strings	*/
						/* get first string	*/
	insert_string = (char*) (long long) check_for_string(s+2, next_char);
	if ((long long) insert_string <= 0)
	    break;
						/* get second string	*/
	string = (char*) (long long) check_for_string(*next_char, next_char);
	if ((long long) string <= 0)
	    break;
					/* make space for result	*/
	if (get_space(&sub_string, strlen(insert_string) + strlen(string) + 1,
								sizeof(char)))
	    break;
	strcpy(sub_string, insert_string);	/* copy first string	*/
	strcat(sub_string, string);		/* append second string	*/
	free(insert_string);			/* release space	*/
	free(string);				/* release space	*/
	rv = (double)(long long) sub_string;
	return_type = STRING;
	break;


      default:						/* oops		*/
	parse_error(s);
	break;
    }
    return(rv);
}


/************************************************************************
 *									*
 * parse_variable() determines if the specified variable exists.  If it	*
 * does not exist a new variable is created.  If the variable is 	*
 * followed by '=' the next command is interpreted and assigned to the	*
 * variable.  The value of the variable is returned.  String variables	*
 * have their string_flag set to 1 and a pointer to the string (cast	*
 * as a double) is returned.						*
 *									*
 ************************************************************************/
double parse_variable(s, next_char)
char *s, **next_char;
{
    char *new_string;
    double rv = -1.0;
    char var_name[20], string[80];
    struct variable *var=NULL;
    int offset, last;

					/* scan in variable name	*/
    if (sscanf(s+1, "%19[_A-Za-z0-9] %n", var_name, &offset)==0) {
	*next_char = s + strlen(s);
	parse_error(s);
	return(rv);
    }

    *next_char = s+1+offset;		/* update pointer to follow	*/

    var = get_var(var_name);		/* locate or create variable	*/
    if (var==NULL)			/* error			*/
	return(rv);

    if (*(s+1 + offset) == '=') {/* assign next value to variable	*/
	rv = set_option(s+1+offset + 1, next_char);
	if (var->string_flag) {
	  if (return_type != STRING) {	/* oops, is a string var	*/
	    printf("error, attempt to assign non-string to string variable\n");
	    return_type = ERROR;
	    parse_err = 1;
	    return(-1);
	  }
	  free((char*) (long long) var->value);	/* release last string	*/
	}
	var->value = rv;		/* assign result of next parse	*/
	if (return_type == STRING) {	/* indicate string returned	*/
	    var->string_flag = 1;
	    rv = get_string((char*)(long long)rv);	/* make copy to return	*/
	}
    }
    else if (var->string_flag) {	/* string variable returned,	*
					 * copy string into new space	*/
	if (get_space(&new_string, strlen((char*)(long long)var->value) + 1,
								sizeof(char)))
	    return(-1);

	strcpy(new_string, (char*)(long long)var->value);

	rv = (double) (long long) new_string;	/* return pointer to new space	*/
	return_type = STRING;
    }
    else {				/* number returned, return it	*/
	rv = var->value;
	return_type = NUMBER;
    }

    return(rv);
} /* end parse_variable() */


/************************************************************************
 *									*
 * get_options() opens up a command file and interprets the lines until	*
 * end of file, ^C, or an error occurs					*
 *									*
 ************************************************************************/
double get_options(s, next_char)
char *s, **next_char;
{
    char *filename;
    char line[4100];		/* must be able to read in MAX_RULE	*/
    FILE *fp;
    double rv;
    struct variable *var, *next_var;
    char *next;

    rv = check_for_string(s+1, next_char);	/* read in filename	*/
    if (rv <= 0)
	return(rv);

    filename = (char*) (long long) rv;
    if (verbose)
	printf("reading commands from file: %s\n", filename);

    option_file_ndx++;		/* prepare to add to list of open files	*/
    if (option_file_ndx >= &option_file[MAX_FILES]) {
	printf("too many files open, did not open %s\n", filename);
	free(filename);
	return_type = ERROR;
	return(-1.0);
    }

						/* open the file	*/
    if ((fp = open_file(filename,"r")) == NULL) {
	free(filename);
	return_type = ERROR;
	return(-1.0);
    }

				/* add filename and pointer to list of	*
				 * open files				*/
    strcpy(option_file_ndx->name, filename);
    free(filename);
    option_file_ndx->fp = fp;

/*    var = var_header;	*/	/* read to end of existing variables	*
				 * all new variables will be deleted	*
				 * when this file is closed		*/
/*    while (var->next)
      var = var->next;	*/

    while (!feof(fp)) {		/* read until end of file		*/
	if (window_count) {	/* kill window processes if exist	*/
	    if (wait3(NULL, WNOHANG | WUNTRACED, NULL) != 0)
		window_count--;
	}

    	if (fgets(line, 4100, fp) == NULL)	/* read up to MAX_RULE	*/
            continue;				/* empty line		*/

	if (int_hit) {				/* ^C pressed		*/
	    break;
	}

	next = line;			/* start at beginning of line	*/
	while (*next) {		/* scan commands until end of line	*/
	    rv = set_option(next, &next);
	    if (return_type)
	       free((char*) (long long) rv);
	}
    }

    if (fclose(fp))			/* close file			*/
        printf("error closing options file \n");

/*    var = var->next;	*/		/* remove all variables added	*
					 * by option file		*/
/*    while (var) {
	next_var = var->next;
	if (var->string_flag)	*/	/* release strings as well	*/
/*	    free((char*) (long long) var->value);
	free(var);
	var = next_var;
    }	*/

    option_file_ndx--;			/* back up open file list	*/

    return_type = NRV;
    return(0);
} /* end get_options() */




/************************************************************************
 *									*
 * set_option() parses an input line from either stdin or a file.  	*
 * Return value is based on operation specified by contents of input	*
 * line.  Pointer *next_char is set to point to first character not 	*
 * examined								*
 * by set_option().  -1.0 is returned in the event of error, pointer to	*
 * string (encoded as a double) may be returned in the event of double	*
 * quotes indicating a string.						*
 *									*
 ************************************************************************/
double set_option(command, next_char)
char	*command, 		/* string containing input command	*/
    	**next_char;		/* address of variable to hold pointer	*
				 * to last character examined		*/
{
    double rv = -1.0;		/* return value, initialize to error	*/
    float f;				/* utility float variable	*/
    int offset;			/* used for calculating next next_char	*/
    double temp;		/* used for temporary storage of double	*/
    char string[MAX_RULES];	/* used for holding string variables	*/
    char var_name[25];		/* name of variable for searching	*/
    struct variable *var_ptr;	/* utility pointer to variable		*/
    int i;			/* utility integer for counting		*/
    char c;			/* utility char for temp storage	*/

		/* if a space allocation or parsing error has occurred	*
		 * read to end of line and return error.  This clears	*
		 * to end of file in the case of input from files, 	*
		 * error indicators are cleared when next input from 	*
		 * stdin is read in.					*/
    if (parse_err || space_err) {
	*next_char = command + strlen(command);
	return(rv);
    }

    fflush(stdout);	/* flush messages to output after each line	*/

    return_type = NRV;

    while(1) {	/* continue scanning input until a command is found	*/

        switch (*command) {
	  case  ' ':				/* scan over spaces	*/
	  case  '\t':				/* scan over tabs	*/
	  case  '(':				/* scan over parens	*/
	  case  ')':				/* scan over parens	*/
	  case  'Z' - 'A':			/* scan over ^Z		*/
	    command++;
	    continue;				/* to top of while(1)	*/


	  case	0:				/* end of string	*/
	  case	'\n':				/* end of line		*/
	  case  ':':			/* unprinted comment line or 	*
					 * target branch for jump	*/
			/* set pointer to end of line and return 0	*/
	    *next_char = command + strlen(command);
	    rv = 0.0;
	    break;


          case	'q':						/* quit	*/
	    if (verbose)
		printf("quit\n");
	    exit(0);


          case	'R':		/* Random, find out which action	*/
	    rv = parse_random(command, next_char);
	    break;


          case	'M':		/* measure statistical property of CA	*/
	    rv = parse_measure(command, next_char);
	    break;


	  case 'C':			/* Cellular Automata rules	*/
	    rv = parse_ca(command, next_char);
	    break;

	  case 'D':				/* Display results	*/
	    rv = parse_display(command, next_char);
	    break;


	  case 'G':				/* GA parameters	*/
	    rv = parse_ga(command, next_char);
	    break;


	  case 'o':				/* read in options file	*/
	    rv = get_options(command, next_char);
	    break;


	  case 'V':		/* 0 argument turns off messages 	*
				 * otherwise turn on messages		*/
	    if ((rv = set_option(command+1, next_char)) != 0)
		verbose = 1;
	    else
		verbose = 0;
	    break;

	  case 'W':
	    rv = set_option(command+1, next_char);
	    if (return_type != NUMBER) {
		parse_error(command);
		break;
	    }
	    if (rv >= 1.0) sleep((int)rv);
	     break;


	  case 'L':			/* logging file commands	*/
	    rv = parse_logging(command, next_char);
	    break;


	  case '$':		/* variable, may set return_type	*/
	      rv = parse_variable(command, next_char);
	      break;


	  case '^':			/* add next two commands	*/
	    temp = set_option(command+1,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    rv = pow(temp, set_option(*next_char,next_char));
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    break;


	  case '+':			/* add next two commands	*/
	    temp = set_option(command+1,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    rv = temp + set_option(*next_char,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    break;


	  case '-':		/* subtract second command from first	*/
	    temp = set_option(command+1,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    rv = temp - set_option(*next_char,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    break;


	  case '*':			/* multiply next two commands	*/
	    temp = set_option(command+1,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    rv = temp * set_option(*next_char,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    break;


	  case '/':		/* divide first command by second	*/
	    temp = set_option(command+1,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    rv = temp / set_option(*next_char,next_char);
	    if (return_type != NUMBER) {
		rv = -1;
		parse_error(command);
		break;
	    }
	    break;


	  case 'E':			/* edit world or rule table	*/
	    rv = parse_edit(command, next_char);
	    break;

	  case 'X':		/* execute CA iteration or GA evolution	*/
	    rv = parse_execute(command, next_char);
	    break;


	  case '.':		/* floating point number, return value	*/
	  case '0':
	  case '1':
	  case '2':
	  case '3':
	  case '4':
	  case '5':
	  case '6':
	  case '7':
	  case '8':
	  case '9':			/* use scanf to do conversion	*/
	    sscanf(command, "%f %n", &f, &offset);
	    *next_char = command + offset;
	    rv = (double) f;
	    return_type = NUMBER;
	    break;


	  case 'J':				/* conditional jump	*/
	    rv = parse_jump(command, next_char);
	    break;


	  case '!':			/* execute shell command	*/
	    if (verbose) {
		printf("executing shell command: %s\n", command+1);
		fflush(stdout);		/* output from shell command	*
					 * can come before this message	*/
	    }
	    *next_char = command + strlen(command);
	    rv = (double) system(command+1);
	    return_type = NUMBER;
	    break;


	  case 'S':				/* index into string	*/
	    rv = parse_string(command, next_char);
	    break;


	  case '"':			/* string constant, used for	*
					 * filenames, rule tables, 	*
					 * world configurations.	*/

	    i = 0;
				/* copy chars into string[] until end	*
				 * of string or newline			*/
	    while ((c = *(command + 1 + i)) && (c != '"') && (c!='\n')) {
		string[i] = c;
		i++;
	    }
	    string[i] = '\0';			/* terminate string	*/
				/* search for string already existing 	*
				 * or create a new permanent place for	*
				 * string				*/
	    rv = get_string(string);


	    if (c=='"') 		/* start next command after	*
					 * closing quote if it exists	*/
	        *next_char = command + i + 2;
	    else
		*next_char = command + i + 1;

	    sscanf(*next_char, " %n", &offset);
	    *next_char += offset;

	    return_type = STRING;
	    break;


	  case '%':		/* print a comment onto standard output	*
				 * print value or string associated w/	*
				 * variables starting w/ $		*/

	    if (!verbose) {		/* only print if in verbose 	*
					 * mode 'V 0/1'			*/
	      *next_char = command + strlen(command);
	      break;
	    }

	    i = 1;
	    while(c = *(command + i)) {		/* do for each char	*/
	      if (c == '$') {		/* if variable, print value	*/

		var_name[0] = '\0';		/* get variable name	*/
		i++;
		sscanf(command+i, "%19[_A-Za-z0-9]%n", var_name, &offset);

		var_ptr = get_var(var_name);	/* look up variable	*/

					/* if string, print in "..."	*/
		if (var_ptr->string_flag)
		    printf("\"%s\"", (char *)(long long)var_ptr->value);
		else				/* just print value	*/
		    printf("%g", var_ptr->value);
		i += offset;			/* skip string name	*/
	      }
	      else if (c == '\n')			/* end of line	*/
		break;
	      else {			/* copy character to output	*/
		printf("%c", *(command+i));
	        i++;					/* next char	*/
	      }
	    }
	    printf("\n");			/* terminate line	*/
	    *next_char = command + strlen(command);
	    rv = 0.0;					/* return 0	*/
	    break;


	  case 'h':				/* print help file	*/
	  case '?':
	    *next_char = command + 1;	/* start next read past this	*/
/* 4-17-91 Thomas E. Kammeyer -- Changed the path name on following line so
 * I could have a help file available!
 */
	    rv = (double) system("more TEK/ca.help");
	    break;


	  default:						/* oops	*/
	    *next_char = command+strlen(command);
	    parse_error(command);
	    break;
	}
	break;
    }
    return(rv);
} /* end set_option() */

