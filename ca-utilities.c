/* file: ca-utilities.c						*/

#define INT_EXT extern 		/* all variables in ca.h are now*\
				 * extern to this file		*/

#include <stdlib.h>
#include <string.h>

#include "ca.h"

/************************************************************************
 *									*
 * interrupt() called when ^C (SIGINT) struck, sets global indication.	*
 *									*
 ************************************************************************/
void interrupt()
{
	signal(SIGINT, SIG_IGN);
   int_hit = 1;
	signal(SIGINT, interrupt);
} /* end interrupt() */


/************************************************************************
 *									*
 * open_file() opens named file in indicated mode, returns pointer to 	*
 * file or NULL in the event of error.  Prints error message to stdout	*
 *									*
 ************************************************************************/
FILE *open_file(filename, mode)
char *filename, *mode;
{
    FILE *fp;

    if ((fp = fopen(filename, mode)) == NULL)
	printf("error opening file %s\n",filename);

    return(fp);
} /* end open_file() */


/************************************************************************
 *									*
 * parse_error() prints message to stdout and sets global indication.	*
 *									*
 ************************************************************************/
void parse_error(s)
char *s;
{
    printf("unable to parse input=> %s\n", s);
    parse_err = 1;
    return_type = ERROR;
} /* end parse_error() */

/************************************************************************
 *									*
 * get_var() searches through linked list of variable names from input	*
 * if passed var_name is found pointer to structure containing it is 	*
 * returned.  If var_name is not found, new structure is added to list	*
 * and initialized to value = 0.0, new structure is placed at end of 	*
 * linked list, and pointer to new structure is returned.		*
 *									*
 ************************************************************************/
struct variable *get_var(var_name)
char *var_name;				/* name of variable to look for	*/
{
    struct variable 
	*var_ptr, *last_ptr;	/* pointers for traversing linked list	*/

    var_ptr = var_header;		/* start at beginning of list	*/
    while (var_ptr) {	/* continue until NULL pointer encountered	*/
      if (strcmp(var_name, var_ptr->name)==0)
	  return(var_ptr);			/* variable names match	*/

      last_ptr = var_ptr;			/* save previous entry	*/
      var_ptr = var_ptr->next;			/* switch to next entry	*/
    }
	/* entire list traversed and variable name was not found create	*
	 * new entry in linked list					*/
    if (get_space((char**)&var_ptr, 1, sizeof(struct variable))) {
      return(NULL);
    }
    last_ptr->next = var_ptr;			/* extend linked list	*/
    strcpy(var_ptr->name, var_name);	/* add name to structure	*/
    return(var_ptr);			/* return pointer to structure	*/
} /* end get_var() */


/************************************************************************
 *									*
 * get_string() allocates enough space to store the string passed.	*
 * The string is copied into the new space and a pointer cast as a	*
 * double is returned.							*
 *									*
 ************************************************************************/
double get_string(str)
char *str;				/* string to make space for	*/
{
    int length;		/* length of str, used for allocating space	*/
    char *rv;

    length = strlen(str) + 1;	/* length of array to hold string	*/
			/* allocate space for structure and for string	*/
    if (get_space(&rv, length, sizeof(char))) {
      return(0);
    }

    strcpy(rv, str);				/* copy string to space	*/
    return((double) (long long) rv);		/* return pointer to space	*/
} /* end get_string() */


/************************************************************************
 *									*
 * get_space() allocates specified space and assigns it to passed 	*
 * pointer.  Prints message and	returns 1 if error, returns 0 o/w	*
 *									*
 ************************************************************************/
int 
get_space(pointer, space_size, object_size)
char **pointer;				/* pointer to assign space to	*/
int space_size, 				/* number of elements	*/
    object_size;				/* size of elements	*/
{
    if ((*pointer = (char *)calloc(space_size, object_size)) == NULL) {
	printf("error allocating space\n");
	space_err = 1;		/* set global indication of error	*/
	return_type = ERROR;
	return(1);
    }
    return(0);
} /* end get_space() */


/************************************************************************
 *									*
 * my_random() returns a double between 0 and 1 with a precision of	*
 * one part in a million.						*
 *									*
 ************************************************************************/
double my_random()
{
    return(((double) (rand() % 1000000))/1000000.0);
}


/************************************************************************
 *									*
 * string_error() prints error message					*
 *									*
 ************************************************************************/
void string_error(s)
char *s;
{
    printf(
    "Error: string constant or NULL returned when expecting number=> %s\n",s);
    parse_err = 1;
    return_type = ERROR;
}

