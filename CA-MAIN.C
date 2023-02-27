/* file: ca-evol.c							*/

#define INT_EXT /* used in ca.h, defined as extern in all other files */

#include "ca.h"


/************************************************************************
 *									*
 *									*
 ************************************************************************/

main(argc, argv) 
int	argc;			/* number of items on command line	*/
char	*argv[];		/* array of pointers to command	line 
				 * items				*/
{
    
    int		i;		/* counter				*/
    long 	j,		/* counter				*/
		total_loops;	/* total epochs completed		*/
	
    char 	s[100];		/* array to hold input string for user	*/
    char	*parameter,*next;
    FILE 	*fp;		/* pointer to input file		*/
    double	rv;		/* value returned from set_option	*/

    if (defaults())		/* set global values to defaults, may	*
    				 * be overridden by set_option()	*/
        exit(1);
        
    signal(SIGINT, interrupt);	/* set vector for ^c			*/

    glob_argc = argc;		/* hold for geometry, etc. for X window	*/
    glob_argv = argv;

    printf("\n");
    while(1) {			/* repeat until user exits from program	*/
	if (window_count) {	/* kill window processes if exist	*/
	    if (wait3(NULL, WNOHANG | WUNTRACED, NULL) != 0)
		window_count--;
	}

	parse_err = space_err = 0;
	int_hit = 0;
	printf("Option(h=help): ");

	if (gets(s) == NULL) break;
	next = s;
	while (*next) {
	    if (int_hit) break;
	    rv = set_option(next,&next);
	    if (return_type == STRING)
		free((char *) (int) rv);
	}
    }   			/* end while(1) to get key input*/
}


defaults()
{
    int_hit = 0;
    log_file = NULL;
    verbose = 1;
    parse_err = 0;
    space_err = 0;
    time_step = 0;
    world_size = 79;
    max_time = 200;
    window_count = 0;

    if (get_space(&var_header, 1, sizeof(struct variable)) ||
	get_space(&ca_states, world_size*MAX_TIME, sizeof(char)) ||
	get_space(&root.world, world_size, sizeof(char))) {
      return(-1);
    }

    root.neighborhood = 1;
    root.states = 2;
    root.lambda = 0.3;

    logging_files_ndx = &logging_files[0];
    option_file_ndx = &option_file[0];
    return 0;
}

