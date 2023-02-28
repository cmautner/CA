
/* file: ca-eval.c							*/
#define INT_EXT extern

#include <stdlib.h>

#include "ca.h"



/************************************************************************
 *									*
 * iterate_ca() copies the world state from the passed CA structure	*
 * into the ca_states table.  The CA world is then iterated through	*
 * max_time time steps (each iteration a new line in the ca_states	*
 * table).  The final world state after max_time iterations is copied 	*
 * back into the world array specified by the passed structure		*
 *									*
 ************************************************************************/
void iterate_ca(source_ca)
struct CA *source_ca;
{
    int i, j, k, rule_index, states, nghbr;
    char *ca_state_ptr, *ptr;

    nghbr = source_ca->neighborhood;		/* neighborhood size	*/
    states = source_ca->states;			/* number of states	*/

			/* copy world state from source_ca to ca_states	*/
    ca_state_ptr = ca_states;
    ptr = source_ca->world;
    for (i=0; i< world_size; i++, ca_state_ptr++, ptr++)
        *ca_state_ptr = *ptr;

					/* do for each time step	*/
    for (time_step=1, ca_state_ptr=ca_states, ptr=&ca_states[world_size]; 
		time_step<max_time; time_step++, ca_state_ptr += world_size) {

				/* do for each position in world	*/
        for (i=0; i<world_size; i++, ptr++) {

			/* use rules to generate next CA from prev CA	*
			 * based on states of neighbors on both sides	*/
	    for (rule_index = 0, j=-nghbr; j<=nghbr; j++) {
	        rule_index *= states;
						/* wrap around world	*/
		if ((k = i + j) < 0)
	            rule_index += ca_state_ptr[world_size + k];
		else if (k >= world_size)
		    rule_index += ca_state_ptr[k - world_size];
		else
		    rule_index += ca_state_ptr[k];
	    }
		    /* look up rule and place new state at current spot	*/
	    *ptr = source_ca->rule_table[rule_index];
        }
    }
					/* copy last state into world	*/
    ca_state_ptr = &ca_states[(time_step-1) * world_size];
    ptr = source_ca->world;
    for (i=0; i< world_size; i++, ca_state_ptr++, ptr++)
       *ptr = *ca_state_ptr ;
}    



/************************************************************************
 *									*
 * random_rule_table() uses the lambda value to generate a random rule	*
 * table.  Returns a string of the new rule table.			*
 *									*
 ************************************************************************/
double random_rule_table(ca_ptr)
struct CA *ca_ptr;
{
    int table_size,i;
    float ran_val, lambda;
    char rule[MAX_RULES];
    double rv;

    table_size = 			/* number of rules in table	*/
	(int) pow((double)ca_ptr->states, (double)(2*ca_ptr->neighborhood + 1));

    ca_ptr->rule_table[0] = 0;	/* quiescent states to quiescence	*/
    rule[table_size - 1] = '0';

    lambda = ca_ptr->lambda;			/* set lambda parameter	*/

    for(i=1; i<table_size; i++) {/* do for each element in rule table	*/
	ran_val = my_random();			/* 0 < ran_val < 1	*/
	if (ran_val > lambda)	/* lambda is % of non-zero rules	*/
	    ca_ptr->rule_table[i] = 0;
	else 		/* fill in random state between 1 and states-1	*/
	    ca_ptr->rule_table[i] = 
			((int)((ca_ptr->states - 1) * ran_val / lambda)) + 1;
					/* add state to return string	*/
	rule[table_size - i - 1] = ca_ptr->rule_table[i] + '0';
    }
    rule[table_size] = '\0';		/* terminate return string	*/
    return_type = STRING;		/* indicate string returned	*/

    rv = get_string(rule);	/* generate new space for string	*/
    if ((rv > 0) && (verbose))
	printf("Random Rule = \"%s\"\n", rule);
    return(rv);
} /* end random_rule_table() */


/************************************************************************
 *									*
 * random_world() generates a random set of states and stores them in	*
 * the passed CA structure's world array.  A string representing the 	*
 * states is returned							*
 *									*
 ************************************************************************/
double random_world(ca_ptr)
struct CA *ca_ptr;
{
    int i;
    char world[MAX_WORLD];				/* string	*/
    double rv;

    for (i=0; i<world_size; i++) {/* do for each position in world	*/
	ca_ptr->world[i] = rand() % ca_ptr->states;
	world[i] = ca_ptr->world[i] + '0';	/* add to string	*/
    }
    world[i] = '\0';				/* terminate string	*/
    return_type = STRING;		/* indicate string returned	*/

    rv = get_string(world);		/* create permanent data	*/
    if ((rv > 0) && verbose)
	printf("Random World = \"%s\"\n", world);
    return (rv);
} /* end random_world() */


void run_evolve()
{ }
