/* file: ca-entropy.c							*/

#define INT_EXT extern
#include "ca.h"


/************************************************************************
 *									*
 * entropy() returns the entropy = -sum(p_i * log(p_i)) of the passed	*
 * parameters.  dist_size is i, prob_dist is the probability of each i	*
 * prob_dist should total to 1 if all i elements are added		*
 *									*
 ************************************************************************/
double entropy(prob_dist, dist_size)
double *prob_dist;
int dist_size;
{
    double rv = 0.0;
    int i;

    for (i=0; i<dist_size; i++)
	if (prob_dist[i] != 0) {
	    rv += prob_dist[i] * log(prob_dist[i]);
	}

    return -rv;
} /* end entropy() */



/************************************************************************
 *									*
 * space_entropy computes the value the entropy of ca_states[].  The	*
 * particular area of ca_states[] to examine is defined as passed	*
 * parameters.  space points to the start of the area, it is assumed	*
 * to be the start of a time step of ca_states[].  size is the number	*
 * of locations to include in the sample, e.g. size=world_size would	*
 * return the entropy of a single time step.  block_size is the number	*
 * adjacent locations to look at in order to determine the probability	*
 * distribution, e.g. block_size=3 in a 2-state CA would yield 2^3=8	*
 * items in the probability distribution.  temporal defines what 	*
 * direction adjacency lies in, i.e. temporal=0 means that successive	*
 * locations in ca_states are used and wraparound at the borders is	*
 * necessary. temporal=1 means that the probability distribution is	*
 * to be generated by looking at a single location at different points	*
 * in time.								*
 *									*
 ************************************************************************/
double space_entropy(space, size, block_size, temporal)
char *space;
int size, block_size, temporal;
{
    int i, state, state_combos;
    int n, p, ndx;
    char *block_ptr;
    double *state_dist;

		/* calculate the size of the probability distribution	*/
    state_combos = (int) pow((double) root.states, (double) block_size);

				/* generate space for the prob dist	*/
    if (get_space(&state_dist, state_combos, sizeof(double))) {
	space_err = 1;
	return(-1);
    }

					/* do for every element in size	*/
    for (i=0; i<size; ) {
				/* do for every element in the world	*/
      for (p=0; (p<world_size) && (i<size); p++, i++, space++) {
				/* start at the beginning of the block 	*
				 * and generate the full state		*/
	for (n = 0, block_ptr = space, ndx = 0; n<block_size; n++) {
	    ndx *= root.states;
	    ndx += *block_ptr;
	    
	    if (temporal)		/* check this world position	*
					 * at next time step		*/
		block_ptr += world_size;
	    else {			/* check next world position	*/
		block_ptr++;
		if ((p + n + 1) == world_size)		/* wrap around	*/
		    block_ptr -= world_size;
	    } /* end not temporal */
	} /* end for each position in block */

	state_dist[ndx]++;	/* increment this prob dist item	*/
      } /* end for each element in world */
    } /* end for each element in size */

    for (i=0; i<state_combos; i++) {		/* normalize prob dist	*/
	state_dist[i] /= (double)size;
    }

    free(state_dist);				/* clean up		*/
    return entropy(state_dist, state_combos);	/* do the calculation	*/
} /* end space_entropy() */



/************************************************************************
 *									*
 * mut_inf_space_time() computes the mutual information of two blocks	*
 * in ca_states separated by distance separation.  i_time is identical	*
 * to start in space_entropy() (above) and block_size and temporal are	*
 * also the same.  size is assumed to be from i_time to the end of 	*
 * ca_states[].	mutual information between two blocks is defined to be	*
 * the sum of the entropies of each block (which will be the same in 	*
 * this case) minus the entropy of the two blocks treated as a single	*
 * unit									*
 *									*
 ************************************************************************/
double mut_inf_space_time(i_time, separation, limit, block_size, temporal)
int i_time, separation, block_size, temporal; 
char *limit;
{
    double rv, *i_prob_dist, *j_prob_dist, *ij_prob_dist, sample_size=0.0;
    int n, m, p, i_ndx, j_ndx, ts /*table_size*/, j_time, j_pos, block_step;
    char *i_source, *j_source, *i_ptr, *j_ptr;

			/* calculate size of prob dist for one block	*/
    ts = (int) pow((double)root.states, (double)block_size);
				/* j_time is line of second block	*/
    j_time = i_time + (separation / world_size);
				/* j_pos is offset within line of 2nd 	*/
    j_pos = separation % world_size;
    
				/* generate space to hold 3 prob dists	*/
    if (get_space(&i_prob_dist, ts, sizeof(double)) ||
        get_space(&j_prob_dist, ts, sizeof(double)) ||
        get_space(&ij_prob_dist, ts*ts, sizeof(double))) {
    	free(i_prob_dist);				/* oops		*/
    	free(j_prob_dist);
    	free(ij_prob_dist);
	space_err = 1;
	return(-1.0);
    }

				/* positions currently being examined	*/
    i_source = &ca_states[i_time * world_size]; 
    j_source = &ca_states[j_time * world_size + j_pos];

    if (temporal) {		/* judge successive time steps		*/
	block_step = world_size;
    }
    else {				/* judge adjacent positions	*/
	block_step = 1;
    }

			/* make sure neither pointer exceeds limit	*/
    while ((i_source<=limit) && (j_source<=limit)) {
			/* look at each block in the current world	*/
      for (p=0; p<world_size; p++, i_source++, j_source++, j_pos++, 
							sample_size++) {

	if (j_pos == world_size) {		/* wrap around j_source */
	    j_pos=0;
	    j_source = &ca_states[j_time * world_size];
	}

					/* do for each element in block	*/
	for (n=0, i_ptr = i_source, j_ptr = j_source, i_ndx = 0, j_ndx = 0; 
		     n<block_size; n++, i_ptr+=block_step, j_ptr+=block_step) {

	    if (!temporal) {
		if ((j_pos + n) == world_size)
		    j_ptr -= world_size;		/* wrap around	*/
		if ((p + n) == world_size)
		    i_ptr -= world_size;
	    }

	    i_ndx *= root.states;		/* increment index	*/
	    i_ndx += *i_ptr;

	    j_ndx *= root.states;		/* increment index	*/
	    j_ndx += *j_ptr;
	} /* end for each element in block */

	i_prob_dist[i_ndx]++;			/* increment prob dist	*/
	j_prob_dist[j_ndx]++;			/* increment prob dist	*/
					/* increment joint prob dist	*/
	ij_prob_dist[(i_ndx * ts) + j_ndx]++;
      } /* end for each position in world */

      j_time++; 				/* go to next time step	*/
      j_source = &ca_states[j_time * world_size + j_pos];	
    } /* end while not at the end of ca_states[] */

    for (n=0; n<ts; n++) {		/* normalize prob dist		*/
	i_prob_dist[n] /= sample_size;
	j_prob_dist[n] /= sample_size;
    }

    for (n=0; n<ts*ts; n++)		/* normalize prob dist		*/
	ij_prob_dist[n] /= sample_size;

    rv = entropy(i_prob_dist, ts);/* calculate individual entropies	*/
    rv += entropy(j_prob_dist, ts);
					/* calculate joint entropy	*/
    rv -= entropy(ij_prob_dist, ts * ts);

    free(i_prob_dist);					/* clean up	*/
    free(j_prob_dist);
    free(ij_prob_dist);
    return(rv);
} /* end mut_inf_space_time() */
