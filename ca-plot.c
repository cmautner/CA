/* UNIX(tm) plot support by Thomas Kammeyer
 */

#define INT_EXT extern
#include "ca.h"
#include <string.h>
#include <math.h>

#define SPACE_X 3120
#define SPACE_Y 3120

#define FILL 1
#define NOFILL 0

void box(x1, y1, x2, y2, fflag)
	int x1, y1, x2, y2, fflag;
{
#if 0  // removed because plot is no longer a curses library (see Makefile)
	move(x1, y1);
	cont(x2, y1);
	cont(x2, y2);
	cont(x1, y2);
	cont(x1, y1);

	if (fflag) {
		/* The following can take to much time, so a pixel
		 * that is filled in just has an X drawn in the box...
		 * for now.
		int i = y1;
		for (; i < y2; i++) line(x1, i, x2, i);
		 */
		line(x1, y1, x2, y2);
		line(x1, y2, x2, y1);
	}
#endif
}

int disp_plot()
{
#if 0  // removed because plot is no longer a curses library (see Makefile)
	int i, j, min_x, max_x, min_y, max_y, xoff, scale;

	/* Open and define plotting space and clear screen */
	openpl();

	int_hit = 0;

	min_x = min_y = 0;
	max_x = SPACE_X * 5 / 4 - 1;	/* 4014 screen is bigger than space...
					 * space() defines a square subarea.
					 */
	max_y = SPACE_Y - 1;

	space(0, SPACE_Y, SPACE_X, 0);
	erase();

	xoff = max_x / 4;
	i = (int) ((double)max_x/(double)world_size/4.0*3.0);
	j = (int) ((double)max_y/(double)time_step);
	scale = (i < j) ? i : j ;
	box(xoff, 0, xoff+scale*world_size, scale*time_step, NOFILL);


	/* Draw the CA states */
	for (i=0; (!int_hit) && (i<time_step); i++) {
		for (j = 0; j < world_size; j++) {
			if (ca_states[i * world_size + j]) {
				box(xoff+scale*j, scale*i,
				    xoff+scale*(j+1), scale*(i+1), FILL);
			}
		}
		fflush(stdout);
	}

	move(0, 0);

	/* Close plotting space. */
	fflush(stdout);
	closepl();

	int_hit = 0;
#endif
	return(0);
}

/* End of file.
 */
