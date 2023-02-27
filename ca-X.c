
#define INT_EXT extern
#include "ca.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "ca_cursor.h"
#include "ca_cursor_mask.h"

char ca_string[] = "CA Evolution";
XRectangle *state_pix;

#define ca_cell_width 2
#define ca_cell_height 2


disp_window()
{
    Display *disp;
    Window ca_window;
    GC gc;
    Cursor curs;
    Pixmap source, mask;

    XEvent event;
    KeySym key;

    XColor black_color, white_color;

    XSizeHints hint;

    int screen;
    unsigned long white, black, foreground;
    int i,j,first;
    int done;
    char text[10];


    if (get_space(&state_pix, world_size, sizeof(XRectangle))) {
	return(1);
    }

    disp = XOpenDisplay("");
    for (i=0; i<world_size; i++) {
	state_pix[i].x = i*ca_cell_width;
	state_pix[i].width = ca_cell_width;
	state_pix[i].height = ca_cell_height;
    }

    screen = DefaultScreen(disp);
    white = WhitePixel (disp, screen);
    black = BlackPixel (disp, screen);

    black_color.red = 0; black_color.green = 0; black_color.blue = 0;
    black_color.flags = DoRed | DoGreen | DoBlue;
    white_color.red = 65535; white_color.green = 65535; white_color.blue = 65535;
    white_color.flags = DoRed | DoGreen | DoBlue;


    hint.x = 200; hint.y = 0;
    hint.width = world_size * ca_cell_width; hint.height = max_time * ca_cell_height;
    hint.flags = PPosition | PSize;

    ca_window = XCreateSimpleWindow (disp, DefaultRootWindow (disp),
	hint.x, hint.y, hint.width, hint.height, 1, white,
	black);

    XSetStandardProperties (disp, ca_window, ca_string, ca_string, 
	None, glob_argv, glob_argc, &hint);
    
    source = XCreateBitmapFromData (disp, ca_window, ca_cursor_bits, 
		ca_cursor_width, ca_cursor_height);
    mask = XCreateBitmapFromData (disp, ca_window, ca_cursor_mask_bits, 
		ca_cursor_mask_width, ca_cursor_mask_height);
    curs = XCreatePixmapCursor (disp, source, mask, &black_color, &white_color, 
		ca_cursor_x_hot, ca_cursor_y_hot);

/*
    curs = XCreateFontCursor(disp, XC_gumby);
*/
    XDefineCursor (disp, ca_window, curs);

    gc = XCreateGC (disp, ca_window, 0, 0);
    XSetBackground (disp, gc, white);

    XSelectInput (disp, ca_window, ButtonPressMask | KeyPressMask | ExposureMask);

    XMapRaised (disp, ca_window);

    done = 0;

    while (done==0) {
	XNextEvent (disp, &event);
	switch (event.type) {
	  case Expose:
	    if (event.xexpose.count == 0)
    		for (i=0; i<time_step; i++) {
		    char state;

		    if (ca_states[i * world_size]==0)
			foreground = black;
		    else
			foreground = white;

		    first = 0;
		    j=0;
		    while (j!=world_size) {
			if (int_hit) {
			    done = 1;
			    break;
			}
			if (ca_states[i * world_size + j])
			  while ((j<world_size) &&
				(ca_states[i * world_size + j] != 0)) {
			    state_pix[j].y = i*ca_cell_height;
			    j++;
			  }
			else 
			  while ((j<world_size) &&
				(ca_states[i * world_size + j] == 0)) {
			    state_pix[j].y = i*ca_cell_height;
			    j++;
			  }

    		        XSetForeground (disp, gc, foreground);
			XFillRectangles(disp, ca_window, gc, &state_pix[first], 
								j - first);
						/* switch foregrounds	*/
			foreground = (foreground==white)? black:white;
			first = j;
		    }
		}

	    break;
	  
	  case MappingNotify:
	    XRefreshKeyboardMapping (&event);
	    break;

	  case ButtonPress:
	    done = 1;
	    break;

	  case KeyPress:
	    i = XLookupString (&event, text, 10, &key, 0);
	    if (i == 1 && text[0] == 'q') done = 1;
	    break;
	} /*switch*/
    } /*while*/

    XFreeGC (disp, gc);
    XDestroyWindow (disp, ca_window);
    XCloseDisplay (disp);

    free(state_pix);
    return(0);
} /**/




