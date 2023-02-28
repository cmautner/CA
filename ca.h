
#include <stdio.h>		/* file i/o routines 		*/
#include <ctype.h>		/* macros: isdigit(),...	*/
#include <math.h>		/* function exp()		*/
#include <signal.h>
#include <sys/time.h>
#include <sys/wait.h>

INT_EXT int world_size;
#define MAX_WORLD 512
#define MAX_TIME 512
#define MAX_RULES 4096
#define MAX_STATES 8

struct CA {
    int	states;
    int neighborhood;
    float lambda;
    char *world;
    char rule_table[MAX_RULES];
};

INT_EXT struct CA root;		/* root CA to return to when	*
				 * exploration gets stuck	*/
INT_EXT int pop_size;		/* population size for GA	*/
INT_EXT float mutation_rate;	/* GA mutation rate		*/
INT_EXT float augmentation_rate;/* GA augmentation rate		*/
INT_EXT int generations;	/* number of generations to run	*
				 * simulations for		*/

INT_EXT FILE *log_file;		/* file to store results in	*/
INT_EXT int int_hit;		/* interrupt key struck		*/

INT_EXT int max_time;
INT_EXT char *ca_states;
INT_EXT int time_step;

INT_EXT int glob_argc;
INT_EXT char **glob_argv;

double set_option();
double get_options();
double space_entropy();
double mut_inf_space_time();
int get_space(char **pointer, int space_size, int object_size);
void parse_error(char *s);
void string_error(char *s);


#define MAX_FILES 32
INT_EXT struct file_struct {
    char name[30];
    FILE *fp;
} option_file[MAX_FILES], *option_file_ndx, 
  logging_files[MAX_FILES], *logging_files_ndx;

INT_EXT struct variable {
    char name[20];
    int string_flag;
    double value;
    struct variable *next;
} *var_header;

struct variable *get_var();
INT_EXT int parse_err;
INT_EXT int space_err;
INT_EXT int verbose;
double get_string();
double random_rule_table();
double random_world();
char *edit_line();

void interrupt();
FILE *open_file();
#define ERROR 0
#define NUMBER 1
#define STRING 2
#define NRV 3			/* no return value */
INT_EXT int return_type;
double my_random();
int window_count;
