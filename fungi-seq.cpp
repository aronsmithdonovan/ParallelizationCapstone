/*******************************************************************************************
 * fungi-seq.cpp
 *******************************************************************************************
 * 
 * sequentially simulates the growth of a mushroom network in a patch of grass
 * 
 * created by Aron Smith-Donovan using code written by Libby Shoop as reference
 * 
 * based on a project description posited in "Introduction to Computational Science:
 *      Modeling and Simulating for the Sciences" by Angela B. Shiflet and George W Shiflet
 *      
 * 
*/

/* LIBRARIES */
    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <unistd.h>
    #include <cstdlib>
    #include <iostream>
    #include <trng/yarn2.hpp>
    #include <trng/uniform01_dist.hpp>
    #include <locale.h>
    #include <wchar.h>
    #include "seq_time.h"  // Libby's timing function that is similar to omp style

/* UNIVERSAL CONSTANTS */
    // probability values for state changes
    #define probSpore 0.001             // probability that a site initially is SPORE
    #define probSporeToYoung 0.25       // probability that a SPORE will become YOUNG at the next time step
    #define probSpread 0.6              // probability that a EMPTY with a neighbor that is YOUNG will become YOUNG at the next time step
    #define probMaturingToMushrooms 0.7 // probability that a MATURING will become MUSHROOMS at the next time step (otherwise it becomes OLDER)
    #define probDepletedToSpore 0.0001  // probability that a DEPLETED will become SPORE at the next time step
    #define probDepletedToEmpty 0.5     // probability that a DEPLETED will become EMPTY at the next time step

    // cell states
    #define EMPTY 0      // empty ground containing no spore or hyphae
    #define SPORE 1      // contains at least one spore
    #define YOUNG 2      // young hyphae that cannot form mushrooms yet 
    #define MATURING 3   // maturing hyphae that cannot form mushrooms yet
    #define MUSHROOMS 4  // older hyphae with mushrooms
    #define OLDER 5      // older hyphae with no mushrooms
    #define DECAYING 6   // decaying hyphae with exhausted nutrients
    #define DEAD 7       // newly dead hyphae with exhausted nutrients
    #define DEADER 8     // hyphae that have been dead for a while
    #define DEPLETED 9   // area whose nutrients have previously been depleted by fungal growth
    #define INERT 10     // inert area where plants cannot grow

/* FUNCTION DECLARATIONS */
void getArguments(int argc, char *argv[], int * ROWS, int * COLUMNS, int * TIME_STEPS);
void allocateGrid(int ***grid, int * ROWS, int * COLUMNS, int * current_row);
void initializeGrid(int ***grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column, double * prob, trng::yarn2 * yarn, trng::uniform01_dist<> * uniform);
void mushrooms(int ***current_grid, int ***next_grid, int * ROWS, int * COLUMNS, int * TIME_STEPS, int * current_row, int * current_column, int * current_time_step, int * neighbor_row, int * neighbor_column, int * current_value, double * prob, trng::yarn2 * yarn, trng::uniform01_dist<> * uniform);
void copyGrid(int ***current_grid, int ***next_grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column);
int check_neighbors(int ***current_grid, int * current_row, int * current_column, int * neighbor_row, int * neighbor_column);
void deallocateGrid(int ***grid, int * ROWS, int * current_row);
void print_number_grid(int ***grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column);
void print_colorful_grid(int ***grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column, int * current_value);
void reset_color();
void black();
void red();
void green();
void brown();
void grey();
void purple();

/* main */
int main(int argc, char **argv){

    // declare variables
    double start_time, end_time, total_time;  // hold timer values
    int ROWS, COLUMNS, TIME_STEPS;  // hold command line arguments
    int **current_grid;  // grid at current time step
    int **next_grid;  // grid at next time step
    int current_row, current_column;  // grid cell counters
    int current_time_step;  // time step counter
    int neighbor_row, neighbor_column;  // check_neighbors() counters
    int current_value;  // hold grid print values
    double prob;  // stores randomly generated probability values

    // initialize random number engine
    trng::yarn2 yarn;  // create engine object
    trng::uniform01_dist<> uniform;  // create distribution fxn

    // parse command line arguments
    getArguments(argc, argv, &ROWS, &COLUMNS, &TIME_STEPS);

    // start timing
    start_time = c_get_wtime();

    // allocate grids
    allocateGrid(&current_grid, &ROWS, &COLUMNS, &current_row);
    allocateGrid(&next_grid, &ROWS, &COLUMNS, &current_row);

    // initialize current_grid
    initializeGrid(&current_grid, &ROWS, &COLUMNS, &current_row, &current_column, &prob, &yarn, &uniform);

    // run the simulation
    mushrooms(&current_grid, &next_grid, &ROWS, &COLUMNS, &TIME_STEPS, &current_row, &current_column, &current_time_step, &neighbor_row, &neighbor_column, &current_value, &prob, &yarn, &uniform);

    // end timing and print result
    end_time = c_get_wtime();
    total_time = end_time - start_time;
    #ifdef DEBUG
        printf("\nruntime: %f seconds\n", total_time);
    #else
        printf("%f", total_time);
    #endif

    // deallocate grids
    deallocateGrid(&current_grid, &ROWS, &current_row);
    deallocateGrid(&next_grid, &ROWS, &current_row);

    // return statement
    return 0;

}

/* getArguments() */
/* fetches and stores command line arguments for # of rows, columns, and time steps */
void getArguments(int argc, char *argv[], int * ROWS, int * COLUMNS, int * TIME_STEPS) {
    
    // declare + initialize variables
    int c;
    int rflag = 0;
    int cflag = 0;
    int sflag = 0;

    // retrieve command line arguments
    while ((c = getopt (argc, argv, "r:c:s:")) != -1) {
        switch (c) {
            case 'r':
                rflag = 1;
                *ROWS = atoi(optarg);
                break;
            
            case 'c':
                cflag = 1;
                *COLUMNS = atoi(optarg);
                break;

            case 's':
                sflag = 1;
                *TIME_STEPS = atoi(optarg);
                break;
            
            case '?':
                if (optopt == 'r') {
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else if (optopt == 'c') {
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else if (optopt == 's') {
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else if (isprint (optopt)) {
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                } else {
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                    exit(EXIT_FAILURE);
                }
        }
    }

    // check command line arguments
    if (rflag == 0) {
        fprintf(stderr, "Usage: %s -r number of rows\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (*ROWS < 1) {
        fprintf(stderr, "Usage: %s -r number of rows must be a positive nonzero integer\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (cflag == 0) {
        fprintf(stderr, "Usage: %s -c number of columns\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (*COLUMNS < 1) {
        fprintf(stderr, "Usage: %s -c number of columns must be a positive nonzero integer\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (sflag == 0) {
        fprintf(stderr, "Usage: %s -s number of time steps\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (*TIME_STEPS < 1) {
        fprintf(stderr, "Usage: %s -s number of time steps must be a positive nonzero integer\n", argv[0]);
        exit(EXIT_FAILURE);
    }
}

/* allocateGrid() */
/* allocates enough space for the input grid to store the rows and columns for the problem */
void allocateGrid(int ***grid, int * ROWS, int * COLUMNS, int * current_row) {
    *grid = new int*[(*ROWS) + 2];  // create pointer array
    for ((*current_row) = 0; (*current_row) <= (*ROWS) + 1; (*current_row)++) {  // at each element in the pointer array...
        (*grid)[*current_row] = new int[(*COLUMNS) + 2];  // ...create another pointer array
    }
}

/* initializeGrid() */
/* initializes the grid with empty spaces and spore spaces to begin the simulation */
void initializeGrid(int ***grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column, double * prob, trng::yarn2 * yarn, trng::uniform01_dist<> * uniform) {
    for ((*current_row) = 1; (*current_row) <= (*ROWS); (*current_row)++) {  // for each row in the grid...
        for ((*current_column) = 1; (*current_column) <= (*COLUMNS); (*current_column)++) {  // for each cell in that row...
            (*prob) = (*uniform)(*yarn);  // get random double between 0 and 1
            if ((*prob) <= probSpore) {  // if prob is less than or equal to probSpore...
                (*grid)[*current_row][*current_column] = SPORE;  // ...then cell starts as SPORE
            } else {  // otherwise...
                (*grid)[*current_row][*current_column] = EMPTY;  // ...cell starts as EMPTY
            }
        }
    }
}

/* mushrooms() */
/* simulates the growth of mushroom networks into fairy rings */
void mushrooms(int ***current_grid, int ***next_grid, int * ROWS, int * COLUMNS, int * TIME_STEPS, int * current_row, int * current_column, int * current_time_step, int * neighbor_row, int * neighbor_column, int * current_value, double * prob, trng::yarn2 * yarn, trng::uniform01_dist<> * uniform) {
    for((*current_time_step) = 0; (*current_time_step) <= (*TIME_STEPS); (*current_time_step)++) {  // for each time step...

        // set up ghost rows
        for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 1; (*current_column)++) {

            // set first row of grid to be the ghost of the second-to-last row
            (*current_grid)[0][*current_column] = (*current_grid)[(*ROWS)][*current_column];

            // set last row of grid to be the ghost of the second row
            (*current_grid)[(*ROWS) + 1][*current_column] = (*current_grid)[1][*current_column];
        }

        // set up ghost columns
        for ((*current_row) = 0; (*current_row) <= (*ROWS) + 1; (*current_row)++) {

            // set left-most column to be the ghost of the second-farthest-right column
            (*current_grid)[*current_row][0] = (*current_grid)[*current_row][*COLUMNS];

            // set right-most column to be the ghost of the second-farthest-left column
            (*current_grid)[*current_row][(*COLUMNS) + 1] = (*current_grid)[*current_row][1];
        }

        // DEBUG: display current grid
        #ifdef DEBUG
            #ifdef COLOR
                setlocale(LC_ALL, "");
                printf("\ntime step %d:\n", (*current_time_step));
                print_colorful_grid(current_grid, ROWS, COLUMNS, current_row, current_column, current_value);
            #else
                printf("\ntime step %d:\n", (*current_time_step));
                print_number_grid(current_grid, ROWS, COLUMNS, current_row, current_column);
            #endif
        #endif

        // determine grid at next time step
        for ((*current_row) = 1; (*current_row) <= (*ROWS); (*current_row)++) {  // for each row in the grid...
            for ((*current_column) = 1; (*current_column) <= (*COLUMNS); (*current_column)++) {  // for each cell in that row...

                (*current_value) = (*current_grid)[*current_row][*current_column];

                switch(*current_value) {
                
                    // if current cell is EMPTY...
                    case 0:
                        if (check_neighbors(current_grid, current_row, current_column, neighbor_row, neighbor_column) == 0) {  // if cell has no YOUNG neighbors...
                            (*next_grid)[*current_row][*current_column] == EMPTY;  // ...cell stays EMPTY in the next time step
                        } else {  // otherwise...
                            (*prob) = (*uniform)(*yarn);  // get random double between 0 and 1
                            if ((*prob) <= probSpread) {  // if prob is less than or equal to probSpread...
                                (*next_grid)[*current_row][*current_column] = YOUNG;  // ...cell becomes YOUNG in the next time step
                            } else {  // otherwise...
                                (*next_grid)[*current_row][*current_column] = EMPTY;  // ...cell stays EMPTY in the next time step
                            }
                        }
                        break;
                    
                    // if current cell is SPORE...
                    case 1:
                        (*prob) = (*uniform)(*yarn);  // get random double between 0 and 1
                        if ((*prob) <= probSporeToYoung) {  // if prob is less than or equal to probSporeToYoung...
                            (*next_grid)[*current_row][*current_column] = YOUNG;  // ...cell becomes YOUNG in the next time step
                        } else {  // otherwise...
                            (*next_grid)[*current_row][*current_column] = SPORE;  // ...cell stays SPORE in the next time step
                        }
                        break;
                    
                    // if current cell is YOUNG...
                    case 2:
                        (*next_grid)[*current_row][*current_column] = MATURING;  // ...cell becomes MATURING in the next time step
                        break;
                    
                    // if current cell is MATURING...
                    case 3:
                        (*prob) = (*uniform)(*yarn);  // get random double between 0 and 1
                        if ((*prob) <= probMaturingToMushrooms) {  // if prob is less than or equal to probMaturingToMushrooms...
                            (*next_grid)[*current_row][*current_column] = MUSHROOMS;  // ...cell becomes MUSHROOMS in the next time step
                        } else {  // otherwise...
                            (*next_grid)[*current_row][*current_column] = OLDER;  // ...cell becomes OLDER in the next time step
                        }
                        break;
                    
                    // if current cell is MUSHROOMS...
                    case 4:
                        (*next_grid)[*current_row][*current_column] = DECAYING;  // ...cell becomes DECAYING in the next time step
                        break;
                    
                    // if current cell is OLDER...
                    case 5:
                        (*next_grid)[*current_row][*current_column] = DECAYING;  // ...cell becomes DECAYING in the next time step
                        break;
                    
                    // if current cell is DECAYING...
                    case 6:
                        (*next_grid)[*current_row][*current_column] = DEAD;  // ...cell becomes DEAD in the next time step
                        break;
                    
                    // if current cell is DEAD...
                    case 7:
                        (*next_grid)[*current_row][*current_column] = DEADER;  // ...cell becomes DEADER in the next time step
                        break;
                    
                    // if current cell is DEADER...
                    case 8:
                        (*next_grid)[*current_row][*current_column] = DEPLETED;  // ...cell becomes DEPLETED in the next time step
                        break;
                    
                    // if current cell is DEPLETED...
                    case 9:
                        (*prob) = (*uniform)(*yarn);  // get random double between 0 and 1
                        if ((*prob) <= probDepletedToSpore) {  // if prob is less than or equal to probDepletedToSpore...
                            (*next_grid)[*current_row][*current_column] = SPORE;  // ...cell becomes SPORE in the next time step
                        } else if ((*prob) <= probDepletedToEmpty) {  // if prob is less than or equal to probDepletedToEmpty...
                            (*next_grid)[*current_row][*current_column] = EMPTY;  // ...cell becomes EMPTY in the next time step
                        } else {  // otherwise...
                            (*next_grid)[*current_row][*current_column] = DEPLETED;  // ...cell stays DEPLETED in the next time step
                        }
                        break;

                    // if current cell is INERT... (not currently used)
                    case 10:
                            // note: there is the potential to initialize the grid with some cells starting out as inert
                                // representing spots where fungi cannot grow (rocks etc.) but this has not been implemented
                        (*next_grid)[*current_row][*current_column] = INERT;  // ...cell stays INERT in the next time step
                        break;
                }
            }
        }
        
        // copy next_grid onto current_grid
        copyGrid(current_grid, next_grid, ROWS, COLUMNS, current_row, current_column);

        // loop simulation for the next time step
    }
}

/* copyGrid() */
/* copies the contents of one grid into another grid of the same size */
void copyGrid(int ***current_grid, int ***next_grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column) {
    for ((*current_row) = 1; (*current_row) <= (*ROWS); (*current_row)++) {  // for each row in the grid (except the ghost rows)...
        for ((*current_column) = 1; (*current_column) <= (*COLUMNS); (*current_column)++) {  // for each cell in that row...
            (*current_grid)[*current_row][*current_column] = (*next_grid)[*current_row][*current_column];  // ...store next_grid value in the same spot in current_grid
        }
    }
}

/* check_neighbors() */
/* checks the neighbors of a cell in the grid; returns 1 if at least one neighbor is YOUNG, otherwise returns 0 */
int check_neighbors(int ***current_grid, int * current_row, int * current_column, int * neighbor_row, int * neighbor_column) {
    for ((*neighbor_row) = (*current_row) - 1; (*neighbor_row) <= (*current_row) + 1; (*neighbor_row)++) {  // for each row in the 3x3 sub-grid...
        for ((*neighbor_column) = (*current_column) - 1; (*neighbor_column) <= (*current_column) + 1; (*neighbor_column)++) {  // for each cell in that row...
            if ( ((*neighbor_row) != (*current_row)) || ((*neighbor_column) != (*current_column)) ) {  // if that cell is a neighbor to the current cell...
                if ((*current_grid)[*neighbor_row][*neighbor_column] == YOUNG) {  // ... and if that neighbor is YOUNG...
                    return 1;  // return 1
                }
            }
        }
    }
    return 0;  // if none of the neighbors are YOUNG, return 0
}

/* deallocateGrid() */
/* deallocates the memory for the input grid */
void deallocateGrid(int ***grid, int * ROWS, int * current_row) {
    for ((*current_row) = 0; (*current_row) <= (*ROWS) + 1; (*current_row)++) {
        delete [] (*grid)[*current_row];
    }
    delete [] (*grid);
}

/* print_number_grid() */
/* prints the values in the input grid as numbers */
void print_number_grid(int ***grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column) {
    for ((*current_row) = 0; (*current_row) <= (*ROWS) + 1; (*current_row)++) {  // for each row in the grid...

        // if current_row is the second row, add a row of dashes (to separate the ghost row)
        if ((*current_row) == 1) {
            for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 1; (*current_column)++) {
                printf("--");
            }
            // new line
            printf("\n");
        }

        for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 1; (*current_column)++) {  // for each cell in that row...
            
            // if current column is the second-from-the-left column, add a column of dashes (to separate the ghost column)
            if ((*current_column) == 1) { printf("| "); }

            // print value of current cell
            printf("%d ", (*grid)[*current_row][*current_column]);

            // if current column is the second-from-the-right columns, add a column of dashes (to separate the ghost column)
            if ((*current_column) == (*COLUMNS)) { printf("| "); }
        }

        // new line
        printf("\n");
        
        // if current row is the second-to-last row, add a row of dashes (to separate the ghost row)
        if ((*current_row) == (*ROWS)) {
            for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 1; (*current_column)++) {
                printf("--");
            }
            // new line
            printf("\n");
        }
    }
    // new line
    printf("\n");
}

/* print_colorful_grid() */
/* prints the values in the input grid as color-coded blocks */
void print_colorful_grid(int ***grid, int * ROWS, int * COLUMNS, int * current_row, int * current_column, int * current_value) {

    // print color key
    printf("\nKEY:\n-----------------------------------------\n");
    printf("|\tEMPTY\t\t|");
        black();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    printf("|\n|\tSPORE\t\t|");
        red();
        printf("\t%lc\t", (wint_t)9547);
        reset_color();
    printf("|\n|\tYOUNG\t\t|");
        red();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    printf("|\n|\tMATURING\t|");
        green();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    printf("|\n|\tMUSHROOMS\t|");
        brown();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    printf("|\n|\tOLDER\t\t|");
        brown();
        printf("\t%lc\t", (wint_t)9619);
        reset_color();
    printf("|\n|\tDECAYING\t|");
        purple();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    printf("|\n|\tDEAD\t\t|");
        grey();
        printf("\t%lc\t", (wint_t)9619);
        reset_color();
    printf("|\n|\tDEADER\t\t|");
        grey();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    printf("|\n|\tDEPLETED\t|");
        black();
        printf("\t%lc\t", (wint_t)9608);
        reset_color();
    // uncomment if using inert
    // printf("|\n|\tINERT\t\t|");
    //     black();
    //     printf("\t%lc\t", (wint_t)9608);
    //     reset_color();
    printf("|\n-----------------------------------------\n\n");


    for ((*current_row) = 0; (*current_row) <= (*ROWS) + 1; (*current_row)++) {  // for each row in the grid...

        // if current_row is the second row, add a row of dashes (to separate the ghost row)
        if ((*current_row) == 1) {
            for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 6; (*current_column)++) {
                printf("-");
            }
            // new line
            printf("\n");
        }

        for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 1; (*current_column)++) {  // for each cell in that row...
            
            // if current column is the second-from-the-left column, add a column of dashes (to separate the ghost column)
            if ((*current_column) == 1) { printf(" | "); }

            // get current cell's value
            (*current_value) = (*grid)[*current_row][*current_column];

            // print current cell's value as color symbol
            switch(*current_value) {
                
                // EMPTY
                case 0:
                    black();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
                
                // SPORE
                case 1:
                    red();
                    printf("%lc", (wint_t)9547);
                    reset_color();
                    break;
                
                // YOUNG
                case 2:
                    red();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
                
                // MATURING
                case 3:
                    green();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
                
                // MUSHROOMS
                case 4:
                    brown();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
                
                // OLDER
                case 5:
                    brown();
                    printf("%lc", (wint_t)9619);
                    reset_color();
                    break;
                
                // DECAYING
                case 6:
                    purple();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
                
                // DEAD
                case 7:
                    grey();
                    printf("%lc", (wint_t)9619);
                    reset_color();
                    break;
                
                // DEADER
                case 8:
                    grey();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
                
                // DEPLETED
                case 9:
                    black();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;

                // INERT (not currently used)
                case 10:
                    black();
                    printf("%lc", (wint_t)9608);
                    reset_color();
                    break;
            }

            // if current column is the second-from-the-right columns, add a column of dashes (to separate the ghost column)
            if ((*current_column) == (*COLUMNS)) { printf(" | "); }
        }

        // new line
        printf("\n");
        
        // if current row is the second-to-last row, add a row of dashes (to separate the ghost row)
        if ((*current_row) == (*ROWS)) {
            for ((*current_column) = 0; (*current_column) <= (*COLUMNS) + 6; (*current_column)++) {
                printf("-");
            }
            // new line
            printf("\n");
        }
    }
    // new line
    printf("\n");
}

/* reset_color() */
/* resets the text color for printf statements */
void reset_color() {
    printf("\033[0m");
}

/* black() */
/* sets the text color for printf statements to black */
void black() {
    printf("\033[0;30m");
}

/* red() */
/* sets the text color for printf statements to red */
void red() {
    printf("\033[0;31m");
}

/* green() */
/* sets the text color for printf statements to green */
void green() {
    printf("\033[0;32m");
}

/* brown() */
/* sets the text color for printf statements to brown */
void brown() {
    printf("\033[0;33m");
}

/* grey() */
/* sets the text color for printf statements to grey */
void grey() {
    printf("\033[1;34m");
}

/* purple() */
/* sets the text color for printf statements to purple */
void purple() {
    printf("\033[1;35m");
}

// end of file