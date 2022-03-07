/*
** Code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the Bhatnagar-Gross-Krook collision step.
**
** The 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** A 2D grid:
**
**           cols
**       --- --- ---
**      | D | E | F |
** rows  --- --- ---
**      | A | B | C |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1D array:
**
**  --- --- --- --- --- ---
** | A | B | C | D | E | F |
**  --- --- --- --- --- ---
**
** Grid indicies are:
**
**          ny
**          ^       cols(ii)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(jj) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** Note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   ./d2q9-bgk input.params obstacles.dat
**
** Be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>

#define NSPEEDS         9
#define FINALSTATEFILE  "final_state.dat"
#define AVVELSFILE      "av_vels.dat"

/* struct to hold the parameter values */
typedef struct
{
  int    nx;            /* no. of cells in x-direction */
  int    ny;            /* no. of cells in y-direction */
  int    maxIters;      /* no. of iterations */
  int    reynolds_dim;  /* dimension for Reynolds number */
  float density;       /* density per link */
  float accel;         /* density redistribution */
  float omega;         /* relaxation parameter */
} t_param;


/* struct to hold the 'speed' values */
// Array of Structure
// typedef struct
// {
//   float speeds[NSPEEDS];
// } t_speed;

// Structure of Array/Pointer
typedef struct
{
  float* restrict speeds_0;
  float* restrict speeds_1;
  float* restrict speeds_2;
  float* restrict speeds_3;
  float* restrict speeds_4;
  float* restrict speeds_5;
  float* restrict speeds_6;
  float* restrict speeds_7;
  float* restrict speeds_8;
} t_speed;

/*
** function prototypes
*/

/* load params, allocate memory, load obstacles & initialise fluid particle densities */
int initialise(const char* paramfile, const char* obstaclefile,
               t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
               int** obstacles_ptr, float** av_vels_ptr);

/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
int timestep(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
int accelerate_flow(const t_param params, t_speed* cells, int* obstacles);
int propagate(const t_param params, t_speed* cells, t_speed* tmp_cells);
int rebound(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
int collision(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
int propa_rebd_collsn(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
int write_values(const t_param params, t_speed* cells, int* obstacles, float* av_vels);

/* finalise, including freeing up allocated memory */
int finalise(const t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
             int** obstacles_ptr, float** av_vels_ptr);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed* cells);

/* compute average velocity */
float av_velocity(const t_param params, t_speed* cells, int* obstacles);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed* cells, int* obstacles);

/* utility functions */
void die(const char* message, const int line, const char* file);
void usage(const char* exe);

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* argv[])
{
  char*    paramfile = NULL;    /* name of the input parameter file */
  char*    obstaclefile = NULL; /* name of a the input obstacle file */
  t_param  params;              /* struct to hold parameter values */
  t_speed* cells     = NULL;    /* grid containing fluid densities */
  t_speed* tmp_cells = NULL;    /* scratch space */
  int*     obstacles = NULL;    /* grid indicating which cells are blocked */
  float* av_vels   = NULL;     /* a record of the av. velocity computed for each timestep */
  struct timeval timstr;                                                             /* structure to hold elapsed time */
  double tot_tic, tot_toc, init_tic, init_toc, comp_tic, comp_toc, col_tic, col_toc; /* floating point numbers to calculate elapsed wallclock time */

  /* parse the command line */
  if (argc != 3)
  {
    usage(argv[0]);
  }
  else
  {
    paramfile = argv[1];
    obstaclefile = argv[2];
  }

  /* Total/init time starts here: initialise our data structures and load values from file */
  gettimeofday(&timstr, NULL);
  tot_tic = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  init_tic=tot_tic;
  initialise(paramfile, obstaclefile, &params, &cells, &tmp_cells, &obstacles, &av_vels);

  /* Init time stops here, compute time starts*/
  gettimeofday(&timstr, NULL);
  init_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  comp_tic=init_toc;


  for (int tt = 0; tt < params.maxIters; tt++)
  {
    timestep(params, cells, tmp_cells, obstacles);
    /*pointer swap here*/
    t_speed *tmp_tmp_cells = cells;
    cells = tmp_cells;
    tmp_cells = tmp_tmp_cells;
    av_vels[tt] = av_velocity(params, cells, obstacles);
#ifdef DEBUG
    printf("==timestep: %d==\n", tt);
    printf("av velocity: %.12E\n", av_vels[tt]);
    printf("tot density: %.12E\n", total_density(params, cells));
#endif
  }
  
  /* Compute time stops here, collate time starts*/
  gettimeofday(&timstr, NULL);
  comp_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  col_tic=comp_toc;

  // Collate data from ranks here 

  /* Total/collate time stops here.*/
  gettimeofday(&timstr, NULL);
  col_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  tot_toc = col_toc;
  
  /* write final values and free memory */
  printf("==done==\n");
  printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, cells, obstacles));
  printf("Elapsed Init time:\t\t\t%.6lf (s)\n",    init_toc - init_tic);
  printf("Elapsed Compute time:\t\t\t%.6lf (s)\n", comp_toc - comp_tic);
  printf("Elapsed Collate time:\t\t\t%.6lf (s)\n", col_toc  - col_tic);
  printf("Elapsed Total time:\t\t\t%.6lf (s)\n",   tot_toc  - tot_tic);
  write_values(params, cells, obstacles, av_vels);
  finalise(&params, &cells, &tmp_cells, &obstacles, &av_vels);

  return EXIT_SUCCESS;
}
 
int timestep(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles)
{
  accelerate_flow(params, cells, obstacles);
  // propagate(params, cells, tmp_cells);
  // rebound(params, cells, tmp_cells, obstacles);
  // collision(params, cells, tmp_cells, obstacles);
  propa_rebd_collsn(params, cells, tmp_cells, obstacles);
  return EXIT_SUCCESS;
}

int accelerate_flow(const t_param params, t_speed* cells, int* obstacles)
{
  /* compute weighting factors */
  float w1 = params.density * params.accel / 9.f;
  float w2 = params.density * params.accel / 36.f;

  /* modify the 2nd row of the grid */
  int jj = params.ny - 2;


  // #pragma omp parallel for //num_threads(28)
  for (int ii = 0; ii < params.nx; ii++)
  { 
    /* if the cell is not occupied and
    ** we don't send a negative density */
    if (!obstacles[ii + jj*params.nx]
        && (cells->speeds_3[ii + jj*params.nx] - w1) > 0.f
        && (cells->speeds_6[ii + jj*params.nx] - w2) > 0.f
        && (cells->speeds_7[ii + jj*params.nx] - w2) > 0.f)
    {
      /* increase 'east-side' densities */
      cells->speeds_1[ii + jj*params.nx] += w1;
      cells->speeds_5[ii + jj*params.nx] += w2;
      cells->speeds_8[ii + jj*params.nx] += w2;
      /* decrease 'west-side' densities */
      cells->speeds_3[ii + jj*params.nx] -= w1;
      cells->speeds_6[ii + jj*params.nx] -= w2;
      cells->speeds_7[ii + jj*params.nx] -= w2;
    }
  }

  return EXIT_SUCCESS;
}

int propa_rebd_collsn(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles)
{
  const float c_sq = 1.f / 3.f; /* square of speed of sound */
  const float w0 = 4.f / 9.f;  /* weighting factor */
  const float w1 = 1.f / 9.f;  /* weighting factor */
  const float w2 = 1.f / 36.f; /* weighting factor */
  /* loop over _all_ cells */
  // printf("Thread ID: %d \t Number of threads:%d\n",omp_get_thread_num(),omp_get_num_threads()); //1 thread before omp parallel
  // #pragma omp parallel //num_threads(28)
  // {
    // printf("Thread ID: %d \t Number of threads:%d\n",omp_get_thread_num(),omp_get_num_threads()); //n threads before omp parallel
    // #pragma omp for collapse(2)
  for (int jj = 0; jj < params.ny; jj++)
  { 
    #pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    { 
      /******************* propagate & rebound *********************/
      /* determine indices of axis-direction neighbours
      ** respecting periodic boundary conditions (wrap around) */
      int y_n = (jj + 1) % params.ny;
      int x_e = (ii + 1) % params.nx;
      int y_s = (jj == 0) ? (jj + params.ny - 1) : (jj - 1);
      int x_w = (ii == 0) ? (ii + params.nx - 1) : (ii - 1);
      /* propagate densities from neighbouring cells, following
      ** appropriate directions of travel and writing into
      ** scratch space grid */

      /*Fuse propagate,collision&rebound. Instead of reading from cells and writing to tmp_cells(propagate), then reading from tmp_cells and 
      writing to cells(rebound&collision), we can do this with one read from cells and one write to tmp_cells, followed by a pointer swap 
      between cells and tmp_cells afterwards*/
      const float cells0 = cells->speeds_0[ii + jj*params.nx]; /* central cell, no movement */
      const float cells1 = cells->speeds_1[x_w + jj*params.nx]; /* east */
      const float cells2 = cells->speeds_2[ii + y_s*params.nx]; /* north */
      const float cells3 = cells->speeds_3[x_e + jj*params.nx]; /* west */
      const float cells4 = cells->speeds_4[ii + y_n*params.nx]; /* south */
      const float cells5 = cells->speeds_5[x_w + y_s*params.nx]; /* north-east */
      const float cells6 = cells->speeds_6[x_e + y_s*params.nx]; /* north-west */
      const float cells7 = cells->speeds_7[x_e + y_n*params.nx]; /* south-west */
      const float cells8 = cells->speeds_8[x_w + y_n*params.nx]; /* south-east */


      /******************* propagate & collision *********************/
      /* if the cell contains an obstacle */
      if (obstacles[jj*params.nx + ii])
      {
        /* called after propagate, so taking values from scratch space
        ** mirroring, and writing into main grid */

        tmp_cells->speeds_0[ii + jj*params.nx] = cells0;
        tmp_cells->speeds_1[ii + jj*params.nx] = cells3;
        tmp_cells->speeds_2[ii + jj*params.nx] = cells4;
        tmp_cells->speeds_3[ii + jj*params.nx] = cells1;
        tmp_cells->speeds_4[ii + jj*params.nx] = cells2;
        tmp_cells->speeds_5[ii + jj*params.nx] = cells7;
        tmp_cells->speeds_6[ii + jj*params.nx] = cells8;
        tmp_cells->speeds_7[ii + jj*params.nx] = cells5;
        tmp_cells->speeds_8[ii + jj*params.nx] = cells6;

      }
      /* don't consider occupied cells */
      if (!obstacles[ii + jj*params.nx])
      {
        /* compute local density total */
        float local_density = 0.f;

        local_density = cells0 
                      + cells1 
                      + cells2
                      + cells3
                      + cells4
                      + cells5
                      + cells6
                      + cells7
                      + cells8;

        /* compute x velocity component */
        float u_x = (cells1
                      + cells5
                      + cells8
                      - (cells3
                        + cells6
                        + cells7))
                    / local_density;
        /* compute y velocity component */
        float u_y = (cells2
                      + cells5
                      + cells6
                      - (cells4
                        + cells7
                        + cells8))
                    / local_density;

        /* velocity squared */
        float u_sq = u_x * u_x + u_y * u_y;

        /* directional velocity components */
        float u[NSPEEDS];
        u[1] =   u_x;        /* east */
        u[2] =         u_y;  /* north */
        u[3] = - u_x;        /* west */
        u[4] =       - u_y;  /* south */
        u[5] =   u_x + u_y;  /* north-east */
        u[6] = - u_x + u_y;  /* north-west */
        u[7] = - u_x - u_y;  /* south-west */
        u[8] =   u_x - u_y;  /* south-east */

        /* equilibrium densities */
        float d_equ[NSPEEDS];
        /* zero velocity density: weight w0 */
        d_equ[0] = 
        w0 * local_density
                  * (1.f - u_sq * 1.5f);
        /* axis speeds: weight w1 */
        d_equ[1] = w1 * local_density * (1.f + u[1] * 3.f
                                        + (u[1] * u[1]) * 4.5f
                                        - u_sq * 1.5f);
        d_equ[2] = w1 * local_density * (1.f + u[2] * 3.f
                                        + (u[2] * u[2]) * 4.5f
                                        - u_sq * 1.5f);
        d_equ[3] = w1 * local_density * (1.f + u[3] * 3.f
                                        + (u[3] * u[3]) * 4.5f
                                        - u_sq * 1.5f);
        d_equ[4] = w1 * local_density * (1.f + u[4] * 3.f
                                        + (u[4] * u[4]) * 4.5f
                                        - u_sq * 1.5f);
        /* diagonal speeds: weight w2 */
        d_equ[5] = w2 * local_density * (1.f + u[5] * 3.f
                                        + (u[5] * u[5]) * 4.5f
                                        - u_sq * 1.5f);
        d_equ[6] = w2 * local_density * (1.f + u[6] * 3.f
                                        + (u[6] * u[6]) * 4.5f
                                        - u_sq * 1.5f);
        d_equ[7] = w2 * local_density * (1.f + u[7] * 3.f
                                        + (u[7] * u[7]) * 4.5f
                                        - u_sq * 1.5f);
        d_equ[8] = w2 * local_density * (1.f + u[8] * 3.f
                                        + (u[8] * u[8]) * 4.5f
                                        - u_sq * 1.5f);

        /* relaxation step */

          tmp_cells->speeds_0[ii + jj*params.nx] = cells0 *(1.f-params.omega) 
                                                  + params.omega*d_equ[0];
          tmp_cells->speeds_1[ii + jj*params.nx] = cells1 *(1.f-params.omega) 
                                                  + params.omega*d_equ[1];
          tmp_cells->speeds_2[ii + jj*params.nx] = cells2 *(1.f-params.omega) 
                                                  + params.omega*d_equ[2];
          tmp_cells->speeds_3[ii + jj*params.nx] = cells3 *(1.f-params.omega) 
                                                  + params.omega*d_equ[3];
          tmp_cells->speeds_4[ii + jj*params.nx] = cells4 *(1.f-params.omega) 
                                                  + params.omega*d_equ[4];
          tmp_cells->speeds_5[ii + jj*params.nx] = cells5 *(1.f-params.omega) 
                                                  + params.omega*d_equ[5];
          tmp_cells->speeds_6[ii + jj*params.nx] = cells6 *(1.f-params.omega) 
                                                  + params.omega*d_equ[6];
          tmp_cells->speeds_7[ii + jj*params.nx] = cells7 *(1.f-params.omega) 
                                                  + params.omega*d_equ[7];
          tmp_cells->speeds_8[ii + jj*params.nx] = cells8 *(1.f-params.omega) 
                                                  + params.omega*d_equ[8];
      }
    }
  }
  // }

  return EXIT_SUCCESS;
}


float av_velocity(const t_param params, t_speed* cells, int* obstacles)
{
  int    tot_cells = 0;  /* no. of cells used in calculation */
  float tot_u;          /* accumulated magnitudes of velocity for each cell */

  /* initialise */
  tot_u = 0.f;

  /* loop over all non-blocked cells */
    // #pragma omp parallel for collapse(2) //This doens't work correctly - to be fixed
    for (int jj = 0; jj < params.ny; jj++)
    {
      for (int ii = 0; ii < params.nx; ii++)
      {
        /* ignore occupied cells */
        if (!obstacles[ii + jj*params.nx])
        {
          /* local density total */
          float local_density = 0.f;

          // for (int kk = 0; kk < NSPEEDS; kk++)
          // {
          //   local_density += cells[ii + jj*params.nx].speeds[kk];
          // }
          local_density = cells->speeds_0[ii + jj*params.nx]
                        + cells->speeds_1[ii + jj*params.nx]
                        + cells->speeds_2[ii + jj*params.nx]
                        + cells->speeds_3[ii + jj*params.nx]
                        + cells->speeds_4[ii + jj*params.nx]
                        + cells->speeds_5[ii + jj*params.nx]
                        + cells->speeds_6[ii + jj*params.nx]
                        + cells->speeds_7[ii + jj*params.nx]
                        + cells->speeds_8[ii + jj*params.nx];

          /* x-component of velocity */
          float u_x = (cells->speeds_1[ii + jj*params.nx]
                        + cells->speeds_5[ii + jj*params.nx]
                        + cells->speeds_8[ii + jj*params.nx]
                        - (cells->speeds_3[ii + jj*params.nx]
                          + cells->speeds_6[ii + jj*params.nx]
                          + cells->speeds_7[ii + jj*params.nx]))
                      / local_density;
          /* compute y velocity component */
          float u_y = (cells->speeds_2[ii + jj*params.nx]
                        + cells->speeds_5[ii + jj*params.nx]
                        + cells->speeds_6[ii + jj*params.nx]
                        - (cells->speeds_4[ii + jj*params.nx]
                          + cells->speeds_7[ii + jj*params.nx]
                          + cells->speeds_8[ii + jj*params.nx]))
                      / local_density;
          /* accumulate the norm of x- and y- velocity components */
          tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
          /* increase counter of inspected cells */
          ++tot_cells;
        }
      }
    }

  return tot_u / (float)tot_cells;
}

int initialise(const char* paramfile, const char* obstaclefile,
               t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
               int** obstacles_ptr, float** av_vels_ptr)
{
  char   message[1024];  /* message buffer */
  FILE*   fp;            /* file pointer */
  int    xx, yy;         /* generic array indices */
  int    blocked;        /* indicates whether a cell is blocked by an obstacle */
  int    retval;         /* to hold return value for checking */

  /* open the parameter file */
  fp = fopen(paramfile, "r");

  if (fp == NULL)
  {
    sprintf(message, "could not open input parameter file: %s", paramfile);
    die(message, __LINE__, __FILE__);
  }

  /* read in the parameter values */
  retval = fscanf(fp, "%d\n", &(params->nx));

  if (retval != 1) die("could not read param file: nx", __LINE__, __FILE__);

  retval = fscanf(fp, "%d\n", &(params->ny));

  if (retval != 1) die("could not read param file: ny", __LINE__, __FILE__);

  retval = fscanf(fp, "%d\n", &(params->maxIters));

  if (retval != 1) die("could not read param file: maxIters", __LINE__, __FILE__);

  retval = fscanf(fp, "%d\n", &(params->reynolds_dim));

  if (retval != 1) die("could not read param file: reynolds_dim", __LINE__, __FILE__);

  retval = fscanf(fp, "%f\n", &(params->density));

  if (retval != 1) die("could not read param file: density", __LINE__, __FILE__);

  retval = fscanf(fp, "%f\n", &(params->accel));

  if (retval != 1) die("could not read param file: accel", __LINE__, __FILE__);

  retval = fscanf(fp, "%f\n", &(params->omega));

  if (retval != 1) die("could not read param file: omega", __LINE__, __FILE__);

  /* and close up the file */
  fclose(fp);

  /*
  ** Allocate memory.
  **
  ** Remember C is pass-by-value, so we need to
  ** pass pointers into the initialise function.
  **
  ** NB we are allocating a 1D array, so that the
  ** memory will be contiguous.  We still want to
  ** index this memory as if it were a (row major
  ** ordered) 2D array, however.  We will perform
  ** some arithmetic using the row and column
  ** coordinates, inside the square brackets, when
  ** we want to access elements of this array.
  **
  ** Note also that we are using a structure to
  ** hold an array of 'speeds'.  We will allocate
  ** a 1D array of these structs.
  */

  /* main grid */
  *cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->ny * params->nx));
  //After changing to SOA we to change how we initialise cells too (need to allocate memory for each of the speed_0, speed_1...)
  (*cells_ptr)->speeds_0 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_1 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_2 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_3 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_4 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_5 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_6 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_7 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*cells_ptr)->speeds_8 = (float*)malloc(sizeof(float) * (params->ny * params->nx));


  if (*cells_ptr == NULL) die("cannot allocate memory for cells", __LINE__, __FILE__);

  /* 'helper' grid, used as scratch space */
  *tmp_cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_0 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_1 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_2 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_3 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_4 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_5 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_6 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_7 = (float*)malloc(sizeof(float) * (params->ny * params->nx));
  (*tmp_cells_ptr)->speeds_8 = (float*)malloc(sizeof(float) * (params->ny * params->nx));

  if (*tmp_cells_ptr == NULL) die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);

  /* the map of obstacles */
  *obstacles_ptr = malloc(sizeof(int) * (params->ny * params->nx));

  if (*obstacles_ptr == NULL) die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  /* initialise densities */
  float w0 = params->density * 4.f / 9.f;
  float w1 = params->density      / 9.f;
  float w2 = params->density      / 36.f;

  // omp_set_num_threads(28);
  // #pragma omp parallel for collapse(2) 
  for (int jj = 0; jj < params->ny; jj++)
  {
    for (int ii = 0; ii < params->nx; ii++)
    {
      /* centre */
      (*cells_ptr)->speeds_0[ii + jj*params->nx] = w0;
      /* axis directions */
      (*cells_ptr)->speeds_1[ii + jj*params->nx] = w1;
      (*cells_ptr)->speeds_2[ii + jj*params->nx] = w1;
      (*cells_ptr)->speeds_3[ii + jj*params->nx] = w1;
      (*cells_ptr)->speeds_4[ii + jj*params->nx] = w1;
      /* diagonals */
      (*cells_ptr)->speeds_5[ii + jj*params->nx] = w2;
      (*cells_ptr)->speeds_6[ii + jj*params->nx] = w2;
      (*cells_ptr)->speeds_7[ii + jj*params->nx] = w2;
      (*cells_ptr)->speeds_8[ii + jj*params->nx] = w2;


  /* first set all cells in obstacle array to zero */
      (*obstacles_ptr)[ii + jj*params->nx] = 0;
    }
  }

  /* open the obstacle data file */
  fp = fopen(obstaclefile, "r");

  if (fp == NULL)
  {
    sprintf(message, "could not open input obstacles file: %s", obstaclefile);
    die(message, __LINE__, __FILE__);
  }

  /* read-in the blocked cells list */
  while ((retval = fscanf(fp, "%d %d %d\n", &xx, &yy, &blocked)) != EOF)
  {
    /* some checks */
    if (retval != 3) die("expected 3 values per line in obstacle file", __LINE__, __FILE__);

    if (xx < 0 || xx > params->nx - 1) die("obstacle x-coord out of range", __LINE__, __FILE__);

    if (yy < 0 || yy > params->ny - 1) die("obstacle y-coord out of range", __LINE__, __FILE__);

    if (blocked != 1) die("obstacle blocked value should be 1", __LINE__, __FILE__);

    /* assign to array */
    (*obstacles_ptr)[xx + yy*params->nx] = blocked;
  }

  /* and close the file */
  fclose(fp);

  /*
  ** allocate space to hold a record of the avarage velocities computed
  ** at each timestep
  */
  *av_vels_ptr = (float*)malloc(sizeof(float) * params->maxIters);

  return EXIT_SUCCESS;
}

int finalise(const t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
             int** obstacles_ptr, float** av_vels_ptr)
{
  /*
  ** free up allocated memory
  */
  free(*cells_ptr);
  *cells_ptr = NULL;

  free(*tmp_cells_ptr);
  *tmp_cells_ptr = NULL;

  free(*obstacles_ptr);
  *obstacles_ptr = NULL;

  free(*av_vels_ptr);
  *av_vels_ptr = NULL;

  return EXIT_SUCCESS;
}


float calc_reynolds(const t_param params, t_speed* cells, int* obstacles)
{
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

  return av_velocity(params, cells, obstacles) * params.reynolds_dim / viscosity;
}

float total_density(const t_param params, t_speed* cells)
{
  float total = 0.f;  /* accumulator */


    // #pragma omp parallel for collapse(2)
    for (int jj = 0; jj < params.ny; jj++)
    {
      for (int ii = 0; ii < params.nx; ii++)
      {
        // for (int kk = 0; kk < NSPEEDS; kk++)
        // {
        //   total += cells[ii + jj*params.nx].speeds[kk];
        // }
        total = cells->speeds_0[ii + jj*params.nx]
              + cells->speeds_1[ii + jj*params.nx]
              + cells->speeds_2[ii + jj*params.nx]
              + cells->speeds_3[ii + jj*params.nx]
              + cells->speeds_4[ii + jj*params.nx]
              + cells->speeds_5[ii + jj*params.nx]
              + cells->speeds_6[ii + jj*params.nx]
              + cells->speeds_7[ii + jj*params.nx]
              + cells->speeds_8[ii + jj*params.nx];
      }
    }

  return total;
}

int write_values(const t_param params, t_speed* cells, int* obstacles, float* av_vels)
{
  FILE* fp;                     /* file pointer */
  const float c_sq = 1.f / 3.f; /* sq. of speed of sound */
  float local_density;         /* per grid cell sum of densities */
  float pressure;              /* fluid pressure in grid cell */
  float u_x;                   /* x-component of velocity in grid cell */
  float u_y;                   /* y-component of velocity in grid cell */
  float u;                     /* norm--root of summed squares--of u_x and u_y */

  fp = fopen(FINALSTATEFILE, "w");

  if (fp == NULL)
  {
    die("could not open file output file", __LINE__, __FILE__);
  }
  // #pragma omp parallel for collapse(2) //this one doens't give the right result-to be fixed
  for (int jj = 0; jj < params.ny; jj++)
  {
    for (int ii = 0; ii < params.nx; ii++)
    {
      /* an occupied cell */
      if (obstacles[ii + jj*params.nx])
      {
        u_x = u_y = u = 0.f;
        pressure = params.density * c_sq;
      }
      /* no obstacle */
      else
      {
        local_density = 0.f;

        // for (int kk = 0; kk < NSPEEDS; kk++)
        // {
        //   local_density += cells[ii + jj*params.nx].speeds[kk];
        // }
        local_density = cells->speeds_0[ii + jj*params.nx]
                      + cells->speeds_1[ii + jj*params.nx]
                      + cells->speeds_2[ii + jj*params.nx]
                      + cells->speeds_3[ii + jj*params.nx]
                      + cells->speeds_4[ii + jj*params.nx]
                      + cells->speeds_5[ii + jj*params.nx]
                      + cells->speeds_6[ii + jj*params.nx]
                      + cells->speeds_7[ii + jj*params.nx]
                      + cells->speeds_8[ii + jj*params.nx];

        /* compute x velocity component */
        u_x = (cells->speeds_1[ii + jj*params.nx]
               + cells->speeds_5[ii + jj*params.nx]
               + cells->speeds_8[ii + jj*params.nx]
               - (cells->speeds_3[ii + jj*params.nx]
                  + cells->speeds_6[ii + jj*params.nx]
                  + cells->speeds_7[ii + jj*params.nx]))
              / local_density;
        /* compute y velocity component */
        u_y = (cells->speeds_2[ii + jj*params.nx]
               + cells->speeds_5[ii + jj*params.nx]
               + cells->speeds_6[ii + jj*params.nx]
               - (cells->speeds_4[ii + jj*params.nx]
                  + cells->speeds_7[ii + jj*params.nx]
                  + cells->speeds_8[ii + jj*params.nx]))
              / local_density;
        /* compute norm of velocity */
        u = sqrtf((u_x * u_x) + (u_y * u_y));
        /* compute pressure */
        pressure = local_density * c_sq;
      }

      /* write to file */
      fprintf(fp, "%d %d %.12E %.12E %.12E %.12E %d\n", ii, jj, u_x, u_y, u, pressure, obstacles[ii * params.nx + jj]);
    }
  }

  fclose(fp);

  fp = fopen(AVVELSFILE, "w");

  if (fp == NULL)
  {
    die("could not open file output file", __LINE__, __FILE__);
  }
  // #pragma omp parallel for //schedule(static) //this one doesn't work yet-to be fixed
  for (int ii = 0; ii < params.maxIters; ii++)
  {
    fprintf(fp, "%d:\t%.12E\n", ii, av_vels[ii]);
  }

  fclose(fp);

  return EXIT_SUCCESS;
}

void die(const char* message, const int line, const char* file)
{
  fprintf(stderr, "Error at line %d of file %s:\n", line, file);
  fprintf(stderr, "%s\n", message);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void usage(const char* exe)
{
  fprintf(stderr, "Usage: %s <paramfile> <obstaclefile>\n", exe);
  exit(EXIT_FAILURE);
}