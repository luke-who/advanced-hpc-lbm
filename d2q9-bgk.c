/*
================================================================
**************** HPC cw by Luke Zhang (az16408) ****************
================================================================

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
#include <mpi.h>

#define NSPEEDS         9
#define FINALSTATEFILE  "final_state.dat"
#define AVVELSFILE      "av_vels.dat"
// #define DEBUG
#define MASTER 0 /* master rank */

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
  int tot_cells;       /* total number of cells that have obstacles==0*/
  /* mpi parameters */
  int kk;                /* index for looping over ranks */  
  int size;              /* number of processes in the communicator */
  int rank;              /* the rank of this process */
  int left;              /* the rank of the process to the left */
  int right;             /* the rank of the process to the right */
  int tag;               /* scope for adding extra information to a message */
  MPI_Status status;     /* struct used by MPI_Recv */
  int local_nrows;       /* number of rows apportioned to this rank */
  int local_ncols;       /* number of columns apportioned to this rank(excluding the extra two halo columns) */
  int start_col;         /* determine the index of the starting column in obstacles for this rank */
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
int initialise(const char* paramfile, const char* obstaclefile, t_param* params, 
               t_speed* restrict cells_ptr, t_speed* restrict tmp_cells_ptr, t_speed* restrict collated_cells,
               int** obstacles_ptr, float** av_vels_ptr, float** sendbuf, float** recvbuf,
               float** send_blockbuf, float** recv_blockbuf);

int cal_tot_cells(t_param params, int* obstacles);
/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
int halo_exchange(t_param params, t_speed* restrict cells, float* sendbuf, float* recvbuf);
float timestep(const t_param params, t_speed* restrict cells, t_speed* restrict tmp_cells, int* obstacles);
int accelerate_flow(const t_param params, t_speed* cells, int* obstacles);
float propa_rebd_collsn_av(const t_param params, t_speed* restrict cells, t_speed* restrict tmp_cells, int* obstacles);
int collate_cells(const t_param params, t_speed* restrict cells, t_speed* restrict collated_cells, float* send_blockbuf, float* recv_blockbuf);
int write_values(const t_param params, t_speed* collated_cells, int* obstacles, float* av_vels);

/* finalise, including freeing up allocated memory */
int finalise(const t_param* params, t_speed* restrict cells_ptr, t_speed* restrict tmp_cells_ptr,
             t_speed* restrict collated_cells, int** obstacles_ptr, float** av_vels_ptr, float** sendbuf, 
             float** recvbuf, float** send_blockbuf, float** recv_blockbuf);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed* cells);

/* compute average velocity */
float av_velocity(const t_param params, t_speed* collated_cells, int* obstacles);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed* collated_cells, int* obstacles);

/* calculate number of columns for each rank */
int calc_ncols_from_rank(int rank, int size, int tot_colmns);

/* calculate the starting column in collated cells for each rank */
int calc_start_columns_from_rank(int rank, int local_ncols, int size, int tot_colmns);

/* utility functions */
void die(const char* message, const int line, const char* file);
void usage(const char* exe);


/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* argv[])
{ 
  MPI_Init(&argc, &argv);

  char*    paramfile = NULL;    /* name of the input parameter file */
  char*    obstaclefile = NULL; /* name of a the input obstacle file */
  t_param  params;              /* struct to hold parameter values */
  t_speed  cells;               /* grid containing fluid densities */
  t_speed  tmp_cells;           /* scratch space */
  t_speed  collated_cells;      /* grid containing fluid densities from all ranks */
  int*     obstacles = NULL;    /* grid indicating which cells are blocked */
  float*   av_vels   = NULL;     /* a record of the av. velocity computed for each timestep */
  struct timeval timstr;                                                             /* structure to hold elapsed time */
  double tot_tic, tot_toc, init_tic, init_toc, comp_tic, comp_toc, col_tic, col_toc; /* floating point numbers to calculate elapsed wallclock time */

  float* sendbuf = NULL;       /* buffer to hold halo values to send */
  float* recvbuf = NULL;       /* buffer to hold received halo values */

  float* send_blockbuf = NULL;      /* buffer to send block values for collating */
  float* recv_blockbuf = NULL;      /* buffer to receive block values for collating */

  /* 
  ** MPI_Init returns once it has started up processes
  ** Get size of cohort and rank for this process
  */
  MPI_Comm_size(MPI_COMM_WORLD, &params.size);
  MPI_Comm_rank(MPI_COMM_WORLD, &params.rank);

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
  initialise(paramfile, obstaclefile, &params, &cells, &tmp_cells, &collated_cells, &obstacles, &av_vels, &sendbuf, &recvbuf, &send_blockbuf, &recv_blockbuf);
  /* calculate totoal number of cells that have obstacles==0, note only the MASTER rank needs this value to compute av_vel[tt] at each timestep */
  params.tot_cells = cal_tot_cells(params, obstacles);

  /* Init time stops here, compute time starts*/
  gettimeofday(&timstr, NULL);
  init_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  comp_tic=init_toc;

  for (int tt = 0; tt < params.maxIters; tt++)
  { 
    halo_exchange(params, &cells, sendbuf, recvbuf);
    float tot_u = timestep(params, &cells, &tmp_cells, obstacles);
    float tot_u_MASTER = 0.f;
    MPI_Reduce(&tot_u, &tot_u_MASTER, 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    if (params.rank == MASTER){
      av_vels[tt] = tot_u_MASTER / (float)params.tot_cells;
    }
    /*pointer swap here*/
    t_speed tmp_tmp_cells = cells;
    cells = tmp_cells;
    tmp_cells = tmp_tmp_cells;
    
#ifdef DEBUG
    float total = total_density(params, &cells);
    float total_MASTER = 0.f;
    MPI_Reduce(&total, &total_MASTER, 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    if (params.rank == MASTER){
      printf("==timestep: %d==\n", tt);
      printf("av velocity: %.12E\n", av_vels[tt]);
      // printf("tot density: %.12E\n", total_density(params, &cells));
      printf("tot density: %.12E\n", total_MASTER);
    }
#endif
  }
  // if (params.rank == MASTER){
  //   printf("av_vels[params.maxIters-1]: %f\n",av_vels[params.maxIters-1]);
  // }
  /* Compute time stops here, collate time starts*/
  gettimeofday(&timstr, NULL);
  comp_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  col_tic=comp_toc;

  // Collate data from all ranks to MASTER rank here 
  //TODO: use MPI_Gather() here? Be aware local_ncols might be different at each rank
  collate_cells(params, &cells, &collated_cells, send_blockbuf, recv_blockbuf);

  /* Total/collate time stops here.*/
  gettimeofday(&timstr, NULL);
  col_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  tot_toc = col_toc;
  if (params.rank == MASTER){
    /* write final values on MASTER rank and free memory on all ranks*/
    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, &collated_cells, obstacles));
    printf("Elapsed Init time:\t\t\t%.6lf (s)\n",    init_toc - init_tic);
    printf("Elapsed Compute time:\t\t\t%.6lf (s)\n", comp_toc - comp_tic);
    printf("Elapsed Collate time:\t\t\t%.6lf (s)\n", col_toc  - col_tic);
    printf("Elapsed Total time:\t\t\t%.6lf (s)\n",   tot_toc  - tot_tic);
    write_values(params, &collated_cells, obstacles, av_vels);
  }
  finalise(&params, &cells, &tmp_cells, &collated_cells, &obstacles, &av_vels, &sendbuf, &recvbuf, &send_blockbuf, &recv_blockbuf);

  MPI_Finalize();
  return EXIT_SUCCESS;
}

int cal_tot_cells(t_param params, int* obstacles){
  params.tot_cells=0;

  __assume(params.nx%8==0);
  __assume(params.ny%8==0);
  for (int jj = 0; jj < params.ny; jj++)
  { 
    #pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    { 
      /* calculate params.tot_cells after initialisation */
      params.tot_cells += (!obstacles[jj*params.nx + ii] ? 1 : 0);
    }
  }
  // printf("params.tot_cells:%d\n",params.tot_cells); //this should be 15876 for 128x128 grid

  return params.tot_cells;
}

float timestep(const t_param params, t_speed* restrict cells, t_speed* restrict tmp_cells, int* obstacles)
{
  accelerate_flow(params, cells, obstacles);

  return propa_rebd_collsn_av(params, cells, tmp_cells, obstacles);
}

int accelerate_flow(const t_param params, t_speed* cells, int* obstacles)
{
  /* compute weighting factors */
  const float w1 = params.density * params.accel / 9.f;
  const float w2 = params.density * params.accel / 36.f;

  /* modify the 2nd row of the grid */
  const int jj = params.ny - 2;

  __assume_aligned(cells->speeds_0, 64);
  __assume_aligned(cells->speeds_1, 64);
  __assume_aligned(cells->speeds_2, 64);
  __assume_aligned(cells->speeds_3, 64);
  __assume_aligned(cells->speeds_4, 64);
  __assume_aligned(cells->speeds_5, 64);
  __assume_aligned(cells->speeds_6, 64);
  __assume_aligned(cells->speeds_7, 64);
  __assume_aligned(cells->speeds_8, 64);

  __assume_aligned(obstacles, 64);

  __assume(params.nx%8==0);

  #pragma omp simd
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

float propa_rebd_collsn_av(const t_param params, t_speed* restrict cells, t_speed* restrict tmp_cells, int* obstacles)
{
  float* restrict cells_speeds_0 = cells->speeds_0;
  float* restrict cells_speeds_1 = cells->speeds_1;
  float* restrict cells_speeds_2 = cells->speeds_2;
  float* restrict cells_speeds_3 = cells->speeds_3;
  float* restrict cells_speeds_4 = cells->speeds_4;
  float* restrict cells_speeds_5 = cells->speeds_5;
  float* restrict cells_speeds_6 = cells->speeds_6;
  float* restrict cells_speeds_7 = cells->speeds_7;
  float* restrict cells_speeds_8 = cells->speeds_8;

  float* restrict tmp_cells_speeds_0 = tmp_cells->speeds_0;
  float* restrict tmp_cells_speeds_1 = tmp_cells->speeds_1;
  float* restrict tmp_cells_speeds_2 = tmp_cells->speeds_2;
  float* restrict tmp_cells_speeds_3 = tmp_cells->speeds_3;
  float* restrict tmp_cells_speeds_4 = tmp_cells->speeds_4;
  float* restrict tmp_cells_speeds_5 = tmp_cells->speeds_5;
  float* restrict tmp_cells_speeds_6 = tmp_cells->speeds_6;
  float* restrict tmp_cells_speeds_7 = tmp_cells->speeds_7;
  float* restrict tmp_cells_speeds_8 = tmp_cells->speeds_8;

  float tot_u = 0.f;          /* accumulated magnitudes of velocity for each cell */

  const float c_sq = 1.f / 3.f; /* square of speed of sound */
  const float w0 = 4.f / 9.f;  /* weighting factor */
  const float w1 = 1.f / 9.f;  /* weighting factor */
  const float w2 = 1.f / 36.f; /* weighting factor */

  __assume_aligned(cells_speeds_0, 64);
  __assume_aligned(cells_speeds_1, 64);
  __assume_aligned(cells_speeds_2, 64);
  __assume_aligned(cells_speeds_3, 64);
  __assume_aligned(cells_speeds_4, 64);
  __assume_aligned(cells_speeds_5, 64);
  __assume_aligned(cells_speeds_6, 64);
  __assume_aligned(cells_speeds_7, 64);
  __assume_aligned(cells_speeds_8, 64);

  __assume_aligned(tmp_cells_speeds_0, 64);
  __assume_aligned(tmp_cells_speeds_1, 64);
  __assume_aligned(tmp_cells_speeds_2, 64);
  __assume_aligned(tmp_cells_speeds_3, 64);
  __assume_aligned(tmp_cells_speeds_4, 64);
  __assume_aligned(tmp_cells_speeds_5, 64);
  __assume_aligned(tmp_cells_speeds_6, 64);
  __assume_aligned(tmp_cells_speeds_7, 64);
  __assume_aligned(tmp_cells_speeds_8, 64);

  __assume_aligned(obstacles, 64);

  __assume(params.nx%2==0);
  __assume(params.nx%4==0);
  __assume(params.nx%8==0);
  __assume(params.nx%16==0);
  __assume(params.nx%32==0);

  __assume(params.ny%2==0);
  __assume(params.ny%4==0);
  __assume(params.ny%8==0);
  __assume(params.ny%16==0);
  __assume(params.ny%32==0);

  /* loop over _all_ cells */
  // #pragma omp parallel for reduction(+:tot_u) //schedule(static)
  for (int jj = 0; jj < params.local_nrows; jj++)
  { 
    // printf("Thread ID: %d \t Number of threads:%d\n",omp_get_thread_num(),omp_get_num_threads());
    #pragma omp simd
    for (int ii = 1; ii < params.local_ncols+1; ii++)
    { 
      /******************* propagate & rebound *********************/
      /* determine indices of axis-direction neighbours
      ** respecting periodic boundary conditions (wrap around) */
      const int y_n = (jj + 1) % params.local_nrows;                        // still wrap around to the north
      const int x_e = (ii + 1);                                             // doesn't wrap around to the east anymore thanks to the halo regions
      const int y_s = (jj == 0) ? (jj + params.local_nrows - 1) : (jj - 1); // still wrap around to the south 
      const int x_w = (ii - 1);                                             // dosen't wrap around to the west anymore thanks to the halo regions
      /* propagate densities from neighbouring cells, following
      ** appropriate directions of travel and writing into
      ** scratch space grid */

      /*Fuse propagate,collision&rebound. Instead of reading from cells and writing to tmp_cells(propagate), then reading from tmp_cells and 
      writing to cells(rebound&collision), we can do this with one read from cells and one write to tmp_cells, followed by a pointer swap 
      between cells and tmp_cells afterwards*/

      const float cells0 = cells_speeds_0[ii + jj*(params.local_ncols+2)]; /* central cell, no movement */
      const float cells1 = cells_speeds_1[x_w + jj*(params.local_ncols+2)]; /* east */
      const float cells2 = cells_speeds_2[ii + y_s*(params.local_ncols+2)]; /* north */
      const float cells3 = cells_speeds_3[x_e + jj*(params.local_ncols+2)]; /* west */
      const float cells4 = cells_speeds_4[ii + y_n*(params.local_ncols+2)]; /* south */
      const float cells5 = cells_speeds_5[x_w + y_s*(params.local_ncols+2)]; /* north-east */
      const float cells6 = cells_speeds_6[x_e + y_s*(params.local_ncols+2)]; /* north-west */
      const float cells7 = cells_speeds_7[x_e + y_n*(params.local_ncols+2)]; /* south-west */
      const float cells8 = cells_speeds_8[x_w + y_n*(params.local_ncols+2)]; /* south-east */

      /* compute local density total */
      const float local_density   = cells0 
                                  + cells1 
                                  + cells2
                                  + cells3
                                  + cells4
                                  + cells5
                                  + cells6
                                  + cells7
                                  + cells8;

      /* compute x velocity component */
      const float u_x = (cells1
                       + cells5
                       + cells8
                      - (cells3
                       + cells6
                       + cells7))
                      / local_density;
      /* compute y velocity component */
      const float u_y = (cells2
                       + cells5
                       + cells6
                      - (cells4
                       + cells7
                       + cells8))
                      / local_density;

      /* velocity squared */
      const float u_sq = u_x * u_x + u_y * u_y;

      /* directional velocity components */
      // float u[NSPEEDS];
      const float u_1 =   u_x;        /* east */
      const float u_2 =         u_y;  /* north */
      const float u_3 = - u_x;        /* west */
      const float u_4 =       - u_y;  /* south */
      const float u_5 =   u_x + u_y;  /* north-east */
      const float u_6 = - u_x + u_y;  /* north-west */
      const float u_7 = - u_x - u_y;  /* south-west */
      const float u_8 =   u_x - u_y;  /* south-east */

      /* equilibrium densities */
      // float d_equ[NSPEEDS];
      /* zero velocity density: weight w0 */
      const float d_equ_0 = 
      w0 * local_density
                * (1.f - u_sq * 1.5f);
      /* axis speeds: weight w1 */
      const float d_equ_1 = w1 * local_density * (1.f + u_1 * 3.f
                                      + (u_1 * u_1) * 4.5f
                                      - u_sq * 1.5f);
      const float d_equ_2 = w1 * local_density * (1.f + u_2 * 3.f
                                      + (u_2 * u_2) * 4.5f
                                      - u_sq * 1.5f);
      const float d_equ_3 = w1 * local_density * (1.f + u_3 * 3.f
                                      + (u_3 * u_3) * 4.5f
                                      - u_sq * 1.5f);
      const float d_equ_4 = w1 * local_density * (1.f + u_4 * 3.f
                                      + (u_4 * u_4) * 4.5f
                                      - u_sq * 1.5f);
      /* diagonal speeds: weight w2 */
      const float d_equ_5 = w2 * local_density * (1.f + u_5 * 3.f
                                      + (u_5 * u_5) * 4.5f
                                      - u_sq * 1.5f);
      const float d_equ_6 = w2 * local_density * (1.f + u_6 * 3.f
                                      + (u_6 * u_6) * 4.5f
                                      - u_sq * 1.5f);
      const float d_equ_7 = w2 * local_density * (1.f + u_7 * 3.f
                                      + (u_7 * u_7) * 4.5f
                                      - u_sq * 1.5f);
      const float d_equ_8 = w2 * local_density * (1.f + u_8 * 3.f
                                      + (u_8 * u_8) * 4.5f
                                      - u_sq * 1.5f);

      /******************* propagate & collision/rebound *********************/
      tmp_cells_speeds_0[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells0 : cells0 * (1.f-params.omega) + params.omega*d_equ_0;
      tmp_cells_speeds_1[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells3 : cells1 * (1.f-params.omega) + params.omega*d_equ_1;
      tmp_cells_speeds_2[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells4 : cells2 * (1.f-params.omega) + params.omega*d_equ_2;
      tmp_cells_speeds_3[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells1 : cells3 * (1.f-params.omega) + params.omega*d_equ_3;
      tmp_cells_speeds_4[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells2 : cells4 * (1.f-params.omega) + params.omega*d_equ_4;
      tmp_cells_speeds_5[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells7 : cells5 * (1.f-params.omega) + params.omega*d_equ_5;
      tmp_cells_speeds_6[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells8 : cells6 * (1.f-params.omega) + params.omega*d_equ_6;
      tmp_cells_speeds_7[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells5 : cells7 * (1.f-params.omega) + params.omega*d_equ_7;
      tmp_cells_speeds_8[ii + jj*(params.local_ncols+2)] = (obstacles[(jj*params.nx+params.start_col) + ii-1]) ? cells6 : cells8 * (1.f-params.omega) + params.omega*d_equ_8;

      /* accumulate the norm of x- and y- velocity components */
      tot_u += (!obstacles[(jj*params.nx+params.start_col) + ii-1]) ? sqrtf((u_x * u_x) + (u_y * u_y)) : 0;

    }
  }

  return tot_u;
}

float av_velocity(const t_param params, t_speed* collated_cells, int* obstacles)
{
  float tot_u = 0.f;          /* accumulated magnitudes of velocity for each cell */
  
  __assume(params.nx%8==0);
  __assume(params.ny%8==0);

  /* loop over all non-blocked cells */
    // #pragma omp parallel for reduction(+:tot_u) //num_threads(28)
    for (int jj = 0; jj < params.ny; jj++)
    { 
      __assume_aligned(collated_cells->speeds_0, 64);
      __assume_aligned(collated_cells->speeds_1, 64);
      __assume_aligned(collated_cells->speeds_2, 64);
      __assume_aligned(collated_cells->speeds_3, 64);
      __assume_aligned(collated_cells->speeds_4, 64);
      __assume_aligned(collated_cells->speeds_5, 64);
      __assume_aligned(collated_cells->speeds_6, 64);
      __assume_aligned(collated_cells->speeds_7, 64);
      __assume_aligned(collated_cells->speeds_8, 64);

      __assume_aligned(obstacles, 64);
      
      #pragma omp simd
      for (int ii = 0; ii < params.nx; ii++)
      {
        /* local density total */
        const float local_density = collated_cells->speeds_0[ii + jj*params.nx]
                                  + collated_cells->speeds_1[ii + jj*params.nx]
                                  + collated_cells->speeds_2[ii + jj*params.nx]
                                  + collated_cells->speeds_3[ii + jj*params.nx]
                                  + collated_cells->speeds_4[ii + jj*params.nx]
                                  + collated_cells->speeds_5[ii + jj*params.nx]
                                  + collated_cells->speeds_6[ii + jj*params.nx]
                                  + collated_cells->speeds_7[ii + jj*params.nx]
                                  + collated_cells->speeds_8[ii + jj*params.nx];

        /* x-component of velocity */
        const float u_x = (collated_cells->speeds_1[ii + jj*params.nx]
                         + collated_cells->speeds_5[ii + jj*params.nx]
                         + collated_cells->speeds_8[ii + jj*params.nx]
                        - (collated_cells->speeds_3[ii + jj*params.nx]
                         + collated_cells->speeds_6[ii + jj*params.nx]
                        +  collated_cells->speeds_7[ii + jj*params.nx]))
                        / local_density;
        /* compute y velocity component */
        const float u_y = (collated_cells->speeds_2[ii + jj*params.nx]
                         + collated_cells->speeds_5[ii + jj*params.nx]
                         + collated_cells->speeds_6[ii + jj*params.nx]
                        - (collated_cells->speeds_4[ii + jj*params.nx]
                         + collated_cells->speeds_7[ii + jj*params.nx]
                         + collated_cells->speeds_8[ii + jj*params.nx]))
                        / local_density;
        /* accumulate the norm of x- and y- velocity components */
        tot_u += (!obstacles[jj*params.nx + ii]) ? sqrtf((u_x * u_x) + (u_y * u_y)) : 0;
        
      }
    }
  return tot_u / (float)params.tot_cells;
}

int initialise(const char* paramfile, const char* obstaclefile, t_param* params, 
               t_speed* restrict cells_ptr, t_speed* restrict tmp_cells_ptr, t_speed* restrict collated_cells,
               int** obstacles_ptr, float** av_vels_ptr, float** sendbuf, float** recvbuf,
               float** send_blockbuf, float** recv_blockbuf)
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
  ** determine process ranks to the left and right of rank
  ** respecting periodic boundary conditions
  */
  params->left = (params->rank == MASTER) ? (params->rank + params->size - 1) : (params->rank - 1);
  params->right = (params->rank + 1) % params->size;

  /* 
  ** determine local grid size
  ** each rank gets all the rows, but a subset of the number of columns
  */
  params->local_nrows = params->ny;
  params->local_ncols = calc_ncols_from_rank(params->rank, params->size, params->nx);
  if (params->local_ncols < 1) {
    fprintf(stderr,"Error: too many processes:- local_ncols < 1\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  /* determine the index of the starting column in obstacles/collated cells for this rank */
  params->start_col = calc_start_columns_from_rank(params->rank, params->local_ncols, params->size, params->nx);

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
  //After changing to SOA we to change how we initialise cells too (need to allocate memory for each of the speed_0, speed_1...)
  cells_ptr->speeds_0 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_1 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_2 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_3 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_4 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_5 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_6 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_7 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  cells_ptr->speeds_8 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);


  // if (*cells_ptr == NULL) die("cannot allocate memory for cells", __LINE__, __FILE__);

  /* 'helper' grid, used as scratch space */
  tmp_cells_ptr->speeds_0 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_1 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_2 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_3 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_4 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_5 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_6 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_7 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);
  tmp_cells_ptr->speeds_8 = (float*)_mm_malloc(sizeof(float) * (params->ny * (params->local_ncols+2)),64);

  // if (*tmp_cells_ptr == NULL) die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);
  
  /* allocate space to store a record of the send and receive columns for this rank's halo regions
  ** at each timestep */
  *sendbuf = (float*)_mm_malloc(sizeof(float) * (params->ny * NSPEEDS),64);
  *recvbuf = (float*)_mm_malloc(sizeof(float) * (params->ny * NSPEEDS),64);


  /* allocate space to store a record of the send and receive blocks(rows&columns) of this rank, for collating cells at the end */
  if (params->rank != MASTER){
    *send_blockbuf = (float*)_mm_malloc(sizeof(float) * (params->ny * params->local_ncols * NSPEEDS),64);
  }else{
    *recv_blockbuf = (float*)_mm_malloc(sizeof(float) * (params->ny * params->local_ncols * NSPEEDS),64);
  }


  /* the map of obstacles */
  *obstacles_ptr = _mm_malloc(sizeof(int) * (params->ny * params->nx),64);

  if (*obstacles_ptr == NULL) die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  /* initialise densities */
  const float w0 = params->density * 4.f / 9.f;
  const float w1 = params->density      / 9.f;
  const float w2 = params->density      / 36.f;

  // __assume_aligned(cells_ptr, 64); 
  __assume_aligned(cells_ptr->speeds_0, 64);
  __assume_aligned(cells_ptr->speeds_1, 64);
  __assume_aligned(cells_ptr->speeds_2, 64);
  __assume_aligned(cells_ptr->speeds_3, 64);
  __assume_aligned(cells_ptr->speeds_4, 64);
  __assume_aligned(cells_ptr->speeds_5, 64);
  __assume_aligned(cells_ptr->speeds_6, 64);
  __assume_aligned(cells_ptr->speeds_7, 64);
  __assume_aligned(cells_ptr->speeds_8, 64);

  __assume_aligned(*obstacles_ptr, 64);

  // __assume((*params)->nx%128==0);
  // __assume((*params)->ny%128==0);
  /*alternatively*/
  __assume(params->nx%8==0);
  __assume(params->ny%8==0);

  // #pragma omp parallel for
  for (int jj = 0; jj < params->ny; jj++)
  { 
    #pragma omp simd
    for (int ii = 1; ii < params->local_ncols+1; ii++)
    {                                  /* main grid for this rank */
      /* centre */
      cells_ptr->speeds_0[ii + jj*(params->local_ncols+2)] = w0;
      /* axis directions */
      cells_ptr->speeds_1[ii + jj*(params->local_ncols+2)] = w1;
      cells_ptr->speeds_2[ii + jj*(params->local_ncols+2)] = w1;
      cells_ptr->speeds_3[ii + jj*(params->local_ncols+2)] = w1;
      cells_ptr->speeds_4[ii + jj*(params->local_ncols+2)] = w1;
      /* diagonals */
      cells_ptr->speeds_5[ii + jj*(params->local_ncols+2)] = w2;
      cells_ptr->speeds_6[ii + jj*(params->local_ncols+2)] = w2;
      cells_ptr->speeds_7[ii + jj*(params->local_ncols+2)] = w2;
      cells_ptr->speeds_8[ii + jj*(params->local_ncols+2)] = w2;
      
      /* Initialise halo regions */
      if (ii==1) {                      /* left halo column */
        /* centre */
        cells_ptr->speeds_0[ii-1 + jj*(params->local_ncols+2)] = w0;
        /* axis directions */
        cells_ptr->speeds_1[ii-1 + jj*(params->local_ncols+2)] = w1;
        cells_ptr->speeds_2[ii-1 + jj*(params->local_ncols+2)] = w1;
        cells_ptr->speeds_3[ii-1 + jj*(params->local_ncols+2)] = w1;
        cells_ptr->speeds_4[ii-1 + jj*(params->local_ncols+2)] = w1;
        /* diagonals */
        cells_ptr->speeds_5[ii-1 + jj*(params->local_ncols+2)] = w2;
        cells_ptr->speeds_6[ii-1 + jj*(params->local_ncols+2)] = w2;
        cells_ptr->speeds_7[ii-1 + jj*(params->local_ncols+2)] = w2;
        cells_ptr->speeds_8[ii-1 + jj*(params->local_ncols+2)] = w2;
      }else if (ii==params->local_ncols){       /* right halo column */
        /* centre */
        cells_ptr->speeds_0[ii+1 + jj*(params->local_ncols+2)] = w0;
        /* axis directions */
        cells_ptr->speeds_1[ii+1 + jj*(params->local_ncols+2)] = w1;
        cells_ptr->speeds_2[ii+1 + jj*(params->local_ncols+2)] = w1;
        cells_ptr->speeds_3[ii+1 + jj*(params->local_ncols+2)] = w1;
        cells_ptr->speeds_4[ii+1 + jj*(params->local_ncols+2)] = w1;
        /* diagonals */
        cells_ptr->speeds_5[ii+1 + jj*(params->local_ncols+2)] = w2;
        cells_ptr->speeds_6[ii+1 + jj*(params->local_ncols+2)] = w2;
        cells_ptr->speeds_7[ii+1 + jj*(params->local_ncols+2)] = w2;
        cells_ptr->speeds_8[ii+1 + jj*(params->local_ncols+2)] = w2;
      }
    }
  }
  
  for (int jj = 0; jj < params->ny; jj++)
  { 
    #pragma omp simd
    for (int ii = 0; ii < params->nx; ii++)
    {
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
  
  /* Allocate memory for collated_cells on MASTER rank */
  if (params->rank == MASTER) {
    
    collated_cells->speeds_0 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_1 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_2 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_3 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_4 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_5 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_6 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_7 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
    collated_cells->speeds_8 = (float*)_mm_malloc(sizeof(float) * (params->ny * params->nx),64);
  }

  /*
  ** allocate space to hold a record of the avarage velocities computed
  ** at each timestep
  */
  *av_vels_ptr = (float*)_mm_malloc(sizeof(float) * params->maxIters,64);
  // printf("params->tot_cells: %d, params->ny: %d, params->nx: %d\n",params->tot_cells,params->ny,params->nx);
  // printf("rank %d * sendbuf size: %d\n",params->rank,(params->ny * NSPEEDS));
  // printf("rank %d <-> rank %d <-> rank %d, size: %d, params->local_ncols: %d, params->start_col: %d\n",params->left,params->rank,params->right,params->size, params->local_ncols, params->start_col);
  return EXIT_SUCCESS;
}

int finalise(const t_param* params, t_speed* restrict cells_ptr, t_speed* restrict tmp_cells_ptr, 
             t_speed* restrict collated_cells, int** obstacles_ptr, float** av_vels_ptr, float** sendbuf, 
             float** recvbuf, float** send_blockbuf, float** recv_blockbuf)
{
  /*
  ** free up allocated memory
  */
  _mm_free(cells_ptr->speeds_0);
  _mm_free(cells_ptr->speeds_1);
  _mm_free(cells_ptr->speeds_2);
  _mm_free(cells_ptr->speeds_3);
  _mm_free(cells_ptr->speeds_4);
  _mm_free(cells_ptr->speeds_5);
  _mm_free(cells_ptr->speeds_6);
  _mm_free(cells_ptr->speeds_7);
  _mm_free(cells_ptr->speeds_8);

  _mm_free(tmp_cells_ptr->speeds_0);
  _mm_free(tmp_cells_ptr->speeds_1);
  _mm_free(tmp_cells_ptr->speeds_2);
  _mm_free(tmp_cells_ptr->speeds_3);
  _mm_free(tmp_cells_ptr->speeds_4);
  _mm_free(tmp_cells_ptr->speeds_5);
  _mm_free(tmp_cells_ptr->speeds_6);
  _mm_free(tmp_cells_ptr->speeds_7);
  _mm_free(tmp_cells_ptr->speeds_8);

  if (params->rank == MASTER){
    _mm_free(collated_cells->speeds_0);
    _mm_free(collated_cells->speeds_1);
    _mm_free(collated_cells->speeds_2);
    _mm_free(collated_cells->speeds_3);
    _mm_free(collated_cells->speeds_4);
    _mm_free(collated_cells->speeds_5);
    _mm_free(collated_cells->speeds_6);
    _mm_free(collated_cells->speeds_7);
    _mm_free(collated_cells->speeds_8);
  }

  _mm_free(*obstacles_ptr);
  *obstacles_ptr = NULL;

  _mm_free(*av_vels_ptr);
  *av_vels_ptr = NULL;

  _mm_free(*sendbuf);
  *sendbuf = NULL;
  _mm_free(*recvbuf);
  *recvbuf = NULL;

  _mm_free(*send_blockbuf);
  *send_blockbuf = NULL;
  _mm_free(*recv_blockbuf);
  *recv_blockbuf = NULL;

  return EXIT_SUCCESS;
}


float calc_reynolds(const t_param params, t_speed* collated_cells, int* obstacles)
{
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

  return av_velocity(params, collated_cells, obstacles) * params.reynolds_dim / viscosity;
}

float total_density(const t_param params, t_speed* cells)
{
  float total = 0.f;  /* accumulator */


  __assume(params.nx%16==0);
  __assume(params.ny%16==0);
  
    // omp_set_num_threads(28);
    // #pragma omp parallel for //collapse(2)
    for (int jj = 0; jj < params.ny; jj++)
    { 
      __assume_aligned(cells->speeds_0, 64);
      __assume_aligned(cells->speeds_1, 64);
      __assume_aligned(cells->speeds_2, 64);
      __assume_aligned(cells->speeds_3, 64);
      __assume_aligned(cells->speeds_4, 64);
      __assume_aligned(cells->speeds_5, 64);
      __assume_aligned(cells->speeds_6, 64);
      __assume_aligned(cells->speeds_7, 64);
      __assume_aligned(cells->speeds_8, 64);
      
      #pragma omp simd
      for (int ii = 1; ii < (params.local_ncols+1); ii++)
      {
        total +=  cells->speeds_0[ii + jj*(params.local_ncols+2)]
                + cells->speeds_1[ii + jj*(params.local_ncols+2)]
                + cells->speeds_2[ii + jj*(params.local_ncols+2)]
                + cells->speeds_3[ii + jj*(params.local_ncols+2)]
                + cells->speeds_4[ii + jj*(params.local_ncols+2)]
                + cells->speeds_5[ii + jj*(params.local_ncols+2)]
                + cells->speeds_6[ii + jj*(params.local_ncols+2)]
                + cells->speeds_7[ii + jj*(params.local_ncols+2)]
                + cells->speeds_8[ii + jj*(params.local_ncols+2)];
      }
    }

  //TODO: use MPI_Reduce() here?
  return total;
}

int write_values(const t_param params, t_speed* collated_cells, int* obstacles, float* av_vels)
{
  FILE* fp;                     /* file pointer */
  const float c_sq = 1.f / 3.f; /* sq. of speed of sound */
  // float local_density;         /* per grid cell sum of densities */
  float pressure;              /* fluid pressure in grid cell */
  float u_x;                   /* x-component of velocity in grid cell */
  float u_y;                   /* y-component of velocity in grid cell */
  float u;                     /* norm--root of summed squares--of u_x and u_y */

  fp = fopen(FINALSTATEFILE, "w");

  if (fp == NULL)
  {
    die("could not open file output file", __LINE__, __FILE__);
  }
  

  for (int jj = 0; jj < params.ny; jj++)
  { 
    for (int ii = 1; ii < (params.nx); ii++)
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
        const float local_density = collated_cells->speeds_0[ii + jj*params.nx]
                                  + collated_cells->speeds_1[ii + jj*params.nx]
                                  + collated_cells->speeds_2[ii + jj*params.nx]
                                  + collated_cells->speeds_3[ii + jj*params.nx]
                                  + collated_cells->speeds_4[ii + jj*params.nx]
                                  + collated_cells->speeds_5[ii + jj*params.nx]
                                  + collated_cells->speeds_6[ii + jj*params.nx]
                                  + collated_cells->speeds_7[ii + jj*params.nx]
                                  + collated_cells->speeds_8[ii + jj*params.nx];

        /* compute x velocity component */
        u_x = (collated_cells->speeds_1[ii + jj*params.nx]
             + collated_cells->speeds_5[ii + jj*params.nx]
             + collated_cells->speeds_8[ii + jj*params.nx]
            - (collated_cells->speeds_3[ii + jj*params.nx]
             + collated_cells->speeds_6[ii + jj*params.nx]
             + collated_cells->speeds_7[ii + jj*params.nx]))
              / local_density;
        /* compute y velocity component */
        u_y = (collated_cells->speeds_2[ii + jj*params.nx]
             + collated_cells->speeds_5[ii + jj*params.nx]
             + collated_cells->speeds_6[ii + jj*params.nx]
            - (collated_cells->speeds_4[ii + jj*params.nx]
             + collated_cells->speeds_7[ii + jj*params.nx]
             + collated_cells->speeds_8[ii + jj*params.nx]))
              / local_density;
        /* compute norm of velocity */
        u = sqrtf((u_x * u_x) + (u_y * u_y));
        /* compute pressure */
        pressure = local_density * c_sq;
      }

      /* write to file */
      fprintf(fp, "%d %d %.12E %.12E %.12E %.12E %d\n", ii, jj, u_x, u_y, u, pressure, obstacles[ii + params.nx * jj]);
    }
  }

  fclose(fp);

  fp = fopen(AVVELSFILE, "w");

  if (fp == NULL)
  {
    die("could not open file output file", __LINE__, __FILE__);
  }
  #pragma omp simd
  for (int ii = 0; ii < params.maxIters; ii++)
  {
    fprintf(fp, "%d:\t%.12E\n", ii, av_vels[ii]);
  }

  fclose(fp);

  return EXIT_SUCCESS;
}

int calc_ncols_from_rank(int rank, int size, int tot_colmns)
{
  int ncols;
  int remainder = 0;

  ncols = tot_colmns / size;       /* integer division */
  remainder = tot_colmns % size;
  if (remainder != 0) {  /* if there is a remainder */
    if (rank < remainder) {
      ncols+=1;           /* distribute the remainder 1 per rank */
    }
    // if (rank == size - 1)
    //   ncols += NCOLS % size;  /* add remainder to last rank , but this can cause severe load imbalance in the worst case (slower, more variable timing) */
  }
  
  return ncols;
}

int calc_start_columns_from_rank(int rank, int local_ncols, int size, int tot_colmns){
  /* determine the index of the starting column in obstacles/collated cells for this rank */
  // int start_col;
  // if (params.nx % params.size == 0){
  //   start_col = params.rank * params.local_ncols;
  // }else {
  //   if (params.local_ncols == (int)(params.nx/params.size) + 1){
  //     start_col = params.rank * params.local_ncols;
  //   }else {
  //     start_col = params.rank * params.local_ncols + params.nx%params.size;
  //   }
  // }
  int start_col = ( (tot_colmns % size == 0) || (local_ncols == (tot_colmns/size) + 1) ) ? 
                       (rank * local_ncols) : (rank * local_ncols + tot_colmns%size);
  return start_col;           
}

int halo_exchange(t_param params, t_speed* restrict cells,
                  float* sendbuf, float* recvbuf) {
  
  //##### packing from cells to sendbuf: pack the 2nd column of the cells to send buffer #####//
  for (int jj = 0; jj < params.local_nrows; jj++){
    sendbuf[0 + (jj*NSPEEDS)] = cells->speeds_0[jj*(params.local_ncols+2) + 1];
    sendbuf[1 + (jj*NSPEEDS)] = cells->speeds_1[jj*(params.local_ncols+2) + 1];
    sendbuf[2 + (jj*NSPEEDS)] = cells->speeds_2[jj*(params.local_ncols+2) + 1];
    sendbuf[3 + (jj*NSPEEDS)] = cells->speeds_3[jj*(params.local_ncols+2) + 1];
    sendbuf[4 + (jj*NSPEEDS)] = cells->speeds_4[jj*(params.local_ncols+2) + 1];
    sendbuf[5 + (jj*NSPEEDS)] = cells->speeds_5[jj*(params.local_ncols+2) + 1];
    sendbuf[6 + (jj*NSPEEDS)] = cells->speeds_6[jj*(params.local_ncols+2) + 1];
    sendbuf[7 + (jj*NSPEEDS)] = cells->speeds_7[jj*(params.local_ncols+2) + 1];
    sendbuf[8 + (jj*NSPEEDS)] = cells->speeds_8[jj*(params.local_ncols+2) + 1];
  }
  // printf("sendbuf size: %d\n",(int)(sizeof(sendbuf)/sizeof(sendbuf[0])));
  /* send to the left, receive from right <- [0] <- [1] <- [2] <- [3] <- ... <- [n] <- */
  MPI_Sendrecv(sendbuf, (params.ny*NSPEEDS), MPI_FLOAT, params.left, params.tag,
               recvbuf, (params.ny*NSPEEDS), MPI_FLOAT, params.right, params.tag, MPI_COMM_WORLD, &params.status);
  // printf("recvbuf size: %d\n",(int)(sizeof(recvbuf)/sizeof(recvbuf[0])));
  //##### unpacking from recvbuf to cells: unpack from recvbuf to the right most(halo) column of cells #####//
  for (int jj = 0; jj < params.local_nrows; jj++){
    cells->speeds_0[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[0 + (jj*NSPEEDS)];
    cells->speeds_1[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[1 + (jj*NSPEEDS)];
    cells->speeds_2[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[2 + (jj*NSPEEDS)];
    cells->speeds_3[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[3 + (jj*NSPEEDS)];
    cells->speeds_4[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[4 + (jj*NSPEEDS)];
    cells->speeds_5[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[5 + (jj*NSPEEDS)];
    cells->speeds_6[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[6 + (jj*NSPEEDS)];
    cells->speeds_7[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[7 + (jj*NSPEEDS)];
    cells->speeds_8[(jj+1)*(params.local_ncols+2) - 1] = recvbuf[8 + (jj*NSPEEDS)];  
  }

  //------------------------------------------------------------------------------------------//

  //##### packing from cells to sendbuf: pack the params.local_ncols+1 column of the cells to send buffer #####//
  for (int jj = 0; jj < params.local_nrows; jj++){
    sendbuf[0 + (jj*NSPEEDS)] = cells->speeds_0[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[1 + (jj*NSPEEDS)] = cells->speeds_1[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[2 + (jj*NSPEEDS)] = cells->speeds_2[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[3 + (jj*NSPEEDS)] = cells->speeds_3[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[4 + (jj*NSPEEDS)] = cells->speeds_4[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[5 + (jj*NSPEEDS)] = cells->speeds_5[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[6 + (jj*NSPEEDS)] = cells->speeds_6[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[7 + (jj*NSPEEDS)] = cells->speeds_7[jj*(params.local_ncols+2) + params.local_ncols];
    sendbuf[8 + (jj*NSPEEDS)] = cells->speeds_8[jj*(params.local_ncols+2) + params.local_ncols];
  }
  
  /* send to the right, receive from left -> [0] -> [1] -> [2] -> [3] -> ... -> [n] -> */
  MPI_Sendrecv(sendbuf, (params.ny*NSPEEDS), MPI_FLOAT, params.right, params.tag,
               recvbuf, (params.ny*NSPEEDS), MPI_FLOAT, params.left, params.tag, MPI_COMM_WORLD, &params.status);

  //##### unpacking from recvbuf to cells: unpack from recvbuf to the left most(halo) column of cells #####//
  for (int jj = 0; jj < params.local_nrows; jj++){
    cells->speeds_0[jj*(params.local_ncols+2)] = recvbuf[0 + (jj*NSPEEDS)];
    cells->speeds_1[jj*(params.local_ncols+2)] = recvbuf[1 + (jj*NSPEEDS)];
    cells->speeds_2[jj*(params.local_ncols+2)] = recvbuf[2 + (jj*NSPEEDS)];
    cells->speeds_3[jj*(params.local_ncols+2)] = recvbuf[3 + (jj*NSPEEDS)];
    cells->speeds_4[jj*(params.local_ncols+2)] = recvbuf[4 + (jj*NSPEEDS)];
    cells->speeds_5[jj*(params.local_ncols+2)] = recvbuf[5 + (jj*NSPEEDS)];
    cells->speeds_6[jj*(params.local_ncols+2)] = recvbuf[6 + (jj*NSPEEDS)];
    cells->speeds_7[jj*(params.local_ncols+2)] = recvbuf[7 + (jj*NSPEEDS)];
    cells->speeds_8[jj*(params.local_ncols+2)] = recvbuf[8 + (jj*NSPEEDS)];  
  } 

  return EXIT_SUCCESS;
}

int collate_cells(t_param params, t_speed* restrict cells, t_speed* restrict collated_cells, float* send_blockbuf, float* restrict recv_blockbuf){
  //collate cells from all ranks into collated_cells on MASTER rank (MPI_Gatherv?)
  if (params.rank != MASTER){
    //##### packing from cells to send_blockbuf: pack the all columns(excluding halo) of the cells to send block buffer #####//
    for (int jj = 0; jj < params.local_nrows; jj++){
      for (int ii = 0; ii < params.local_ncols; ii++){
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 0)] = cells->speeds_0[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 1)] = cells->speeds_1[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 2)] = cells->speeds_2[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 3)] = cells->speeds_3[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 4)] = cells->speeds_4[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 5)] = cells->speeds_5[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 6)] = cells->speeds_6[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 7)] = cells->speeds_7[(ii+1) + jj*(params.local_ncols+2)];
        send_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 8)] = cells->speeds_8[(ii+1) + jj*(params.local_ncols+2)];        
      }
    }
    MPI_Send(send_blockbuf, (params.local_nrows * params.local_ncols * NSPEEDS), MPI_FLOAT, MASTER, params.tag, MPI_COMM_WORLD);

  }else {
    /* collate the cells value from MASTER first */
    for (int jj = 0; jj < params.local_nrows; jj++){
      for (int ii = 0; ii < params.local_ncols; ii++){
        collated_cells->speeds_0[ii + jj*params.nx] = cells->speeds_0[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_1[ii + jj*params.nx] = cells->speeds_1[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_2[ii + jj*params.nx] = cells->speeds_2[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_3[ii + jj*params.nx] = cells->speeds_3[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_4[ii + jj*params.nx] = cells->speeds_4[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_5[ii + jj*params.nx] = cells->speeds_5[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_6[ii + jj*params.nx] = cells->speeds_6[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_7[ii + jj*params.nx] = cells->speeds_7[(ii+1) + jj*(params.local_ncols+2)];
        collated_cells->speeds_8[ii + jj*params.nx] = cells->speeds_8[(ii+1) + jj*(params.local_ncols+2)];   
      }
    }
    /* then collate the cells value from other ranks */
    for (int source=1; source<params.size; source++) {
      /* recieving cells values from all other ranks other than MASTER.. */
      int local_ncols = calc_ncols_from_rank(source, params.size, params.nx);
      int start_col = calc_start_columns_from_rank(source, local_ncols, params.size, params.nx);

      MPI_Recv(recv_blockbuf, (params.local_nrows * local_ncols * NSPEEDS), MPI_FLOAT, source, params.tag, MPI_COMM_WORLD, &params.status);
      for (int jj = 0; jj < params.local_nrows; jj++){
        for (int ii = 0; ii < local_ncols; ii++){
          collated_cells->speeds_0[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 0)];
          collated_cells->speeds_1[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 1)];
          collated_cells->speeds_2[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 2)];
          collated_cells->speeds_3[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 3)];
          collated_cells->speeds_4[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 4)];
          collated_cells->speeds_5[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 5)];
          collated_cells->speeds_6[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 6)];
          collated_cells->speeds_7[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 7)];
          collated_cells->speeds_8[ii + start_col + jj*params.nx] = recv_blockbuf[ii + jj*(params.local_ncols) + (params.local_nrows * params.local_ncols * 8)];        
        }
      } 
    }    
  }
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
