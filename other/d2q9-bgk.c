/*
**======================================================
**
** Coursework done by            Marc Goulding (mg17752)
**
**======================================================
**
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
typedef struct
{
  float* speed0;
  float* speed1;
  float* speed2;
  float* speed3;
  float* speed4;
  float* speed5;
  float* speed6;
  float* speed7;
  float* speed8;
} t_speed;

/*
** function prototypes
*/
/* load params, allocate memory, load obstacles & initialise fluid particle densities */
int initialise(const char* restrict paramfile,
              const char* restrict obstaclefile,
               t_param* restrict params,
               int** restrict obstacles_ptr,
               float** restrict av_vels_ptr);
float timestep(const t_param params, t_speed* restrict cells, t_speed* restrict tmp_cells, int* restrict obstacles);
/* compute average velocity */
float av_velocity(const t_param params, const t_speed* restrict cells, const int* restrict obstacles);
int write_values(const t_param params, t_speed* restrict cells, int* restrict obstacles, float* restrict av_vels);
/* finalise, including freeing up allocated memory */
int finalise(const t_param* restrict params, t_speed* restrict cells_ptr, t_speed* restrict tmp_cells_ptr,
             int** restrict obstacles_ptr, float** restrict av_vels_ptr);
/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed* restrict cells, int* restrict obstacles);
/* utility functions */
void die(const char* restrict message, const int line, const char* restrict file);
void usage(const char* restrict exe);

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* restrict argv[])
{
  char*    paramfile = NULL;    /* name of the input parameter file */
  char*    obstaclefile = NULL; /* name of a the input obstacle file */
  t_param  params;              /* struct to hold parameter values */
  t_speed cells;    /* grid containing fluid densities */
  t_speed tmp_cells;    /* scratch space */
  int*     obstacles = NULL;    /* grid indicating which cells are blocked */
  float* av_vels   = NULL;     /* a record of the av. velocity computed for each timestep */
  struct timeval timstr;        /* structure to hold elapsed time */
  struct rusage ru;             /* structure to hold CPU time--system and user */
  double tic, toc;              /* floating point numbers to calculate elapsed wallclock time */
  double usrtim;                /* floating point number to record elapsed user CPU time */
  double systim;                /* floating point number to record elapsed system CPU time */

  /* parse the command line */
  if (argc != 3) usage(argv[0]); else
  {
    paramfile = argv[1];
    obstaclefile = argv[2];
  }

  /* initialise our data structures and load values from file */
  initialise(paramfile, obstaclefile, &params, &obstacles, &av_vels);
  cells    .speed0 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed1 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed2 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed3 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed4 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed5 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed6 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed7 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  cells    .speed8 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed0 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed1 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed2 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed3 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed4 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed5 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed6 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed7 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  tmp_cells.speed8 = (float*)_mm_malloc(sizeof(float) * params.ny * params.nx, 64);
  /* initialise densities */
  const float w0 = params.density * 4.f / 9.f;
  const float w1 = params.density       / 9.f;
  const float w2 = params.density      / 36.f;

  #pragma omp parallel for
  for (int jj = 0; jj < params.ny; jj++)
  {
    __assume_aligned(cells.speed0, 64);
    __assume_aligned(cells.speed1, 64);
    __assume_aligned(cells.speed2, 64);
    __assume_aligned(cells.speed3, 64);
    __assume_aligned(cells.speed4, 64);
    __assume_aligned(cells.speed5, 64);
    __assume_aligned(cells.speed6, 64);
    __assume_aligned(cells.speed7, 64);
    __assume_aligned(cells.speed8, 64);
    #pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    {
      /* centre */
      cells.speed0[ii + jj*params.nx] = w0;
      /* axis directions */
      cells.speed1[ii + jj*params.nx] = w1;
      cells.speed2[ii + jj*params.nx] = w1;
      cells.speed3[ii + jj*params.nx] = w1;
      cells.speed4[ii + jj*params.nx] = w1;
      /* diagonals */
      cells.speed5[ii + jj*params.nx] = w2;
      cells.speed6[ii + jj*params.nx] = w2;
      cells.speed7[ii + jj*params.nx] = w2;
      cells.speed8[ii + jj*params.nx] = w2;
    }
  }

  /* iterate for maxIters timesteps */
  gettimeofday(&timstr, NULL);
  tic = timstr.tv_sec + (timstr.tv_usec / 1000000.0);

  for (int tt = 0; tt < params.maxIters; tt+=2)
  {
  __assume(tt%2==0);

    av_vels[tt]   = timestep(params, &cells, &tmp_cells, obstacles);
    av_vels[tt+1] = timestep(params, &tmp_cells, &cells, obstacles);

  }

  gettimeofday(&timstr, NULL);
  toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  getrusage(RUSAGE_SELF, &ru);
  timstr = ru.ru_utime;
  usrtim = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  timstr = ru.ru_stime;
  systim = timstr.tv_sec + (timstr.tv_usec / 1000000.0);

  /* write final values and free memory */
  printf("==done==\n");
  printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, &cells, obstacles));
  printf("Elapsed time:\t\t\t%.6lf (s)\n", toc - tic);
  printf("Elapsed user CPU time:\t\t%.6lf (s)\n", usrtim);
  printf("Elapsed system CPU time:\t%.6lf (s)\n", systim);
  write_values(params, &cells, obstacles, av_vels);
  finalise(&params, &cells, &tmp_cells, &obstacles, &av_vels);

  return EXIT_SUCCESS;
}

float timestep(const t_param params, t_speed* restrict cells, t_speed* restrict tmp_cells, int* restrict obstacles)
{
  /* ==================== ACCELERATE FLOW ====================*/
  
  /* compute weighting factors */
  const float w_1 = params.density * params.accel / 9.f;
  const float w_2 = params.density * params.accel / 36.f;
  /* modify the 2nd row of the grid */
  int jj = params.ny - 2;

  __assume_aligned(cells    ->speed0, 64);
  __assume_aligned(cells    ->speed1, 64);
  __assume_aligned(cells    ->speed2, 64);
  __assume_aligned(cells    ->speed3, 64);
  __assume_aligned(cells    ->speed4, 64);
  __assume_aligned(cells    ->speed5, 64);
  __assume_aligned(cells    ->speed6, 64);
  __assume_aligned(cells    ->speed7, 64);
  __assume_aligned(cells    ->speed8, 64);
  __assume_aligned(tmp_cells->speed0, 64);
  __assume_aligned(tmp_cells->speed1, 64);
  __assume_aligned(tmp_cells->speed2, 64);
  __assume_aligned(tmp_cells->speed3, 64);
  __assume_aligned(tmp_cells->speed4, 64);
  __assume_aligned(tmp_cells->speed5, 64);
  __assume_aligned(tmp_cells->speed6, 64);
  __assume_aligned(tmp_cells->speed7, 64);
  __assume_aligned(tmp_cells->speed8, 64);
  __assume_aligned(obstacles        , 64);
  #pragma omp parallel for simd
  for (int ii = 0; ii < params.nx; ii++)
  {
    /* if the cell is not occupied and
    ** we don't send a negative density */
    if (!obstacles[ii + jj*params.nx]
        && (cells->speed3[ii + jj*params.nx] - w_1) > 0.f
        && (cells->speed6[ii + jj*params.nx] - w_2) > 0.f
        && (cells->speed7[ii + jj*params.nx] - w_2) > 0.f)
    {
      cells->speed1[ii + jj*params.nx] += w_1;
      cells->speed3[ii + jj*params.nx] -= w_1;
      cells->speed5[ii + jj*params.nx] += w_2;
      cells->speed6[ii + jj*params.nx] -= w_2;
      cells->speed7[ii + jj*params.nx] -= w_2;
      cells->speed8[ii + jj*params.nx] += w_2;
    }
  }

  /* ====================== LBM ====================== */
  const float w0 = 4.f / 9.f;  /* weighting factor */
  const float w1 = 1.f / 9.f;  /* weighting factor */
  const float w2 = 1.f / 36.f; /* weighting factor */

  int   tot_cells = 0;        /* no. of cells used in calculation */
  float tot_u = 0.f;          /* accumulated magnitudes of velocity for each cell */

  __assume_aligned(cells    ->speed0, 64);
  __assume_aligned(cells    ->speed1, 64);
  __assume_aligned(cells    ->speed2, 64);
  __assume_aligned(cells    ->speed3, 64);
  __assume_aligned(cells    ->speed4, 64);
  __assume_aligned(cells    ->speed5, 64);
  __assume_aligned(cells    ->speed6, 64);
  __assume_aligned(cells    ->speed7, 64);
  __assume_aligned(cells    ->speed8, 64);
  __assume_aligned(tmp_cells->speed0, 64);
  __assume_aligned(tmp_cells->speed1, 64);
  __assume_aligned(tmp_cells->speed2, 64);
  __assume_aligned(tmp_cells->speed3, 64);
  __assume_aligned(tmp_cells->speed4, 64);
  __assume_aligned(tmp_cells->speed5, 64);
  __assume_aligned(tmp_cells->speed6, 64);
  __assume_aligned(tmp_cells->speed7, 64);
  __assume_aligned(tmp_cells->speed8, 64);
  __assume_aligned(obstacles        , 64);

  const int step=512;
  #pragma omp parallel for reduction(+:tot_cells,tot_u)
 for (int k = 0; k < params.ny*params.nx; k+=step)
 {
    #pragma omp simd
    for (int ii = k; ii < k+step; ii++)
    {

      const int y_n = (ii>(params.ny-1)*params.nx) ? -(params.ny-1)*params.nx : (params.nx);
      const int x_e = ((ii+1)%params.nx==0) ? (1-params.nx) : (1);
      const int y_s = (ii<params.nx) ? (params.nx*(params.ny-1)) : (-params.nx);
      const int x_w = (ii%params.nx==0) ? (params.nx - 1) : (- 1);

      __declspec(align(64)) float sp[NSPEEDS];

      sp[0] = cells->speed0[ii];
      sp[1] = cells->speed1[ii+x_w];
      sp[2] = cells->speed2[ii  + y_s];
      sp[3] = cells->speed3[x_e + ii ];
      sp[4] = cells->speed4[ii  + y_n];
      sp[5] = cells->speed5[ii+ x_w + y_s];
      sp[6] = cells->speed6[ii+ x_e + y_s];
      sp[7] = cells->speed7[ii+ x_e + y_n];
      sp[8] = cells->speed8[ii+ x_w + y_n];

      /* compute local density total */
      const float local_density = sp[0]
                                + sp[3]
                                + sp[4]
                                + sp[1]
                                + sp[2]
                                + sp[7]
                                + sp[8]
                                + sp[5]
                                + sp[6];

      /* compute x velocity component */
      const float u_x =(sp[1]
                      + sp[5]
                      + sp[8]
                      -(sp[3]
                      + sp[6]
                      + sp[7]))
                      / local_density;
      /* compute y velocity component */
      const float u_y =(sp[2]
                      + sp[5]
                      + sp[6]
                      -(sp[4]
                      + sp[7]
                      + sp[8]))
                      / local_density;

      tot_u += (!obstacles[ii]) ? sqrtf((u_x * u_x) + (u_y * u_y)) : 0;
      tot_cells += (!obstacles[ii]) ? 1 : 0;

      // /* accumulate the norm of x- and y- velocity components */
      // tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
      // /* increase counter of inspected cells */
      // ++tot_cells;

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

      /* relaxation step */
      const float op = (u_x * u_x + u_y * u_y) * 1.5f;
      tmp_cells->speed0[ii] = (obstacles[ii]) ? sp[0] : sp[0] + params.omega * (w0 * local_density * (1.f                             - op) - sp[0]);
      tmp_cells->speed1[ii] = (obstacles[ii]) ? sp[3] : sp[1] + params.omega * (w1 * local_density * (1.f + u[1] * (3.0f + u[1]*4.5f) - op) - sp[1]);
      tmp_cells->speed2[ii] = (obstacles[ii]) ? sp[4] : sp[2] + params.omega * (w1 * local_density * (1.f + u[2] * (3.0f + u[2]*4.5f) - op) - sp[2]);
      tmp_cells->speed3[ii] = (obstacles[ii]) ? sp[1] : sp[3] + params.omega * (w1 * local_density * (1.f + u[3] * (3.0f + u[3]*4.5f) - op) - sp[3]);
      tmp_cells->speed4[ii] = (obstacles[ii]) ? sp[2] : sp[4] + params.omega * (w1 * local_density * (1.f + u[4] * (3.0f + u[4]*4.5f) - op) - sp[4]);
      tmp_cells->speed5[ii] = (obstacles[ii]) ? sp[7] : sp[5] + params.omega * (w2 * local_density * (1.f + u[5] * (3.0f + u[5]*4.5f) - op) - sp[5]);
      tmp_cells->speed6[ii] = (obstacles[ii]) ? sp[8] : sp[6] + params.omega * (w2 * local_density * (1.f + u[6] * (3.0f + u[6]*4.5f) - op) - sp[6]);
      tmp_cells->speed7[ii] = (obstacles[ii]) ? sp[5] : sp[7] + params.omega * (w2 * local_density * (1.f + u[7] * (3.0f + u[7]*4.5f) - op) - sp[7]);
      tmp_cells->speed8[ii] = (obstacles[ii]) ? sp[6] : sp[8] + params.omega * (w2 * local_density * (1.f + u[8] * (3.0f + u[8]*4.5f) - op) - sp[8]);
    }
  }
  return tot_u / (float)tot_cells;
}

float av_velocity(const t_param params, const t_speed* restrict cells, const int* restrict obstacles)
{

  int   tot_cells = 0;  /* no. of cells used in calculation */
  float tot_u;          /* accumulated magnitudes of velocity for each cell */

  /* initialise */
  tot_u = 0.f;

  __assume_aligned(cells        , 64);
  __assume_aligned(cells->speed0, 64);
  __assume_aligned(cells->speed1, 64);
  __assume_aligned(cells->speed2, 64);
  __assume_aligned(cells->speed3, 64);
  __assume_aligned(cells->speed4, 64);
  __assume_aligned(cells->speed5, 64);
  __assume_aligned(cells->speed6, 64);
  __assume_aligned(cells->speed7, 64);
  __assume_aligned(cells->speed8, 64);
  __assume_aligned(obstacles    , 64);
  /* loop over all non-blocked cells */
  for (int jj = 0; jj < params.ny; jj++)
  {
    #pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    {
      /* ignore occupied cells */
      if (!obstacles[ii + jj*params.nx])
      {
        /* local density total */
        const float local_density = cells->speed0[ii + jj *params.nx]
                                  + cells->speed3[ii + jj *params.nx]
                                  + cells->speed4[ii + jj *params.nx]
                                  + cells->speed1[ii + jj *params.nx]
                                  + cells->speed2[ii + jj *params.nx]
                                  + cells->speed7[ii + jj *params.nx]
                                  + cells->speed8[ii + jj *params.nx]
                                  + cells->speed5[ii + jj *params.nx]
                                  + cells->speed6[ii + jj *params.nx];

        /* x-component of velocity */
        const float u_x =(cells->speed1[ii + jj *params.nx]
                        + cells->speed5[ii + jj *params.nx]
                        + cells->speed8[ii + jj *params.nx]
                        -(cells->speed3[ii + jj *params.nx]
                        + cells->speed6[ii + jj *params.nx]
                        + cells->speed7[ii + jj *params.nx]))
                        / local_density;
        /* y-component of velocity */
        const float u_y =(cells->speed2[ii + jj *params.nx]
                        + cells->speed5[ii + jj *params.nx]
                        + cells->speed6[ii + jj *params.nx]
                        -(cells->speed4[ii + jj *params.nx]
                        + cells->speed7[ii + jj *params.nx]
                        + cells->speed8[ii + jj *params.nx]))
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

int initialise(const char* restrict paramfile, const char* restrict obstaclefile,
               t_param* restrict params,
               int** restrict obstacles_ptr, float** restrict av_vels_ptr)
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


  /* the map of obstacles */
  *obstacles_ptr = _mm_malloc(sizeof(int) * (params->ny * params->nx) ,64);

  if (*obstacles_ptr == NULL) die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  /* first set all cells in obstacle array to zero */
  #pragma omp parallel for
  for (int jj = 0; jj < params->ny; jj++)
  {
    #pragma omp simd
    for (int ii = 0; ii < params->nx; ii++)
    {
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
  *av_vels_ptr = (float*)_mm_malloc(sizeof(float) * params->maxIters,64);

  return EXIT_SUCCESS;
}

int finalise(const t_param* restrict params, t_speed* restrict cells_ptr, t_speed* restrict tmp_cells_ptr,
             int** restrict obstacles_ptr, float** restrict av_vels_ptr)
{
  /*
  ** free up allocated memory
  */
  _mm_free(cells_ptr->speed0);
  _mm_free(cells_ptr->speed1);
  _mm_free(cells_ptr->speed2);
  _mm_free(cells_ptr->speed3);
  _mm_free(cells_ptr->speed4);
  _mm_free(cells_ptr->speed5);
  _mm_free(cells_ptr->speed6);
  _mm_free(cells_ptr->speed7);
  _mm_free(cells_ptr->speed8);
  _mm_free(tmp_cells_ptr->speed0);
  _mm_free(tmp_cells_ptr->speed1);
  _mm_free(tmp_cells_ptr->speed2);
  _mm_free(tmp_cells_ptr->speed3);
  _mm_free(tmp_cells_ptr->speed4);
  _mm_free(tmp_cells_ptr->speed5);
  _mm_free(tmp_cells_ptr->speed6);
  _mm_free(tmp_cells_ptr->speed7);
  _mm_free(tmp_cells_ptr->speed8);
  _mm_free(*obstacles_ptr);
  _mm_free(*av_vels_ptr);

  *obstacles_ptr = NULL;
  *av_vels_ptr = NULL;



  return EXIT_SUCCESS;
}


float calc_reynolds(const t_param params, t_speed* restrict cells, int* restrict obstacles)
{
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

  return av_velocity(params, cells, obstacles) * params.reynolds_dim / viscosity;
}

int write_values(const t_param params, t_speed* restrict cells, int* restrict obstacles, float* restrict av_vels)
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
        local_density = cells->speed0[ii + jj * params.nx]
                      + cells->speed1[ii + jj * params.nx]
                      + cells->speed2[ii + jj * params.nx]
                      + cells->speed3[ii + jj * params.nx]
                      + cells->speed4[ii + jj * params.nx]
                      + cells->speed5[ii + jj * params.nx]
                      + cells->speed6[ii + jj * params.nx]
                      + cells->speed7[ii + jj * params.nx]
                      + cells->speed8[ii + jj * params.nx];

        /* compute x velocity component */
        u_x =(cells->speed1[ii + jj*params.nx]
            + cells->speed5[ii + jj*params.nx]
            + cells->speed8[ii + jj*params.nx]
            -(cells->speed3[ii + jj*params.nx]
            + cells->speed6[ii + jj*params.nx]
            + cells->speed7[ii + jj*params.nx]))
            / local_density;
        /* compute y velocity component */
        u_y =(cells->speed2[ii + jj*params.nx]
            + cells->speed5[ii + jj*params.nx]
            + cells->speed6[ii + jj*params.nx]
            -(cells->speed4[ii + jj*params.nx]
            + cells->speed7[ii + jj*params.nx]
            + cells->speed8[ii + jj*params.nx]))
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

  for (int ii = 0; ii < params.maxIters; ii++)
  {
    fprintf(fp, "%d:\t%.12E\n", ii, av_vels[ii]);
  }

  fclose(fp);

  return EXIT_SUCCESS;
}

void die(const char* restrict message, const int line, const char* restrict file)
{
  fprintf(stderr, "Error at line %d of file %s:\n", line, file);
  fprintf(stderr, "%s\n", message);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void usage(const char* restrict exe)
{
  fprintf(stderr, "Usage: %s <paramfile> <obstaclefile>\n", exe);
  exit(EXIT_FAILURE);
}
