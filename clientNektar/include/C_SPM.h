
#ifndef C_SPM_H
#define C_SPM_H

#include <mpi.h>
#include "C_SHAPE.h"
#include "mapping.h"

class C_SPM{

public:

   int num_partls;   // number of solid bodies
   double buffer_width; //parameter for computing the indicator function

   double ***Xp, ***Yp, ***Zp, 
          ***Up, ***Vp, ***Wp; //position and velocities of the solid bodies
   double **Xp_init, **Yp_init; //initial position of each particle

   int ***XPCOUNT;   //ask what is it

   double ***Ang, ***ALFA;  //angular velocity and phase

   double **gaussian;       //indicator function \in [0 1]
   double ***bound_box;     //box containing particle + 2*ksi

   double **tens, **stif;

   double ***FHx, ***FHy;
   double ***NH;
   double ***FHphix, ***FHphiy;

   int *particle_type;    
   /* types:
    0 - particle is static and haz zero translational and rotational velocities
    1 - particle is moving, no rotation
    2 - rotation only
    3 - translation + rotation
   */
   
   int *mobile;
   /*
    0 - particle is not moved by the fluid
    1 - particle is moved by the fliud
   */

   FILE **file_particleout; 

   //shape that can be represented analitically
   C_SHAPE  *shape;
  /*
   0 - not used
   1 - circle
  */

   Mapping  *mapx,  *mapy; /* Mapping in x and y                */

   //special:
   double **FHspan;
   double ***Fxext, ***Fyext;
   double ***Text;
  
#if 0

   // for indicator function computed externally we need 
   // a list of ranks (in MPI_COMM_WORLD) which will communicate the data to Nektar 
   int *List_of_partners_external;

   //index of indicator function computed externally
   int **offset_external;  
   int **mapping_unique2original_local; //mapping of points  i = 0 -> Nq_local-1;
                                       //   PHI_local[i] = PHI_unique_local[mapping_unique2original_local[i]];

   //next array is important on root only
   int **mapping_unique2original_global; //mapping of points  i = 0 -> Nq_local-1;
                                       //   PHI_semi_unique_global[i] = PHI_unique_global[mapping_unique2original_local[i]];


   int *Nq_local;            //number of quadrature points in this partition 
   int *Nq_unique_local;     //number of unique quadrature points in this partition
   int **Nq_unique_per_rank; //number of unique quadrature points per rank (of L4 communicator);

   int *Nq_global;            //total number of quadrature points (after local filtering)  
   int *Nq_unique_global;     //total number of unique quadrature points 
   MPI_Comm  *EXTERNAL_communicator;
   int  *EXTERNAL_communicator_size; //size of L4
   int  *EXTERNAL_communicator_rank; //rank in L4

      
   double **indicator_unique_global;
   double **indicator_global;
   double **indicator_unique_local;
   double **indicator_local

#endif


   int  SPM_init(int nparts, int ngrid_points);
   int  read_num_parts(char *Fname);
   int  read_init_conditions(char *name);
   int  set_init_conditions_external();
   int  read_in_particle_geometry(char *name);
   int  create_file_particle_out(char *name);

   void hist_data_out(double t);

   void write_restart(char *Fname);

   int  ReadMap   (char *name);
   void WriteMap  (char *name);
   void UpdateMap (char *name);
   void writedog2 ();
   void WriteInit(char *name);

   void set_extern_partners(int *List_of_partners_external);
   void set_bound_box_external(double *bound_box);
};

#endif
  
