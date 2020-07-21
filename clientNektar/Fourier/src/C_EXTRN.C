#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


#include "C_EXTRN.h"



void C_EXTRN::post_receive_indicator_function(){

  if (EXTERNAL_communicator_rank != 0) return;

  fprintf(stdout,"C_EXTRN::post_receive_indicator_function: posting irecv, source = %d, tag = 111\n",rank_L1_partner_external);
  recv_rqst = MPI_REQUEST_NULL;  
  MPI_Irecv(indicator_unique_global, Nq_unique_global*(1+3*vel_provided), MPI_DOUBLE,
    rank_L1_partner_external,111,MPI_COMM_WORLD,&recv_rqst); 

}

void C_EXTRN::process_received_data(){
 
   if (EXTERNAL_communicator == MPI_COMM_NULL) return;
 
   int i;

   if (EXTERNAL_communicator_rank == 0){
     MPI_Wait(&recv_rqst, MPI_STATUS_IGNORE);
     fprintf(stdout,"Nektar: C_EXTRN::process_received_data: data has been arrived\n");
  
     //process indicator function 
     for (i = 0; i < Nq_global; ++i)
        indicator_global[i] = indicator_unique_global[ mapping_unique2original_global[i] ];
     
     if (vel_provided == 1){
       double *uvw_unique_global = indicator_unique_global + Nq_unique_global;
       double *uvw_global =  indicator_global + Nq_global;
       int k;
       for (i = 0; i < Nq_global; ++i){
         k = mapping_unique2original_global[i];
         uvw_global[i*3] = uvw_unique_global[k*3];
         uvw_global[i*3+1] = uvw_unique_global[k*3+1];        
         uvw_global[i*3+2] = uvw_unique_global[k*3+2];
       }
     }

   }
   //scatter to other ranks
   MPI_Scatterv (indicator_global, recvcnt, displs, MPI_DOUBLE,
                 indicator_unique_local, Nq_unique_local, MPI_DOUBLE,
                 0, EXTERNAL_communicator);  

   for (i = 0; i < Nq_local; ++i)
     indicator_local[i] = indicator_unique_local[ mapping_unique2original_local[i] ];

   if (vel_provided == 0)  return;

   if (EXTERNAL_communicator_rank == 0){
 
       int recvcnt_uvw[EXTERNAL_communicator_size];
       int displs_uvw[EXTERNAL_communicator_size];
       displs_uvw[0] = 0;
       recvcnt_uvw[0] = recvcnt[0]*3;  
       for (i = 1; i < EXTERNAL_communicator_size; ++i){
         displs_uvw[i] = displs_uvw[i-1]+recvcnt_uvw[i-1];
         recvcnt_uvw[i] = recvcnt[i]*3;
       }
       MPI_Scatterv (indicator_global+Nq_global, recvcnt_uvw, displs_uvw, MPI_DOUBLE,
                     indicator_unique_local+Nq_unique_local, Nq_unique_local*3, MPI_DOUBLE,
                     0, EXTERNAL_communicator);
   }
   else
       MPI_Scatterv (NULL, NULL, NULL, MPI_DOUBLE,
                     indicator_unique_local+Nq_unique_local, Nq_unique_local*3, MPI_DOUBLE,
                     0, EXTERNAL_communicator);
   

   double *uvw_unique_local = indicator_unique_local  + Nq_unique_local;
   double *uvw_local  =  indicator_local + Nq_local;
   int k;
   for (i = 0; i < Nq_local; ++i){
     k = mapping_unique2original_local[i];
     uvw_local[i*3]   = uvw_unique_local[k*3];
     uvw_local[i*3+1] = uvw_unique_local[k*3+1];
     uvw_local[i*3+2] = uvw_unique_local[k*3+2];
   }
}

void C_EXTRN::update_indicator(double *ind_func){

   if (EXTERNAL_communicator == MPI_COMM_NULL) return;

   for (int i = 0; i < Nq_local; ++i)
     ind_func[offset_external[i]] = indicator_local[i];
}


void C_EXTRN::update_upf(double *u_p, double *v_p, double *w_p){

   if (EXTERNAL_communicator == MPI_COMM_NULL) return;

   if (vel_provided == 0)  return;

   int i, k;
   double *uvw_local = indicator_local + Nq_local;
   for (i = 0; i < Nq_local; ++i){
     k = offset_external[i];
     u_p[k] = uvw_local[i*3];
     v_p[k] = uvw_local[i*3+1];
     w_p[k] = uvw_local[i*3+2];
   }
}





