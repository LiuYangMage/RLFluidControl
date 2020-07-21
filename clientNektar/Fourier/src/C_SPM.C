
#include <hotel.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> 
#include "math.h"
#include "veclib.h"

#include "nektarF.h"
#include "map.h"

#include "C_SPM.h"

static int forcx, forcy;
void writeheader (FILE *fp, char *name, int lmax, double lz, int nz, 
			 int nel, int step, double t, char *typelist);

int C_SPM::SPM_init(int nparts, int ngridpoints){

  int  i;
  int  nztot = option("NZTOT");

  Xp = (double ***) malloc(nztot*sizeof(double**));
  Yp = (double ***) malloc(nztot*sizeof(double**));

  Xp_init = (double **) malloc(nztot*sizeof(double*));
  Yp_init = (double **) malloc(nztot*sizeof(double*));

  num_partls = nparts;
  for(i=0; i<nztot; ++i)
   {
    Xp[i] = dmatrix(0,num_partls-1,0,2);
    Yp[i] = dmatrix(0,num_partls-1,0,2);

    Xp_init[i] = dvector(0,num_partls-1);
    Yp_init[i] = dvector(0,num_partls-1);

    memset(Xp[i][0],'\0',num_partls*3*sizeof(double));
    memset(Yp[i][0],'\0',num_partls*3*sizeof(double));

    memset(Xp_init[i],'\0',num_partls*sizeof(double));
    memset(Yp_init[i],'\0',num_partls*sizeof(double));
   }
  

  Up = (double ***) malloc(nztot*sizeof(double**));
  Vp = (double ***) malloc(nztot*sizeof(double**));

  for(i=0; i<nztot; ++i)
   {
    Up[i] = dmatrix(0,num_partls-1,0,2);
    Vp[i] = dmatrix(0,num_partls-1,0,2);

    memset(Up[i][0],'\0',num_partls*3*sizeof(double));
    memset(Vp[i][0],'\0',num_partls*3*sizeof(double));
   }

//  XPCOUNT = imatrix(0,num_partls-1,0,2);
  XPCOUNT = (int ***) malloc(nztot*sizeof(int**));
  for(i=0; i<nztot; ++i)
   {
    XPCOUNT[i] = imatrix(0,num_partls-1,0,2);
    memset(XPCOUNT[i][0],'\0',num_partls*3*sizeof(int));
   }
//angular velocity on z direction only !
  Ang = (double ***) malloc(nztot*sizeof(double**));
  ALFA = (double ***) malloc(nztot*sizeof(double**));
 
  for(i=0; i<nztot; i++){
     Ang[i] = dmatrix(0,num_partls,0,2);
     memset(Ang[i][0],'\0',num_partls*3*sizeof(double));
  }

  for(i=0; i<nztot; i++){
     ALFA[i] = dmatrix(0,num_partls,0,2);
     memset(ALFA[i][0],'\0',num_partls*3*sizeof(double));
  } 
   
  FHx = (double ***) malloc(nztot*sizeof(double**));
  FHy = (double ***) malloc(nztot*sizeof(double**));

  FHphix = (double ***) malloc(nztot*sizeof(double**));
  FHphiy = (double ***) malloc(nztot*sizeof(double**));

  NH = (double ***) malloc(nztot*sizeof(double**));

  Fxext = (double ***) malloc(nztot*sizeof(double**));
  Fyext = (double ***) malloc(nztot*sizeof(double**));
  Text = (double ***) malloc(nztot*sizeof(double**));

  for(i=0; i<nztot; i++){

   FHx[i] = dmatrix(0,num_partls-1,0,2);
   FHy[i] = dmatrix(0,num_partls-1,0,2);
   FHphix[i] = dmatrix(0,num_partls-1,0,2);
   FHphiy[i] = dmatrix(0,num_partls-1,0,2);

   NH[i] = dmatrix(0,num_partls-1,0,2);

   Fxext[i] = dmatrix(0,num_partls-1,0,2);
   Fyext[i] = dmatrix(0,num_partls-1,0,2);
   Text[i] = dmatrix(0,num_partls-1,0,2);

 
   memset(FHx[i][0],'\0',num_partls*3*sizeof(double));
   memset(FHy[i][0],'\0',num_partls*3*sizeof(double));
   memset(FHphix[i][0],'\0',num_partls*3*sizeof(double));
   memset(FHphiy[i][0],'\0',num_partls*3*sizeof(double));

   memset(NH[i][0],'\0',num_partls*3*sizeof(double));
   memset(Fxext[i][0],'\0',num_partls*2*sizeof(double)); 
   memset(Fyext[i][0],'\0',num_partls*2*sizeof(double)); 
   memset(Text[i][0],'\0',num_partls*2*sizeof(double));

  }


  stif = (double **) malloc(num_partls*sizeof(double*));
  tens = (double **) malloc(num_partls*sizeof(double*));
  for(i=0; i<num_partls; i++){
   stif[i] = dvector(0,nztot-1);
   tens[i] = dvector(0,nztot-1);
   memset(stif[i],'\0',nztot*sizeof(double));
   memset(tens[i],'\0',nztot*sizeof(double));
  }

  gaussian = dmatrix(0,num_partls-1,0,ngridpoints-1);
  memset(gaussian[0],'\0',num_partls*ngridpoints*sizeof(double));

  /* types:
  0 - particle is static and haz zero translational and rotational velocities
  1 - particle is moving, no rotation
  2 - rotation only
  3 - translation + rotation 
  */

  /*  
   0 - PLOT3D 
   1 - BACKGROUNDGRID
   2 - SHAPE
   3 - External Field (LAMMPS)
   4 - FWMAV (Flipping Wing Micro Air Vehicle) 
   5 - FARING
  */

  shape = new C_SHAPE[num_partls];

  mapx  = new Mapping[num_partls];
  mapy  = new Mapping[num_partls];
//  mapx    = (Mapping *) calloc (num_partls, sizeof(Mapping));
//  mapy    = (Mapping *) calloc (num_partls, sizeof(Mapping));
  

  for(i=0; i<num_partls; ++i)
   {
     mapx[i].NZ = option("NZTOT");
     mapy[i].NZ = option("NZTOT");

     mapx[i].d    = (double *) calloc (nztot, sizeof(double)); /* displacement     */
     mapx[i].t    = (double *) calloc (nztot, sizeof(double)); /* t   - derivative */
     mapx[i].tt   = (double *) calloc (nztot, sizeof(double)); /* tt  - derivative */
     mapx[i].f    = (double *) calloc (nztot, sizeof(double)); /* forcing          */

     mapy[i].d    = (double *) calloc (nztot, sizeof(double)); /* displacement     */
     mapy[i].t    = (double *) calloc (nztot, sizeof(double)); /* t   - derivative */
     mapy[i].tt   = (double *) calloc (nztot, sizeof(double)); /* tt  - derivative */
     mapy[i].f    = (double *) calloc (nztot, sizeof(double)); /* forcing          */
   }

  mobile = ivector(0,num_partls-1);
  /*
    0 - particle is not moved by the fluid 
    1 - particle is moved by the fliud
  */
  memset(mobile,'\0',num_partls*sizeof(int));

  file_particleout = new FILE*[num_partls];  


  bound_box = new double**[num_partls];

  forcx = (int) dparam("FORCX");
  forcy = (int) dparam("FORCY");
  
  return 0;
}

int C_SPM::create_file_particle_out(char *name){


   int procid = option("PROCID");
   if (procid != 0) 
     return 0;

   //check if file already exists
   struct stat buffer ;
   char f_name[BUFSIZ];
   int part;

   for (part = 0; part < num_partls; ++part){

     sprintf(f_name, "%s_particle_%d.out",name,part);

     if ( stat( f_name, &buffer ) == 0){
      
       file_particleout[part] = fopen(f_name,"a");
       if (file_particleout[part] == NULL){
         fprintf(stderr,"failed to open file %s\n",f_name);  
         return -1;
       }
     }
     else{ //file does not exist

       file_particleout[part] = fopen(f_name,"w");
       if (file_particleout[part] == NULL){
         fprintf(stderr,"failed to open file %s\n",f_name);
         return -1;
       }
      fprintf(file_particleout[part],"VARIABLES=\"Time\",  \"Xp\",  \"Yp\", \"ALFAZ\",  \"Up\", \"Vp\", \"Angz\", \"Fhx\", \"Fhy\", \"Nhz\", \"Fhnewx\",  \"Fhnewy\"\n");

     } 
  }
  return 0;
}


int C_SPM::read_init_conditions(char *name){

  int i,j;
  int  nztot = option("NZTOT");
  char Fname[FILENAME_MAX];
  sprintf(Fname, "%s_init.dat", name); // PARTICLE TRAJECTORY FILE
  FILE * file_particlein = fopen(Fname, "r");
  if ( file_particlein == NULL){
    fprintf(stderr,"C_SPM::read_init_conditions: file %s is not found\n",Fname);
    return -1;
  }

 for(j=0; j<nztot; j++) 
  for(i=0; i<num_partls; i++) 
     fscanf(file_particlein,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",Xp[j][i],Yp[j][i],ALFA[j][i],Up[j][i],Vp[j][i],Ang[j][i],Fxext[j][i],Fyext[j][i], Text[j][i]);

 fclose(file_particlein);

 for(j=0; j<nztot; j++){ 
  for(i=0; i<num_partls; i++)  {
    Up[j][i][2] = Up[j][i][1] = Up[j][i][0]; 
    Vp[j][i][2] = Vp[j][i][1] = Vp[j][i][0]; 
    Ang[j][i][2] = Ang[j][i][1] = Ang[j][i][0];
//it could be replaceb by SPM_readMap
    mapx[i].d[j] = Xp[j][i][0];
    mapx[i].t[j] = Up[j][i][0];
    mapy[i].d[j] = Yp[j][i][0];
    mapy[i].t[j] = Vp[j][i][0];

    Xp_init[j][i] = Xp[j][i][0];
    Yp_init[j][i] = Yp[j][i][0];
  }
 }
  return 0;
}

#if 0
int C_SPM::read_num_parts(char *Fname){
//first line in the input file contains number of particles 

  FILE *DB_File;
  DB_File = fopen(Fname,"r");
  if (DB_File == NULL){
    fprintf(stderr,"C_SPM::read_num_parts: file %s is not found\n",Fname);
    return -1;
  }
  int Nshapes;
  fscanf(DB_File,"%d",&Nshapes);
  fclose(DB_File);
  return Nshapes;
}


void C_SPM::write_restart(char *Fname){

  FILE *pFILE;
  pFILE = fopen(Fname,"w");
  int shp_ID;
   
  fprintf(pFILE,"%d\n",num_partls); 
  for (shp_ID = 0; shp_ID < num_partls; ++shp_ID){
    if (SURF_MODEL[shp_ID] == 2){
      double *data = shape[shp_ID].data;
      if (strcmp(shape[shp_ID].type,"cyl") == 0){
        fprintf(pFILE,"CYLINDER\n");
        fprintf(pFILE,"%d\n",8);
        fprintf(pFILE,"NORMAL\n");
        fprintf(pFILE,"%2.14f %2.14f %2.14f\n",data[0],data[1],data[2]);
        fprintf(pFILE,"FACECENTER\n");
        fprintf(pFILE,"%2.14f %2.14f %2.14f\n",data[3],data[4],data[5]);
        fprintf(pFILE,"LENGTH\n");
        fprintf(pFILE,"%2.14f\n",data[6]);
        fprintf(pFILE,"RADIUS\n");
        fprintf(pFILE,"%2.14f\n",data[7]);
        fprintf(pFILE,"THICKNESS\n");
        fprintf(pFILE,"%2.14f\n",data[8]); 
         fprintf(pFILE,"VELOCITY\n");
        fprintf(pFILE,"%2.14f %2.14f %2.14f\n",data[9],data[10],data[11]);
        fprintf(pFILE,"MASSRATIO\n");
        fprintf(pFILE,"%2.14f\n",data[12]);
        fprintf(pFILE,"MOBILE\n");
        fprintf(pFILE,"%d\n",mobile[shp_ID]);
      }
    }
  } 
  fclose(pFILE);   
}


#endif


int C_SPM::read_in_particle_geometry(char *name){


  char buf[BUFSIZ];
  char tag;
  char F_DB_name[BUFSIZ];
  int  Nshapes,Nplanes;
  int  nztot = option("NZTOT");
  double LZ = dparam("LZ");

  char Fname[FILENAME_MAX];
  sprintf(Fname, "%s_particle.geo", name); // PARTICLE TRAJECTORY FILE
  FILE *IN_File;
  IN_File = fopen(Fname,"r");
  if (IN_File == NULL){
    fprintf(stderr,"C_SPM::read_in_particle_geometry: file %s is not found\n",Fname);
    return -1;
  }

  fscanf(IN_File,"%d",&Nplanes);
//  fscanf(IN_File,"%d",&Nshapes);

  if(Nplanes != nztot)
    fprintf(stderr,"Warning nz != Nplanes \n");

  int shp_ID;
  int nparams;
  char param[512];

//  for (shp_ID = 0; shp_ID < Nshapes; ++shp_ID){
  for (shp_ID = 0; shp_ID < num_partls; ++shp_ID){

    fscanf(IN_File,"%s",buf);
  
    if (strcmp(buf,"CIRCLE")== 0){

       sprintf(shape[shp_ID].type,"cir"); 
       shape[shp_ID].data = new double[9*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',9*nztot*sizeof(double));
 //read initial particle geometry for each plane     
      for(int k = 0; k<Nplanes; ++k)
       {
        int offset = k*9;
        fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
        
//         if (strcmp(param,"FACECENTER")==0) 
//            fscanf(IN_File,"%lf %lf",&data[offset+0],&data[offset+1]); /*center of the bottom face */
         if (strcmp(param,"RADIUS")==0)
            fscanf(IN_File,"%lf",&data[offset+2]);  /* radius  */
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[offset+3]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf",&data[offset+4],&data[offset+5]); /* velocity of the center in X-,Y- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[offset+6]); /* mass ratio */
         else if (strcmp(param,"TENSION")==0)
             fscanf(IN_File,"%lf",&data[offset+7]); /* tension */
         else if (strcmp(param,"STIFFNESS")==0)
             fscanf(IN_File,"%lf",&data[offset+8]); /* stiffness */
        }
         data[offset+0] = Xp[k][shp_ID][0];
         data[offset+1] = Yp[k][shp_ID][0];
         
//         tens[k][shp_ID] = dparam("WNC")*dparam("WNC"); 
//         stif[k][shp_ID] = dparam("WNB")*dparam("WNB");
         tens[shp_ID][k] = data[offset+7]; 
         stif[shp_ID][k] = data[offset+8];
       }
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */
 
//         if(mobile[shp_ID])
//           iparam_set("MOVING_PARTLS", 1);
//         data[0] = Xp[0][shp_ID][0];
//         data[1] = Yp[0][shp_ID][0];
/*
       for(int k = 1; k<nztot; ++k)
        {
          offset = k*7;
          data[offset+0] = data[0];
          data[offset+1] = data[1];
          data[offset+2] = data[2];
          data[offset+3] = data[3];
          data[offset+4] = data[4];
          data[offset+5] = data[5];
          data[offset+6] = data[6];
        }
*/
//         if (mobile[shp_ID] == 1)
//            particle_type[shp_ID] = 1;

//      fprintf(stdout,"part %d  center = %g, %g  radius = %g  thickness %g  velocity = %g,%g massratio = %g \n",
//                         shp_ID,data[0],data[1],data[2],data[3],data[4],data[5],data[6]); 

      shape[shp_ID].bound_box = new double[4*nztot];
      
      for(int k = 0; k<nztot; ++k)
       {
          int offset1 = k*4;
          int offset2 = k*7;
          shape[shp_ID].bound_box[offset1+0] = data[offset2+0] - 5.0*data[offset2+2]; //x_min
          shape[shp_ID].bound_box[offset1+1] = data[offset2+0] + 5.0*data[offset2+2]; //x_max
          shape[shp_ID].bound_box[offset1+2] = data[offset2+1] - 5.0*data[offset2+2]; //y_min
          shape[shp_ID].bound_box[offset1+3] = data[offset2+1] + 5.0*data[offset2+2]; //y_max
       }

    }
    else if (strcmp(buf,"RECTANGLE")== 0){
       sprintf(shape[shp_ID].type,"rect"); 
       shape[shp_ID].data = new double[12*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',12*nztot*sizeof(double));
       
       fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"VERTICES")==0) 
            fscanf(IN_File,"%lf %lf %lf %lf %lf %lf %lf %lf",
                &data[0],&data[1],&data[2],&data[3],&data[4],&data[5],&data[6],&data[7]); /*vertices of the rectange */
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[8]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf",&data[9],&data[10]); /* velocity of the center in X-,Y- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[11]); /* mass ratio */
        }
        
        fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */


       for(int k = 1; k<nztot; ++k)
        {
          int offset = k*12;
          data[offset+0] = data[0];
          data[offset+1] = data[1];
          data[offset+2] = data[2];
          data[offset+3] = data[3];
          data[offset+4] = data[4];
          data[offset+5] = data[5];
          data[offset+6] = data[6];
          data[offset+7] = data[7];
          data[offset+8] = data[8];
          data[offset+9] = data[9];
          data[offset+10] = data[10];
          data[offset+11] = data[11];
        }

    }
    else if (strcmp(buf,"BUOYANCY")== 0){

       sprintf(shape[shp_ID].type,"buoy"); 
       shape[shp_ID].data = new double[9*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',9*nztot*sizeof(double));
 //read initial particle geometry for each plane     
      for(int k = 0; k<Nplanes; ++k)
       {
        int offset = k*9;
        fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
        
         if (strcmp(param,"RADIUS")==0)
            fscanf(IN_File,"%lf",&data[offset+2]);  /* radius  */
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[offset+3]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf",&data[offset+4],&data[offset+5]); /* velocity of the center in X-,Y- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[offset+6]); /* mass */
         else if (strcmp(param,"TENSION")==0)
             fscanf(IN_File,"%lf",&data[offset+7]); /* tension */
         else if (strcmp(param,"STIFFNESS")==0)
             fscanf(IN_File,"%lf",&data[offset+8]); /* stiffness */
        }

         data[offset+0] = Xp[k][shp_ID][0];
         data[offset+1] = Yp[k][shp_ID][0];
       
//        ROOT fprintf(stderr, " radius = %g thikness = %g velocity =[ %g  %g] mass = %g  tension = %g  stifness = %g \n",
//                     data[offset+2],data[offset+3],data[offset+4],data[offset+5],data[offset+6],data[offset+7],data[offset+8]);

 // insert tension and stiffness
 // zwang 20170523
 // not sure if it is right, check it later. 
        tens[shp_ID][k] = data[offset+7]; 
        stif[shp_ID][k] = data[offset+8];
       }
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */
        
 //      if(mobile[shp_ID])
//         iparam_set("MOVING_PARTLS", 1);

//      fprintf(stdout,"part %d  center = %g, %g  radius = %g  thickness %g  velocity = %g,%g massratio = %g \n",
//                         shp_ID,data[0],data[1],data[2],data[3],data[4],data[5],data[6]); 

      shape[shp_ID].bound_box = new double[4*nztot];
      
      for(int k = 0; k<nztot; ++k)
       {
          int offset1 = k*4;
          int offset2 = k*5;
          shape[shp_ID].bound_box[offset1+0] = data[offset2+0] - 5.0*data[offset2+2]; //x_min
          shape[shp_ID].bound_box[offset1+1] = data[offset2+0] + 5.0*data[offset2+2]; //x_max
          shape[shp_ID].bound_box[offset1+2] = data[offset2+1] - 5.0*data[offset2+2]; //y_min
          shape[shp_ID].bound_box[offset1+3] = data[offset2+1] + 5.0*data[offset2+2]; //y_max
       }

    }
    else if (strcmp(buf,"STRAKE")== 0){

       sprintf(shape[shp_ID].type,"strake"); 
       shape[shp_ID].data = new double[11*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',11*nztot*sizeof(double));
 //read initial particle geometry for each plane     
      for(int k = 0; k<Nplanes; ++k)
       {
        int offset = k*11;
        fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
        
         if (strcmp(param,"RADIUS")==0)
            fscanf(IN_File,"%lf",&data[offset+2]);  /* radius  */
         else if (strcmp(param,"HEIGHT")==0)
            fscanf(IN_File,"%lf",&data[offset+3]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"ANGLE")==0)
            fscanf(IN_File,"%lf",&data[offset+4]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[offset+5]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf",&data[offset+6],&data[offset+7]); /* velocity of the center in X-,Y- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[offset+8]); /* mass */
         else if (strcmp(param,"TENSION")==0)
             fscanf(IN_File,"%lf",&data[offset+9]); /* tension */
         else if (strcmp(param,"STIFFNESS")==0)
             fscanf(IN_File,"%lf",&data[offset+10]); /* stiffness */
        }

         data[offset+0] = Xp[k][shp_ID][0];
         data[offset+1] = Yp[k][shp_ID][0];
       
//        ROOT fprintf(stderr, " radius = %g thikness = %g velocity =[ %g  %g] mass = %g  tension = %g  stifness = %g \n",
//                     data[offset+2],data[offset+3],data[offset+4],data[offset+5],data[offset+6],data[offset+7],data[offset+8]);

 // insert tension and stiffness
 // zwang 20170523
 // not sure if it is right, check it later. 
        tens[shp_ID][k] = data[offset+9]; 
        stif[shp_ID][k] = data[offset+10];
       }
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */
        
 //      if(mobile[shp_ID])
//         iparam_set("MOVING_PARTLS", 1);

//      fprintf(stdout,"part %d  center = %g, %g  radius = %g  thickness %g  velocity = %g,%g massratio = %g \n",
//                         shp_ID,data[0],data[1],data[2],data[3],data[4],data[5],data[6]); 

      shape[shp_ID].bound_box = new double[4*nztot];
      
      for(int k = 0; k<nztot; ++k)
       {
          int offset1 = k*4;
          int offset2 = k*5;
          shape[shp_ID].bound_box[offset1+0] = data[offset2+0] - 5.0*data[offset2+2]; //x_min
          shape[shp_ID].bound_box[offset1+1] = data[offset2+0] + 5.0*data[offset2+2]; //x_max
          shape[shp_ID].bound_box[offset1+2] = data[offset2+1] - 5.0*data[offset2+2]; //y_min
          shape[shp_ID].bound_box[offset1+3] = data[offset2+1] + 5.0*data[offset2+2]; //y_max
       }

    }
    else if (strcmp(buf,"SPHERE")== 0){

       sprintf(shape[shp_ID].type,"sph"); 
       shape[shp_ID].data = new double[9*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',9*nztot*sizeof(double));
 //read initial particle geometry for each plane     
      for(int k = 0; k<Nplanes; ++k)
       {
        int offset = k*9;
        fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
        
         if (strcmp(param,"RADIUS")==0)
            fscanf(IN_File,"%lf",&data[offset+3]);  /* radius  */
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[offset+4]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf %lf",&data[offset+5],&data[offset+6], &data[offset+7]); /* velocity of the center in X-,Y-,Z- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[offset+8]); /* mass ratio */
        }

         data[offset+0] = Xp[k][shp_ID][0];
         data[offset+1] = Yp[k][shp_ID][0];
         data[offset+2] = LZ/2.; // z is allways LZ/2 !

 // insert tension and stiffness
 // zwang 20170516
 // not sure if it is right, check it later. 
        tens[shp_ID][k] = dparam("WNC")*dparam("WNC"); 
        stif[shp_ID][k] = dparam("WNB")*dparam("WNB");
       }
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */
        
 //      if(mobile[shp_ID])
//         iparam_set("MOVING_PARTLS", 1);

//      fprintf(stdout,"part %d  center = %g, %g  radius = %g  thickness %g  velocity = %g,%g massratio = %g \n",
//                         shp_ID,data[0],data[1],data[2],data[3],data[4],data[5],data[6]); 

      shape[shp_ID].bound_box = new double[4*nztot];
/*      
      for(int k = 0; k<nztot; ++k)
       {
          int offset1 = k*4;
          int offset2 = k*7;
          shape[shp_ID].bound_box[offset1+0] = data[offset2+0] - 5.0*data[offset2+2]; //x_min
          shape[shp_ID].bound_box[offset1+1] = data[offset2+0] + 5.0*data[offset2+2]; //x_max
          shape[shp_ID].bound_box[offset1+2] = data[offset2+1] - 5.0*data[offset2+2]; //y_min
          shape[shp_ID].bound_box[offset1+3] = data[offset2+1] + 5.0*data[offset2+2]; //y_max
       }
*/

    }
    else if (strcmp(buf,"CABLE2D")== 0){
       sprintf(shape[shp_ID].type,"cable2d"); 
       shape[shp_ID].data = new double[6*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',6*nztot*sizeof(double));
       
      for(int k = 0; k<Nplanes; ++k)
      {
        int offset = k*6;
       fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"POSITION")==0) 
            fscanf(IN_File,"%lf %lf",
                &data[offset+0],&data[offset+1]); /*vertices of the rectange */
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[offset+2]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf",&data[offset+3],&data[offset+4]); /* velocity of the center in X-,Y- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[offset+5]); /* mass ratio */
        }
      }
        
        fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */


       for(int k = 1; k<nztot; ++k)
        {
          int offset = k*6;
          data[offset+0] = data[0];
          data[offset+1] = data[1];
          data[offset+2] = data[2];
          data[offset+3] = data[3];
          data[offset+4] = data[4];
          data[offset+5] = data[5];
        }

    }
    else if (strcmp(buf,"WALL")== 0){

       sprintf(shape[shp_ID].type,"wall"); 
       shape[shp_ID].data = new double[10*nztot];
       double *data = shape[shp_ID].data;
       memset(data,'\0',10*nztot*sizeof(double));
 //read initial particle geometry for each plane     
      for(int k = 0; k<Nplanes; ++k)
       {
        int offset = k*10;
        fscanf(IN_File,"%d",&nparams);
        for (int i = 0; i < nparams; ++i){
         fscanf(IN_File,"%s",&param[0]);
        
         if (strcmp(param,"ZPOSITION")==0)
            fscanf(IN_File,"%lf",&data[offset+2]);  /* radius  */
         else if (strcmp(param,"WALL_THICK")==0)
            fscanf(IN_File,"%lf",&data[offset+3]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"THICKNESS")==0)
            fscanf(IN_File,"%lf",&data[offset+4]);  /* thickness of indicator function (ksi)*/
         else if (strcmp(param,"VELOCITY")==0)
            fscanf(IN_File,"%lf %lf",&data[offset+5],&data[offset+6]); /* velocity of the center in X-,Y- directions */
         else if (strcmp(param,"MASS")==0)
             fscanf(IN_File,"%lf",&data[offset+7]); /* mass */
         else if (strcmp(param,"TENSION")==0)
             fscanf(IN_File,"%lf",&data[offset+8]); /* tension */
         else if (strcmp(param,"STIFFNESS")==0)
             fscanf(IN_File,"%lf",&data[offset+9]); /* stiffness */
        }

         data[offset+0] = Xp[k][shp_ID][0];
         data[offset+1] = Yp[k][shp_ID][0];

//        ROOT fprintf(stderr, " z_position = %g wall_thick = %g  tension = %g  stifness = %g \n",
//                     data[offset+2],data[offset+3],data[offset+7],data[offset+8]);

        tens[shp_ID][k] = data[offset+8]; 
        stif[shp_ID][k] = data[offset+9];
       }
         fscanf(IN_File,"%s",&param[0]);
         if (strcmp(param,"MOBILE")==0)
             fscanf(IN_File,"%d",&mobile[shp_ID]);  /* 1 for flow induced motion, zero for the rest */


      shape[shp_ID].bound_box = new double[4*nztot];
      
      for(int k = 0; k<nztot; ++k)
       {
          int offset1 = k*4;
          int offset2 = k*5;
          shape[shp_ID].bound_box[offset1+0] = data[offset2+0] - 5.0*data[offset2+2]; //x_min
          shape[shp_ID].bound_box[offset1+1] = data[offset2+0] + 5.0*data[offset2+2]; //x_max
          shape[shp_ID].bound_box[offset1+2] = data[offset2+1] - 5.0*data[offset2+2]; //y_min
          shape[shp_ID].bound_box[offset1+3] = data[offset2+1] + 5.0*data[offset2+2]; //y_max
       }

    }
    else{
      fprintf(stderr,"C_SPM::read_in_particle_geometry: data file format is unrecognized [ %s ]\n",buf);
      return -1;
    }
 }
 fclose(IN_File); 

 return 0;
}

void  C_SPM::hist_data_out( double time){

  int i,k; 
  int  nztot = option("NZTOT");

  for(k=0; k<nztot; ++k)
   for(i=0; i < num_partls; i++)  {
    fprintf(file_particleout[i],"%lf   %lf %lf   %lf   %lf %lf   %lf   %lf %lf   %lf   %lf %lf  %d\n",
          time, Xp[k][i][0], Yp[k][i][0], 
          ALFA[k][i][0],
          Up[k][i][0],       Vp[k][i][0],
          Ang[k][i][0],
          FHx[k][i][0],      FHy[k][i][0],
          NH[k][i][0],      
          FHphix[k][i][0],   FHphiy[k][i][0],k);
    fflush(file_particleout[i]);
   }
 }

void SPM_readMap (FILE *fp, Mapping *map);

int C_SPM::ReadMap (char *name)
{
  char  fname[FILENAME_MAX];
  FILE *fp;
  int i,j,k;
  
  for(i=0; i<num_partls;++i)
  {
   sprintf (fname, "%s.map.%d.rst", name,i);

   if ((fp = fopen (fname,"r")) == (FILE *) NULL) {
    ErrorHandler ("ReadMap", "failed to open the map file", WARNING);
    return 0;
   } else {
    SPM_readMap (fp, &mapx[i]);
    SPM_readMap (fp, &mapy[i]);
    fclose  (fp);
  }
   // check for consistent times 

  if ((mapx[i].time != mapy[i].time)) {
      ROOT fprintf(stderr, "Different start times for mapx and mapy!\n");
    // make them agree
    mapy[i].time = mapx[i].time;
 	
   }
 }

  return 1;
}


void SPM_readMap (FILE *fp, Mapping *map)
{
  char   buf[BUFSIZ];
  int    k, NZ, nZ;
  double time, dummy;

  /* read up to -- and then past -- next set of comment lines */
  
  while ((*buf = getc(fp)) != '#') fgets (buf, BUFSIZ, fp); ungetc (*buf, fp); 
  while ((*buf = getc(fp)) == '#') fgets (buf, BUFSIZ, fp); ungetc (*buf, fp); 
  
  fscanf (fp, "%lf", &time);
  map->time = time;

  fscanf (fp, "%d",  &NZ);
  nZ = min(NZ,map->NZ); // To avoid overwrites
  
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->d   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy); // skip over rest of line
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->t   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->tt  + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->f   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);

  if (option("nomean")) {
    /* We explicitly force all of the variables below to be zero */
    map->d[0] = map->t[0] = map->tt[0] = 0.;
  }

  return;
}

static FILE  **Fp, **Fp2;
static int step = 0;

void C_SPM::UpdateMap (char *name)
{
  static double ampx, freqx, phizx, phitx, ampy, freqy, phizy, phity,
                betax, betay, wn, wnc, wnb, mass_, zeta_, dt, ramp, zramp;
  double        mass, zeta;
  
  int   nel = iparam("NUM_ELEMENTS");
  int   i,j,k;
  int   nz = option("NZ"); 
  int   nztot = option("NZTOT"); 
  int   mstep = max(iparam("MSTEP"),5);

//  if (step == 0) {                               /* first time through */
//   Fp   = (FILE **) malloc(sizeof(FILE*)*num_partls);
//   Fp2  = (FILE **) malloc(sizeof(FILE*)*num_partls);
//  }

  for(i = 0; i<num_partls; ++i)
  {
   if (step == 0) {                               /* first time through */
    char fname[FILENAME_MAX];
    ampx = dparam("AMPX"); freqx = dparam("FREQX"); phizx = dparam("PHIZX"); phitx = dparam("PHITX"); betax = dparam("BETAX");
    ampy = dparam("AMPY"); freqy = dparam("FREQY"); phizy = dparam("PHIZY"); phity = dparam("PHITY"); betay = dparam("BETAY");

    wn   = dparam("WN");   wnc   = dparam("WNC");   wnb   = dparam("WNB");   mass_ = dparam("ZMASS");
    zeta_= dparam("ZETA"); dt    = dparam("DELT");  ramp  = dparam("RAMP");  zramp = dparam("ZRAMP");
   

/*
    ROOT {
      sprintf (fname, "%s.%d.dog", name,i);
      if ((Fp[i] = fopen (fname,"w")) == (FILE *) NULL)
	   ErrorHandler ("UpdateMap", "failed to open the force file", ERROR);
      writeheader (Fp[i], name, LGmax, dparam("LZ"), nztot,
		   nel, step, mapx[i].time, "txyDLpk");
      sprintf (fname, "%s.cab", name);
      if ((Fp2[i] = fopen (fname,"w")) == (FILE *) NULL)
	   ErrorHandler ("UpdateMap", "failed to open the cable energy file", ERROR);
      writeheader (Fp2[i], name, LGmax, dparam("LZ"), nztot,
		  nel, step, mapx[i].time, "t Ex Ey");
    }
*/
  }
 
  if (i == 0)
  step++;

  mapx[i].time += dt;
  mapy[i].time += dt;


  mass = mass_;                                   /* variable mass/zeta at startup */
  zeta = zeta_;
  if (step*dt < ramp) {
    mass = mass_ / (0.5 * (1.0 - cos(M_PI*step*dt/ramp)) + 1e-6);
    ROOT printf ("ramp: mass = %g\n", mass);
  }
  if (step*dt < zramp) {
    double expon = exp(-7.0*step*dt/zramp);
    zeta = 1.0 * expon + zeta_ * (1.0 - expon);
    ROOT printf ("ramp: zeta = %g\n", zeta);
  }
  /* update motion */
  if (forcx > 0) 
    update_forc (&mapx[i], ampx, freqx, phizx, phitx, betax, forcx); 
  else {
//    filter      (mapx, 'x');
  if(option("oldupdate"))
   {
    fprintf(stderr,"update_free() not works, please don't use option '-ou' \n");
    exit(-1);
    update_free (&mapx[i], dt, mass, wn, wnc, wnb, zeta, forcx);
   }
  else
    update_free (&mapx[i], dt, mass, wn, wnc, wnb, zeta, forcx, tens[i], stif[i]);	
  }
  if (forcy > 0)
    update_forc (&mapy[i], ampy, freqy, phizy, phity, betay, forcy);
  else {
//    filter      (mapy, 'y');
  if(option("oldupdate"))
   {
    fprintf(stderr,"update_free() not works, please don't use option '-ou' \n");
    exit(-1);
    update_free (&mapy[i], dt, mass, wn, wnc, wnb, zeta, forcy);
  }
  else
   update_free (&mapy[i], dt, mass, wn, wnc, wnb, zeta, forcy, tens[i], stif[i]);
  }

 // copy displacement, velocity to shape.data
   for(k=0; k<nztot; ++k)
   {
       int offset = k*9;
    if ( (strcmp(shape[i].type,"cir") == 0) || (strcmp(shape[i].type,"buoy") == 0) )
       {
         shape[i].data[offset+0] = mapx[i].d[k];
         shape[i].data[offset+1] = mapy[i].d[k];
         
//         ROOT fprintf(stderr,"dx = %lf dy = %lf  %d \n",shape[i].data[offset+0],shape[i].data[offset+1],k);
       }
    else
       {
         fprintf(stderr,"Not Implemented Yet ! Updatemap ... \n");
         exit(-1);
       }


      Xp[k][i][0] = mapx[i].d[k];
      Yp[k][i][0] = mapy[i].d[k];

      Up[k][i][0] = mapx[i].t[k];
      Vp[k][i][0] = mapy[i].t[k];

      if(Xp_init[k][i] != 0.)
        {
          Xp[k][i][0]  += Xp_init[k][i];
          mapx[i].d[k] += Xp_init[k][i];
        }

      if(Yp_init[k][i] != 0.)
        {
          Yp[k][i][0]  += Yp_init[k][i];
          mapy[i].d[k] += Yp_init[k][i];
        }

   }

  } // end of all particls

  return;
}


void SPM_writeMap (FILE *fp, Mapping *map);

/* ------------------------------------------------------------------------ *
 * WriteMap() - write the map file (for restarts etc.)
 * ------------------------------------------------------------------------ */

void C_SPM::WriteMap (char *name)
{
  char  fname[FILENAME_MAX];
  FILE *fp;
  int   i,j,k;
  int   nel = iparam("NUM_ELEMENTS");
  static int mapstep = 0, map_number = 0;

  for(i=0; i<num_partls; ++i)
  {

   if (option("SLICES")) {
    sprintf (fname, "%s_%d.%d.map",  name,map_number,i);

    if(i == 0)
    ++map_number;

    if ((fp = fopen (fname,"w")) == (FILE *) NULL)
      ErrorHandler ("WriteMap", "failed to open the map file", ERROR);
   } else {
    sprintf (fname, "%s.%d.map", name,i);
    int err = backup (fname);
    if (err != 0 && mapstep != 0) 
      ErrorHandler ("WriteMap", "failed to backup the map file", WARNING);
    if ((fp = fopen (fname,"w")) == (FILE *) NULL)
      ErrorHandler ("WriteMap", "failed to open the map file", ERROR);
   }
  /* Write the header just to be safe */

  if (i==0) // first particle oonly
  mapstep = min(iparam("NSTEPS"),mapstep+iparam("IOSTEP"));

  writeheader (fp,
	       name,
	       LGmax, 
	       dparam("LZ"),
	       option("NZTOT"),
	       nel,
	       mapstep,                   /* ignore number of steps */
	       mapx[i].time,
	       "maps");
  
  /* Write additional information */
  
  fprintf (fp, "%-25.6g Forcx\n", dparam("FORCX"));
  fprintf (fp, "%-25.6g Ampx \n", dparam("AMPX"));
  fprintf (fp, "%-25.6g Freqx\n", dparam("FREQX"));
  fprintf (fp, "%-25.6g Phizx\n", dparam("PHIZX"));
  fprintf (fp, "%-25.6g Phitx\n", dparam("PHITX"));
  fprintf (fp, "%-25.6g Forcy\n", dparam("FORCY"));
  fprintf (fp, "%-25.6g Ampy \n", dparam("AMPY"));
  fprintf (fp, "%-25.6g Freqy\n", dparam("FREQY"));
  fprintf (fp, "%-25.6g Phizy\n", dparam("PHIZY"));
  fprintf (fp, "%-25.6g Phity\n", dparam("PHITY"));
 
  fprintf (fp, "%-25.6g Natural Frequency  \n", dparam("WN"));
  fprintf (fp, "%-25.6g Mass               \n", dparam("ZMASS"));
  fprintf (fp, "%-25.6g Damping Coefficient\n", dparam("ZETA"));
  fprintf (fp, "%-25.6g Phase Speed cable  \n", dparam("WNC"));
  fprintf (fp, "%-25.6g Phase Speed beam   \n", dparam("WNB"));

  
  /* Write out the maps */
  
  fprintf (fp, "#\n# mapx \n#\n"); SPM_writeMap (fp, &mapx[i] );
  fprintf (fp, "#\n# mapy \n#\n"); SPM_writeMap (fp, &mapy[i] );
 
  fclose (fp);

 }//end of all partils
  
  return;
}

/* ------------------------------------------------------------------------ *
 * writeMap() - write a single map
 * ------------------------------------------------------------------------ */

void SPM_writeMap (FILE *fp, Mapping *map)
{
  int k, NZ = map->NZ;
  
  fprintf (fp, "%-25.6g\n", map->time);
  fprintf (fp, "%-25d  \n", map->NZ);

  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->d[k]);   fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->t[k]);   fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->tt[k]);  fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->f[k]);   fprintf (fp, "\n");

  return;
}


/* ------------------------------------------------------------------------ *
 * writedog2() - write displacement and forces and base pressure to session.dog
 * ------------------------------------------------------------------------ */

void C_SPM::writedog2 ()
{
  int pid, i, NZ = mapx->NZ;
  
  for(pid=0; pid<num_partls; ++pid)
   {
   for (i = 0; i < NZ; i++)
	  fprintf (Fp[pid], "%#10.5f  %#10.6f %#10.6f %#10.6f  %#10.6f %#2d\n",
	     mapx[pid].time, mapx[pid].d[i], mapy[pid].d[i], mapx[pid].f[i], mapy[pid].f[i], i);
    fflush (Fp[pid]);
  }

  return;
}

static int check_number = 0;
void C_SPM::WriteInit(char *name){

  int i,j;
  int  nztot = option("NZTOT");
  char fname[FILENAME_MAX];
  FILE *fp;
  
  if (option("SLICES")) {
    sprintf (fname, "%s_%d.init",  name,check_number);
    check_number++;
  }
  else
  {
    sprintf (fname, "%s.init",  name);
    int err = backup (fname);
    if (err != 0 ) 
      ErrorHandler ("Writeinit", "failed to backup the *.init file", WARNING);
  }
  
  if ((fp = fopen (fname,"w")) == (FILE *) NULL)
      ErrorHandler ("WriteInit", "failed to open the *.init file", ERROR);

 for(j=0; j<nztot; j++) 
  for(i=0; i<num_partls; i++) 
     fprintf(fp,"%#10.6f %#10.6f %#10.6f %#10.6f %#10.6f %#10.6f %#10.6f %#10.6f %#10.6f \n",Xp[j][i][0],Yp[j][i][0],ALFA[j][i][0],
                  Up[j][i][0],Vp[j][i][0],Ang[j][i][0],Fxext[j][i][0],Fyext[j][i][0], Text[j][i][0]);
 fclose(fp);
}
