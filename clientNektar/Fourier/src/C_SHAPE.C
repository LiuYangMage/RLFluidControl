#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nektarF.h"

#include <hotel.h>
#include "C_SHAPE.h"

#define max(a,b) ( (b) < (a) ? (a) : (b) )
#define min(a,b) ( (b) > (a) ? (a) : (b) )

static double indicator(double d, double ksi);
static double FindDist(const C_VERT &p, const C_VERT &v0, const C_VERT &v1);
static int    FindSign(const C_VERT &p, const C_VERT &v0, const C_VERT &v1, const C_VERT &v2, const C_VERT &v3);

static int    strake_num = 100;
static double **strake_profile = (double**)0;
static double xx0,yy0,xx1,yy1,xx2,yy2,xx3,yy3,kk;

void C_SHAPE::compute_indicator(int nqt, int plane, double *x, double *y, double *indicator_value){


/* MEMO

   for type "circle"
   data[0],data[1]     - center of the circle 
   data[2]             - radius  
   data[3]             - interface thickness 
   data[4],data[5]     - velocity of center of particle (data[4-5])
   data[6]             - mass ratio = fluid density / particle density

*/
  

//  if ( (strcmp(type,"cir") == 0) || (strcmp(type,"buoy") == 0) )
  if ( strcmp(type,"cir") == 0)
    {

    int offset = plane*9;
    double lx1, lx2, lx, radius; 
    int LoR, i;

    double Rc = data[offset+2];
    double thick = data[offset+3];

    for (i = 0; i < nqt; ++i){
      double dist=pow(pow(x[i]-data[offset+0],2.0)+pow( (y[i] - data[offset+1]),2.0),0.5);
      indicator_value[i]=(tanh((Rc-dist)/thick)+1.)/2.;
       
    } // end of "for (i = 0; i < Np; ++i)"

  }//end of "if (strcmp(type,"cir") == 0)" 
  else if ( strcmp(type,"buoy") == 0)
  {
    double lzz = dparam("LZ");
    int    nzz = option("NZTOT");
    double dzz = lzz/(double)nzz;
    double zz = plane*dzz;
    double scale = 100.;

    int offset = plane*9;
    double lx1, lx2, lx, radius; 
    int LoR, i,j;

    double Rc = data[offset+2];
    double thick = data[offset+3];

    for (i = 0; i < nqt; ++i){
      double dist=pow(pow(x[i]-data[offset+0],2.0)+pow( (y[i] - data[offset+1]),2.0),0.5);
      indicator_value[i]=(tanh((Rc-dist)/thick)+1.)/2.;
  /*
      double min_dist = 1e10;
      double rr = sqrt(x[i]*x[i]+y[i]*y[i]);

      for (j = 0; j < nzz; ++j){
         int offset1 = j*9;
         double dist= pow(pow(x[i] - data[offset1+0],2.0)
                         +pow(y[i] - data[offset1+1],2.0)
                         +pow(zz   - dzz*j,          2.0),0.5);

        if(dist <= min_dist)
         {
           min_dist = dist;
           Rc = data[offset1+2];
         }

        }

      indicator_value[i]=(tanh((Rc-min_dist)/thick)+1.)/2.;
*/
    } // end of "for (i = 0; i < Np; ++i)"

  }//end of "if (strcmp(type,"buoy") == 0)" 
  else if ( strcmp(type,"strake") == 0)
  {
    static int strake_step = 1;

    double lzz = dparam("LZ");
    int    nzz = option("NZTOT");
    double dzz = lzz/(double)nzz;
    double zz = plane*dzz;
    double scale = 100.;

    double dphi = 2.0*M_PI/(double)nzz*plane;

    int offset = plane*11;
    double lx1, lx2, lx, radius; 
    int LoR, i,j;

    double Rc = data[offset+2];
    double Hc = data[offset+3]; //strake height
    double Ac = data[offset+4]; //half of the strake angle 

	  if(strake_step == 1)
     {
      strake_profile = dmatrix(0, 1, 0, strake_num-1);
      
      xx0 = -0.5*sqrt(3.0)*(Rc+Hc); 
      yy0 = -0.5*(Rc+Hc);

      kk= tan(M_PI/6.0-Ac);

      double AA = 1.0+kk*kk;
      double BB = 2.0*kk*(yy0-kk*xx0);
      double CC = kk*kk*xx0*xx0+yy0*yy0-2.0*yy0*kk*xx0-Rc*Rc;

      xx1 = (-BB-sqrt(BB*BB-4.0*AA*CC))/(2.0*AA);
      yy1 = kk*(xx1-xx0)+yy0;

      xx2 = -xx1;
      yy2 = yy1;

      xx3 = -xx0;
      yy3 = yy0;

      double dx = (xx3-xx0)/(double) strake_num;

      for (i=0; i<strake_num; ++i)
       {
         double xx = xx0+i*dx;
         strake_profile[0][i] = xx;

         if((xx>=xx0)&&(xx<xx1))
          strake_profile[1][i] = yy0+kk*(xx-xx0);
         else if ((xx>xx1) && (xx<xx2))
           strake_profile[1][i] = -sqrt(Rc*Rc-xx*xx);
         else
          strake_profile[1][i] = yy3-kk*(xx-xx3);
       }

      strake_step ++;
     }

    double thick = data[offset+5];


    for (i = 0; i < nqt; ++i)
     {
       double xx = x[i]*cos(dphi)-y[i]*sin(dphi);
       double yy = y[i]*cos(dphi)+x[i]*sin(dphi);
//       double xx = x[i]; double yy = y[i];

       double xx_prime = xx;
       double yy_prime = yy;
       double rotating_angle = 2./3.*M_PI;

       double theta = atan2(yy,xx);

       if( (theta>-M_PI/6.) && (theta<1./2.*M_PI)) //rotate -2./3.*pi
        {
          xx_prime = xx*cos(-rotating_angle)-yy*sin(-rotating_angle);
          yy_prime = yy*cos(-rotating_angle)+xx*sin(-rotating_angle);
        }
       else if((theta>=M_PI/2.0)||(theta<-5.0*M_PI/6.0)) //rotate 2./3.*pi
        {
          xx_prime = xx*cos(rotating_angle)-yy*sin(rotating_angle);
          yy_prime = yy*cos(rotating_angle)+xx*sin(rotating_angle);
        }
      else
        {
          xx_prime = xx;
          yy_prime = yy;
        }
       
        double min_dist = 1e6;
       for(int j = 0; j<strake_num; ++j)
        {
         double xx_s = strake_profile[0][j];
         double yy_s = strake_profile[1][j];
         double dist_tmp = sqrt((xx_prime-xx_s)*(xx_prime-xx_s)+(yy_prime-yy_s)*(yy_prime-yy_s));

         min_dist = min(dist_tmp,min_dist);
        }
        
      double theta_1 = atan2(yy1,xx1);
      double theta_2 = atan2(yy2,xx2);

 //     fprintf(stderr,"theta1 = %lf  theta2 = %lf \n", theta_1, theta_2);

      double theta_prime = atan2(yy_prime,xx_prime);
      

       if(theta_prime<theta_1) 
         {
//           double kk_prime = yy_prime/xx_prime; //xx_prime could not be zero...
//           double xx_intersect = kk_prime*(yy0-kk*xx0)/(kk_prime-kk);
//           double yy_intersect = kk_prime*xx_intersect;
//           double dist_prime = sqrt(xx_prime*xx_prime+yy_prime*yy_prime);
//           double dist_intersect = sqrt(xx_intersect*xx_intersect+yy_intersect*yy_intersect);
//
//           if(dist_prime <= dist_intersect)
//              min_dist = -min_dist;
          if (yy_prime>(yy0+kk*(xx_prime-xx0))) //above the line
              min_dist = -min_dist;
//          if( ( sqrt(xx_prime*xx_prime+yy_prime*yy_prime)-Rc)<0.) //inside the circle
//             min_dist = -min_dist;

         }
       else if((theta_prime>=theta_1)&&(theta_prime<theta_2))
         {
          if( ( sqrt(xx_prime*xx_prime+yy_prime*yy_prime)-Rc)<0.) //inside the circle
             min_dist = -min_dist;
         }
       else
         {
//           double kk_prime = yy_prime/xx_prime; //xx_prime could not be zero...
//           double xx_intersect = kk_prime*(yy3+kk*xx3)/(kk_prime+kk);
//           double yy_intersect = kk_prime*xx_intersect;
//
//           double dist_prime = sqrt(xx_prime*xx_prime+yy_prime*yy_prime);
//           double dist_intersect = sqrt(xx_intersect*xx_intersect+yy_intersect*yy_intersect);
//
//           if(dist_prime <= dist_intersect)
//              min_dist = -min_dist;
           if(yy_prime>(yy3-kk*(xx_prime-xx3))) //above the line
               min_dist = -min_dist;
//          if( ( sqrt(xx_prime*xx_prime+yy_prime*yy_prime)-Rc)<0.) //inside the circle
//             min_dist = -min_dist;
         }


      indicator_value[i]=(tanh((-min_dist)/thick)+1.)/2.;

    } // end of "for (i = 0; i < Np; ++i)"

  }//end of "if (strcmp(type,"strake") == 0)" 
  else if ( strcmp(type,"sph") == 0  )
    {
    
    double lzz = dparam("LZ");
    int    nzz = option("NZTOT");
    double dzz = lzz/(double)nzz;

    double zz = plane*dzz;

    int offset = plane*9;
    double lx1, lx2, lx, radius; 
    int LoR, i;

    double Rc = data[offset+3];
    double thick = data[offset+4];

    for (i = 0; i < nqt; ++i){
      double dist=pow(pow(x[i]-data[offset+0],2.0)+pow( (y[i] - data[offset+1]),2.0),0.5);
      indicator_value[i]=(tanh((Rc-dist)/thick)+1.)/2.;
    } // end of "for (i = 0; i < Np; ++i)"

  }//end of "if (strcmp(type,"strake") == 0)" 
  else if ( strcmp(type,"sph") == 0  )
    {
    
    double lzz = dparam("LZ");
    int    nzz = option("NZTOT");
    double dzz = lzz/(double)nzz;

    double zz = plane*dzz;

    int offset = plane*9;
    double lx1, lx2, lx, radius; 
    int LoR, i;

    double Rc = data[offset+3];
    double thick = data[offset+4];

    for (i = 0; i < nqt; ++i){
      double dist=pow(pow(x[i] - data[offset+0],2.0)
                     +pow(y[i] - data[offset+1],2.0)
                     +pow(zz   - data[offset+2],2.0), 0.5);
      indicator_value[i]=(tanh((Rc-dist)/thick)+1.)/2.;

    } // end of "for (i = 0; i < Np; ++i)"

  }//end of "if (strcmp(type,"cir") == 0)" 
  else if (strcmp(type,"rect") == 0){

    int i;
    int offset = plane*12;
    double thick = data[offset+8];
    C_VERT v0,v1,v2,v3; 
    C_VERT p; 

    v0.x = data[offset+0];
    v0.y = data[offset+1];
    v1.x = data[offset+2];
    v1.y = data[offset+3];
    v2.x = data[offset+4];
    v2.y = data[offset+5];
    v3.x = data[offset+6];
    v3.y = data[offset+7];


    for (i = 0; i < nqt; ++i){
       p.x = x[i];
       p.y = y[i];

       double d0 = FindDist(p,v0,v1);
       double d1 = FindDist(p,v1,v2);
       double d2 = FindDist(p,v2,v3);
       double d3 = FindDist(p,v3,v0);

       d0 = min(d0,d1);
       d1 = min(d2,d3);

       double dist=min(d0,d1);

      int csign = FindSign(p,v0,v1,v2,v3); 

      dist *= (double) csign;

      indicator_value[i]=(tanh(dist/thick)+1.)/2.;

      if(indicator_value[i]>1.05)
        indicator_value[i] = 1.;
      if(indicator_value[i]<-0.05)
        indicator_value[i] = 0.;

    } // end of "for (i = 0; i < Np; ++i)"


  }//end of "if (strcmp(type,"cir") == 0)" 
  else if (strcmp(type,"cable2d") == 0){

    int i;
    int offset = plane*6;
    double thick = data[offset+2];
    
    double width = data[offset+1]; 
    double center_line = data[offset+0]; 

    for (i = 0; i < nqt; ++i){
      double dist= fabs(x[i]-center_line);
      indicator_value[i]=(tanh((width-dist)/thick)+1.)/2.;

    } // end of "for (i = 0; i < Np; ++i)"


  }//end of "if (strcmp(type,"cir") == 0)" 
  else if ( strcmp(type,"wall") == 0)
  {
    double lzz = dparam("LZ");
    int    nzz = option("NZTOT");
    double dzz = lzz/(double)nzz;
    double zz = plane*dzz;

    int offset = plane*10;
    double lx1, lx2, lx, radius; 
    int LoR, i,j;

    double z_position = data[offset+2];
    double wall_thick = data[offset+3];
    double thick = data[offset+4];

    double dist=fabs(zz-z_position);
    
    fprintf(stderr," dist = %g  z_position = %g  zz = %g \n",dist,z_position,zz);

    for (i = 0; i < nqt; ++i){
      indicator_value[i]=(tanh((wall_thick-dist)/thick)+1.)/2.;
    } // end of "for (i = 0; i < Np; ++i)"

  }//end of "if (strcmp(type,"wall") == 0)" 
  else
  {
    fprintf(stderr,"compute_indicator() Not implemeted yet !");
    exit(-1);
  }
}


void C_SHAPE::update_position(double dt){
   if ((strcmp(type,"cyl") == 0) || (strcmp(type,"buoy") == 0)) { 
/*
     Xp[plane] += dt*Up[plane];  
     data[4] += dt*Vp[plane];
     data[5] += dt*data[11];
     //update bounding box
     bound_box[0] = data[3]-5.0*data[8];
     bound_box[1] = data[3] + data[6] + 5.0*data[8];
     bound_box[2] = data[4] + data[7] + 5.0*data[8];
     bound_box[3] = data[4] - data[7] - 5.0*data[8];
     bound_box[4] = data[5] + data[7] + 5.0*data[8];
     bound_box[5] = data[5] - data[7] - 5.0*data[8];
*/
   }
   else{
    fprintf(stderr,"update_position() Not implemeted yet !");
    exit(-1);

   }
}

void C_SHAPE::update_position_and_velocity(double fx, double fy, double fz, double dt){

  if (strcmp(type,"cyl") == 0) {
    //LG: do not allow motion in the direction of the axis!
    //LG:  hardcoded value of cylinder volume (0.2*pi*r^2)
    data[9] = 0.0;

    double u_half, dt_inv_mass_ratio_volume, xp = data[4]; 
    double spring_constant;

    dt_inv_mass_ratio_volume = 1.0/(data[12]*0.2*M_PI*data[7]*data[7]); 
    spring_constant = 0.0*(0.215*2.0*M_PI)*(0.215*2.0*M_PI);    
    
    u_half =  data[10] + 0.5*dt*(-spring_constant*xp + fy*dt_inv_mass_ratio_volume);
    data[4] += dt*u_half;    
    data[10] = data[10] + dt*(-spring_constant*0.5*(xp+data[4]) + fy*dt_inv_mass_ratio_volume);

    data[11] = 0.0; 
  }
  else if (strcmp(type,"sph") == 0){
    double xp, u_half, dt_inv_mass_ratio_volume, spring_constant;
    dt_inv_mass_ratio_volume = 1.0/(data[8]*4.0/3.0*M_PI*data[4]*data[4]*data[4]);
    spring_constant = (0.215*2.0*M_PI)*(0.215*2.0*M_PI);
   
    xp = data[0];
    u_half =  data[5] + 0.5*dt*(-spring_constant*xp + fx*dt_inv_mass_ratio_volume);
    data[0] += dt*u_half;
    data[5] = data[5] + dt*(-spring_constant*0.5*(xp+data[0]) + fx*dt_inv_mass_ratio_volume);

    xp = data[1];
    u_half =  data[6] + 0.5*dt*(-spring_constant*xp + fy*dt_inv_mass_ratio_volume);
    data[1] += dt*u_half;
    data[6] = data[6] + dt*(-spring_constant * 0.5*(xp+data[1]) + fy*dt_inv_mass_ratio_volume);

    xp = data[2];
    u_half =  data[6] + 0.5*dt*(-spring_constant*xp + fz*dt_inv_mass_ratio_volume);
    data[2] += dt*u_half;
    data[6] = data[6] + dt*(-spring_constant * 0.5*(xp+data[2]) + fz*dt_inv_mass_ratio_volume);
  }

}


void C_SHAPE::update_velocity(double f_x, double f_y, double f_z, double dt){


  if (strcmp(type,"cyl") == 0) {
    //LG: do not allow motion in the direction of the axis!
    //LG:  hardcoded value of cylinder volume (0.2*pi*r^2)
    data[9]  = 0.0;
    double dt_inv_mass_ratio = dt/data[12]/0.157079633; 
    data[10] += f_y*dt_inv_mass_ratio;
    data[11] = 0.0;
  }
  else if (strcmp(type,"sph") == 0){
    double dt_inv_mass_ratio = dt/data[8];
    data[5] += f_x*dt_inv_mass_ratio;
    data[6] += f_y*dt_inv_mass_ratio;
    data[7] += f_z*dt_inv_mass_ratio;
  }
}
void C_SHAPE::get_position(double *x, double *y, double *z){
   if (strcmp(type,"cyl") == 0) {
     x[0] = data[3];
     y[0] = data[4];
     z[0] = data[5];
   }
   else if (strcmp(type,"sph") == 0){
     x[0] = data[0];
     y[0] = data[1];
     z[0] = data[2];
   }
}


void C_SHAPE::get_velocity(double *u, double *v, double *w){
   if (strcmp(type,"cyl") == 0) {
     u[0] = data[9];
     v[0] = data[10];
     w[0] = data[11];
   }
   else if (strcmp(type,"sph") == 0){
     u[0] = data[5];
     v[0] = data[6];
     w[0] = data[7];
   }
}

static double indicator(double d, double ksi){
  //fprintf(stdout,"indicator: ksi = %f\n",ksi);
  return 0.5*(tanh(-d/ksi)+1.0);
  //return exp(-(d*d)/ksi); 
}

//
//http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
//
static double FindDist(const C_VERT &p, const C_VERT &v0, const C_VERT &v1){
  
  C_VERT p1;

  double len = sqrt((v0.x-v1.x)*(v0.x-v1.x)+(v0.y-v1.y)*(v0.y-v1.y));

  if(len == 0.0) return sqrt((p.x-v0.x)*(p.x-v0.x)+(p.y-v0.y)*(p.y-v0.y));


  double t = sqrt((p.x-v0.x)*(p.x-v0.x)+(v1.y-v0.y)*(v1.y-v0.y))/len;

  t = max(0,min(1,t));

  p1.x = v0.x+t*(v1.x-v0.x);
  p1.y = v0.y+t*(v1.y-v0.y);

  return sqrt((p.x-p1.x)*(p.x-p1.x)+(p.y-p1.y)*(p.y-p1.y));

 }
//
//
//http://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
static int FindSign(const C_VERT &p, const C_VERT &v0, const C_VERT &v1, const C_VERT &v2, const C_VERT &v3){
  C_VERT d0,d1,d2;

  d0.x = p.x-v3.x;
  d0.y = p.y-v3.y;

  d1.x = v0.x-v3.x;
  d1.y = v0.y-v3.y;

  d2.x = v2.x-v3.x;
  d2.y = v2.y-v3.y;
  
  double x0 = d0.x*d1.x+d0.y*d1.y;
  double x1 = d0.x*d2.x+d0.y*d2.y;

  double y0= d1.x*d1.x+d1.y*d1.y;
  double y2= d2.x*d2.x+d2.y*d2.y;

  if((x0>0.)&&(x0<y0)&&(x1>0.)&&(x1<y2))
    return 1;
  else
    return -1;

 }














