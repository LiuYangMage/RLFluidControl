//#include "nektarF.h"


#ifndef NSE_UTIL_H
#define NSE_UTIL_H


void NSE_Iprod_div_Quad(Element *E,
			double *u, double *v);
void NSE_Iprod_div_Tri(Element *E,
		       double *u, double *v);
void NSE_Iprod_div(Element *E,
		   double *u, double *v);
void NSE_Iprod_div(Element_List *UL, 
		   double *u, double *v);
void FourierList_Iprod_div(Element_List *UL, double *u, double *v);
void FourierList_Iprod_div(Element_List *UL, 
			   double *u, double *v, double *w,
			   double *work);
void FourierList_Grad_Iprod_div(Element_List *EL,
				double *f,
				double *work_1, double *work_2,
				double *work_3);
void FourierList_Grad_Iprod_div(Element_List *EL,
				double *work_1, double *work_2,
				double *work_3);


void NSE_VarDensity_MakeFlux_div(double alpha, double *mu, 
				 Element *E, Bndry *B, double *u, double *v, int ifac);
void NSE_VarDensity_MakeFlux_div_A(double alpha, double *mu, 
				   Element *E, Bndry *B, double *u, double *v, int ifac);

void MakeFlux_div(Element *E, Bndry *B, double *u, double *v, int ifac);
void MakeFlux_div_Tri(Tri *E, Bndry *B, double *u, double *v, int ifac);
void MakeFlux_div_Quad(Quad *E, Bndry *B, double *u, double *v, int ifac);


void NSE_Imag_Multiply(double alpha, int nz, int htot, double *data, double *work);
void NSE_Imag_Multiply_BetaK(double alpha, int nz, int htot, double *data, double *work);

void FList_Calc_Vorticity(Element_List *U,
			  double *u, double *v, double *w,
			  double *vort_1, double *vort_2, double *vort_3,
			  double *work_1, double *work_2);


#endif // NSE_UTIL_H
