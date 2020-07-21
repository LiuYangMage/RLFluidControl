
#include <string.h>
#include "nektarF.h"
#include "dmath.h"
#include "nse_util.h"


//******************************************************************************//
//
void NSE_Iprod_div(Element_List *UL, 
		   double *u, double *v);

void FourierList_Iprod_div(Element_List *UL, double *u, double *v)
{
  // compute \int_{volume} [ V \cdot \grad(basis_p)],
  //   where V = (u,v), and basis_p is 2D basis function
  //
  // On input:
  //   UL : a Fourier_List structure, UL->(base_h,base_hj) to be used as workspace
  //   (u,v) contains V in (Q,F)-space, dimension: [UL->nz*UL->htot]
  //         Note: assume the velocity data is in Quadrature/Fourier space, not in Quadrature/Physical space
  //               therefore no de-aliasing for nz here
  //
  // On exit:
  //   UL->base_hj contains \int_{volume} [ V \cdot \grad(basis_p)], in J-space
  //     ->base_h contains junk
  //   (u,v) unchanged from before the call
  //

  int nz = UL->nz;
  int htot = UL->htot;
  int i;

  double *tu, *tv;

  for(i=0;i<nz;i++) {
    tu = u + i*htot;
    tv = v + i*htot;
    NSE_Iprod_div(UL->flevels[i], tu, tv);
  }

}



//******************************************************************************//
// Utility routines, variants of Iprod, 

void NSE_Iprod_div(Element *E,
		   double *u, double *v);

void NSE_Iprod_div(Element_List *UL, 
		   double *u, double *v)
{
  // compute \int (u*d_base_p/dx + v*d_base_p/dy + w*d_base_p/dz), where base_p
  //   is the basis function, for all base_p functions
  // result stored in UL->base_hj
  //
  // On input:
  //   (u,v) contains velocity, which cannot be NULL pointers
  //   UL->(base_h,base_hj) is workspace
  //
  // On exit:
  //   UL->base_hj contains results
  //   (u,v) unchanged
  //   UL->base_h contains garbage

  Element *E;
  double *tu = u,
    *tv = v;
  
  for(E=UL->fhead;E;E=E->next) {
    NSE_Iprod_div(E, tu, tv);
    tu = tu + E->qtot;
    tv = tv + E->qtot;
  }
}


void NSE_Iprod_div_Quad(Element *E,
			double *u, double *v);
void NSE_Iprod_div_Tri(Element *E,
		       double *u, double *v);

void NSE_Iprod_div(Element *E,
		   double *u, double *v)
{
  // elemental version of Iprod_div

  switch(E->identify()) {
  case Nek_Quad:
    NSE_Iprod_div_Quad(E, u, v);
    break;
  case Nek_Tri:
    NSE_Iprod_div_Tri(E, u, v);
    break;
  default:
    error_msg(NSE_Iprod_div(): unsupported element types);
    break;
  }

}


void NSE_Iprod_div_Quad(Element *Ef,
			double *u, double *v)
{
  // use E->h as workspace
  // On exit:
  //   E->vert->hj contains results
  //   (u,v) unchanged

  Quad *E = (Quad*)Ef;
  int qtot = E->qtot, Nm = E->Nmodes;
  Basis *b = E->getbasis(); // basis
  Basis *db = E->derbasis(); // derivative basis
  Geom *G = E->geom;
  int i;

  double *store = DMath::newD(Nm);

  // form_diprod
  // u*dr/dx + v*dr/dy  ---> E->h[0][*]
  for(i=0;i<qtot;i++)
    E->h[0][i] = u[i]*G->rx.p[i] + v[i]*G->ry.p[i];
  
  E->Iprod_d(E, db, b);
  dcopy(Nm, E->vert->hj, 1, store, 1); // save this part

  // u*ds/dx + v*ds/dy ---> E->h[0][*]
  for(i=0;i<qtot;i++)
    E->h[0][i] = u[i]*G->sx.p[i] + v[i]*G->sy.p[i];
  
  E->Iprod_d(E, b, db);
  daxpy(Nm, 1.0, store, 1, E->vert->hj, 1); // sum of the 2 parts
  // Now E->vert->hj[Nm] contains the results

  DMath::del(store);
}


void NSE_Iprod_div_Tri(Element *Ef,
		       double *u, double *v)
{
  // use E->h as workspace
  // On exit:
  //   E->vert->hj contains results
  //   (u,v) unchanged
  
  Tri *E = (Tri*) Ef;
  int qtot = E->qtot, Nm = E->Nmodes;
  Basis *b = E->getbasis(); // basis
  Basis *db = E->derbasis(); // derivative basis
  
  double *tv = DMath::newD(qtot+Nm);
  double *store = tv + qtot;

  dcopy(qtot, u, 1, E->h[0], 1);
  dcopy(qtot, v, 1, tv, 1);

  E->form_diprod(E->h[0], tv, NULL, b->vert+1); // from HelmHoltz routine
  // Now E->h[0] contains: u*drdx + v*dr/dy
  //     tv contains: u*ds/dx + v*ds/dy

  E->Iprod_d(E, db, b);
  dcopy(Nm, E->vert->hj, 1, store, 1); // save first part

  dcopy(qtot, tv, 1, E->h[0], 1);
  E->Iprod_d(E, b, db); // second part

  daxpy(Nm, 1.0, store, 1, E->vert->hj, 1); 
  // now E->vert->hj contains result


  // clean up
  DMath::del(tv);
}


//***************************************************************************//


void FourierList_Iprod_div(Element_List *UL, 
			   double *u, double *v, double *w,
			   double *work)
{
  // compute: 
  //        \int_{volume_2.5D} [ (u,v,w) \cdot grad{basis_p_k}]
  //    = \int_{volume_2D} [(u,v)_k \cdot grad_2D{basis_p}] - i*Beta(k)*\int_{volume_2D} [w_k * basis_p]
  // 
  // On input:
  //   (u,v,w) contains velocity in (Q,F)-space, note: data in Quadrature-Fourier space, 
  //           dimension: [UL->nz*UL->htot]
  //   UL is Element_List structure, workspace
  //   work is workspace, dimension: [UL->nz*UL->hjtot]
  //
  // On exit:
  //   UL->base_hj contains integral result in (J,F)-space
  //     ->base_h contains junk
  //   (u,v,w) unchanged
  //   work contains junk
  
  int nz = UL->nz,
    htot = UL->htot,
    hjtot = UL->hjtot;
  int htot_nz = nz*htot;
  int hjtot_nz = nz*hjtot;

  dcopy(htot_nz, w, 1, UL->base_h, 1);
  UL->Iprod(UL);
  // Now UL->base_hj contains \int_{volume_2D} [ w_k * basis_p], in (J,F)-space

  NSE_Imag_Multiply_BetaK(-1.0, UL->nz, UL->hjtot, UL->base_hj, work);
  // Now UL->base_hj contains -i*Beta(k)*\int_{volume_2D} [ w_k * basis_p], in (J,F)-space

  // save this integral to work
  dcopy(hjtot_nz, UL->base_hj, 1, work, 1); 
  // Now work contains -i*Beta(k)*\int_{volume_2D} [ w_k * basis_p], in (J,F)-space

  FourierList_Iprod_div(UL, u, v);
  // Now UL->base_hj contains \int_{volume_2D} [(u,v)_k \cdot \grad(basis_p)], in (J,F)-space

  // add two integrals
  daxpy(hjtot_nz, 1.0, work, 1, UL->base_hj, 1);

  // Now UL->base_hj contains 
  //        \int_{volume_2D} [(u,v)_k \cdot \grad(basis_p)] - i*Beta(k)*\int_{volume_2D} [ w_k * basis_p]
  //            in (J,F)-space
  //       ->base_h contains junk
  //     work contains nothing important
  //     (u,v,w) unchanged

}




//************************************************************************//
// Routines for computing surface integral with divergence

void NSE_VarDensity_MakeFlux_div(double alpha, double *mu, 
				 Element *E, Bndry *B, double *u, double *v, int ifac)
{
  // Compute 
  //    \int_{boundary} [alpha * mu * [V \cdot grad(basis_p)] ] 
  //   for all basis_p where V=(u,v) is a vector field, B stores boundary information,
  //   and result will be stored in E->vert[0].hj[Nmodes]. Element E
  //   must be the element with id the same as B->elmt->id, i.e., have
  //   the geometry that B corresponds to.
  //
  // On input:
  //   B contains the boundary information
  //   (u,v) contains vector field data on this boundary, dimension : [qa] 
  //         Note: (u,v) already interpolated to grids consistent with face "ifac"
  //   ifac is id of face whose grids are consistent with grids on which (u,v), mu data are given
  //   alpha is coefficient
  //   mu is a scalar field, dimension : [qa] , on grids consistent with face "ifac"
  //   E->vert->hj and E->h[qb][qa] contains nothing important, will be overwritten
  //
  // On exit:
  //   E->vert->hj[Nmodes] contains integral in J-space
  //   E->h[qb][qa] contains junk
  //   B unchanged
  //   (u,v) unchanged
  //   alpha and mu unchanged

  int qa = E->qa,
    qb = E->qb;
  int qm = qa > qb ? qa : qb;
  int i;

  double *u_tmp = DMath::newD(qm*2);
  double *v_tmp = &u_tmp[qm];
  
  for(i=0;i<qm;i++) {
    u_tmp[i] = alpha*mu[i]*u[i];
    v_tmp[i] = alpha*mu[i]*v[i];
  }

  MakeFlux_div(E, B, u_tmp, v_tmp, ifac);

  // clean up
  DMath::del(u_tmp);
}

//***********************************************************//

void NSE_VarDensity_MakeFlux_div_A(double alpha, double *mu, 
				   Element *E, Bndry *B, double *u, double *v, int ifac)
{
  // Compute 
  //    \int_{boundary} [alpha * mu * [V \cdot grad(basis_p)] ] 
  //   for all basis_p where V=(u,v) is a vector field, B stores boundary information,
  //   and result will be stored in E->vert[0].hj[Nmodes]. Element E
  //   must be the element with id the same as B->elmt->id, i.e., have
  //   the geometry that B corresponds to.
  //
  // Note:
  //   if pointer mu is NULL, will assume that mu does not appear in the above integral
  //   This is the difference from function NSE_VarDensity_MakeFlux_div()
  //
  // On input:
  //   B contains the boundary information
  //   (u,v) contains vector field data on this boundary, dimension : [qa] 
  //         Note: (u,v) already interpolated to grids consistent with face "ifac"
  //   ifac is id of face whose grids are consistent with grids on which (u,v), mu data are given
  //   alpha is coefficient
  //   mu is a scalar field, dimension : [qa] , on grids consistent with face "ifac"
  //   E->vert->hj and E->h[qb][qa] contains nothing important, will be overwritten
  //
  // On exit:
  //   E->vert->hj[Nmodes] contains integral in J-space
  //   E->h[qb][qa] contains junk
  //   B unchanged
  //   (u,v) unchanged
  //   alpha and mu unchanged

  int qa = E->qa,
    qb = E->qb;
  int qm = qa > qb ? qa : qb;
  int i;

  double *u_tmp = DMath::newD(qm*2);
  double *v_tmp = &u_tmp[qm];
  
  if(mu) {
    for(i=0;i<qm;i++) {
      u_tmp[i] = alpha*mu[i]*u[i];
      v_tmp[i] = alpha*mu[i]*v[i];
    }
  }
  else {
    dscal(qm, alpha, u, 1);
    dscal(qm, alpha, v, 1);
  }

  MakeFlux_div(E, B, u_tmp, v_tmp, ifac);

  // clean up
  DMath::del(u_tmp);
}


//******************************************************************//

void MakeFlux_div(Element *E, Bndry *B, double *u, double *v, int ifac)
{
  // Compute \int_{boundary} [ V \cdot grad(basis_p) ] for all basis_p
  //   where V=(u,v) is a vector field, B stores boundary information,
  //   and result will be stored in E->vert[0].hj[Nmodes]. Element E
  //   must be the element with id the same as B->elmt->id, i.e., have
  //   the geometry that B corresponds to.
  //
  // On input:
  //   B contains the boundary information
  //   (u,v) contains vector field data on this boundary, 
  //         already interpolated to grids consistent with face "ifac". This is important
  //   ifac is the id of the face whose grids are consistent with the grids on which (u,v) data are given!
  //   E->vert->hj and E->h[qb][qa] contains nothing important, will be overwritten
  //
  // On exit:
  //   E->vert->hj[Nmodes] contains integral in J-space
  //   E->h[qb][qa] contains junk
  //   B unchanged
  //   (u,v) unchanged
  //

  switch(E->identify()) {
  case Nek_Quad:
    MakeFlux_div_Quad((Quad*)E, B, u, v, ifac);
    break;
  case Nek_Tri:
    MakeFlux_div_Tri((Tri*)E, B, u, v, ifac);
    break;
  default:
    error_msg(MakeFlux_div(): unsupported element types);
    break;
  }

}


void MakeFlux_div_Quad(Quad *E, Bndry *B, double *u, double *v, int ifac)
{
  // Compute \int_{boundary} [ V \cdot grad(basis_p) ] for all basis_p
  //   where V=(u,v) is a vector field, B stores boundary information,
  //   and result will be stored in E->vert[0].hj[Nmodes]. Element E
  //   must be the element with id the same as B->elmt->id, i.e., have
  //   the geometry that B corresponds to.
  //
  // On input:
  //   B contains the boundary information
  //   (u,v) contains vector field data on this boundary, dimension: [qa]
  //         Note: (u,v) data already interpolated to grids consistent with the grids of face "ifac"!
  //   ifac is id of face whose grids are consistent with grids on which (u,v) data are given
  //   E->vert->hj and E->h[qb][qa] contains nothing important, will be overwritten
  //
  // On exit:
  //   E->vert->hj[Nmodes] contains integral in J-space
  //   E->h[qb][qa] contains junk
  //   B unchanged
  //   (u,v) unchanged
  //
  // Note:
  //   V \cdot grad(basis_p) =  [V\cdot grad(r)]*d(basis_p)/dr
  //                           +[V\cdot grad(s)]*d(basis_p)/ds
  //

  int qa = E->qa,
    qb = E->qb,
    qtot = E->qtot;
  int qm = qa>qb ? qa : qb;
  int lmax = E->lmax;
  int face = B->face;
  size_t d = sizeof(double);
  int i;

  Edge *e = E->edge;
  Face *f = E->face;

  Geom *G = E->geom;
  Basis *bb = E->getbasis(); // basis
  Basis *db = E->derbasis(); // derivative basis
  double *phi_b, *phi_a, *phi_db, *phi_da;

  double *ui = DMath::newD(2*qa);
  double *vi = &ui[qa];
  
  double *rx = DMath::newD(5*qm);
  double *ry = &rx[qm],
    *sx = &rx[2*qm],
    *sy = &rx[3*qm],
    *tmp = &rx[4*qm];

  double *H_data, *wk, *vec_wk_1, *vec_wk_2, *H; // workspaces
  H_data = DMath::newD(lmax*lmax + lmax*2); 
  vec_wk_1 = &H_data[lmax*lmax];
  vec_wk_2 = &H_data[lmax*lmax+lmax];
  wk = DMath::newD(lmax*qm);

  double *wa, *sjac;
  getzw(qa, &wa, &wa, 'a'); // weight data
  sjac = B->sjac.p; // surface jacobian for this boundary

  // first get (u,v) on face 0
  // Note (u,v) data already on grids consistent with face "ifac"!!
  E->InterpToFace1(ifac, u, ui); // Note: ifac in the parameter, should not be "face" here
  E->InterpToFace1(ifac, v, vi);
  // Now (ui,vi)[0:qa-1] contains (u,v) data on grid points consistent with face 0

  // get grad(r) and grad(s) data
  E->GetFace(G->rx.p, face, tmp); 
  E->InterpToFace1(face, tmp, rx);
  E->GetFace(G->ry.p, face, tmp);
  E->InterpToFace1(face, tmp, ry);
  E->GetFace(G->sx.p, face, tmp);
  E->InterpToFace1(face, tmp, sx);
  E->GetFace(G->sy.p, face, tmp);
  E->InterpToFace1(face, tmp, sy);
  // Now (rx,ry,sx,sy)[0:qa-1] contains gradient(r) and gradient(s) boundary data on
  //     grid points consistent with face 0
  
  //****************************************************************//
  // Compute \int_{boundary} [V\cdot grad(r)]*d(basis_p)/dr

  for(i=0;i<qa;i++)
    tmp[i] = (ui[i]*rx[i] + vi[i]*ry[i])*sjac[i]*wa[i];
  // Now tmp[0:qa-1] contains [V \cdot grad(r)]*Jacobian*Weight

  // compute d(basis_p)/dr = phi_a^{prime}(r_i)*phi_b(s_j)
  switch(face) {

    //----------------------------------------------------//
    // face 0
  case 0:
    // for face 0, phi_b(s_j) = phi_b(s_j=-1)
    phi_b = bb->vert[2].b; // Now phi_b[lmax*qb]
    dcopy(lmax, phi_b, qb, vec_wk_1, 1);
    // Now vec_wk_1[0:lmax-1] contains phi_b(s_j=-1)
    
    phi_da = db->vert[1].a; // Now phi_da[lmax*qa] contains derivative basis
    for(i=0;i<lmax;i++) 
      E->InterpToFace1(face, &phi_da[i*qa], &wk[i*qa]);
    // Now wk[lmax*qa] contains derivative basis in a direction on grid points consistent with face 1
    
    // compute phi_da(xi_i)*F(xi_i), where F(xi_i) is contained in tmp[0:qa-1]
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains \sum_i [phi_da[lmax][i] * F[i] ]

    // Now that vec_wk_1 contains phi_b[lmax] and vec_wk_2 contains phi_da[lmax]
    // compute the matrix H_data[b][a]_{lmax x lmax}, Note: first index is b, second index is a
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p

    break;

    //---------------------------------------------------//
    // face 1
  case 1:

    // for face 1, phi_da(r_i) = phi_da(r_i=1)
    phi_da = db->vert[1].a; // Now phi_da[lmax*qa]
    dcopy(lmax, &phi_da[qa-1], qa, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains phi_da[0:lmax-1][qa-1], i.e. phi_da(r_i=1)

    phi_b = bb->vert[2].b; // Now phi_b[lmax*qb]
    // need to interpolate to grid points consistent with face 1
    for(i=0;i<lmax;i++)
      E->InterpToFace1(face, &phi_b[i*qb], &wk[i*qa]);
    // Now wk[lmax*qa] contains basis in b direction on grid points consistent with face 1

    // compute phi_b(xi_i)*F(xi_i), where tmp[0:qa-1] contains F(xi_i)
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_1, 1);
    // Now vec_wk_1[0:lmax-1] contains \sum_i [ phi_b[lmax][i] * F[i] ]

    // Now that vec_wk_1 contains phi_b[lmax] and vec_wk_1 contains phi_da[lmax]
    // compute matrix H_data[lmax][lmax]
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p
    
    break;

    //---------------------------------------------------//
    // face 2
  case 2:

    // for face 2, phi_b(s_j) = phi_b(s_j=1)
    phi_b = bb->vert[2].b; // Now phi_b[lmax*qb]
    dcopy(lmax, &phi_b[qb-1], qb, vec_wk_1, 1);
    // Now vec_wk_1[0:lmax-1] contains phi_b(s_j=1)

    phi_da = db->vert[1].a; // Now phi_da[lmax*qa] contains derivative basis
    for(i=0;i<lmax;i++) 
      E->InterpToFace1(face, &phi_da[i*qa], &wk[i*qa]);
    // Now wk[lmax*qa] contains derivative basis in a direction on grid points consistent with face 1
    
    // compute phi_da(xi_i)*F(xi_i), where F(xi_i) is contained in tmp[0:qa-1]
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains \sum_i [phi_da[lmax][i] * F[i] ]

    // Now that vec_wk_1 contains phi_b[lmax] and vec_wk_2 contains phi_da[lmax]
    // compute the matrix H_data[b][a]_{lmax x lmax}, Note: first index is b, second index is a
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p
    
    break;

    //---------------------------------------------------//
    // face 3
  case 3:

    // for face 3, phi_da(r_i) = phi_da(r_i=-1)
    phi_da = db->vert[1].a; // Now phi_da[lmax*qa]
    dcopy(lmax, phi_da, qa, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains phi_da[0:lmax-1][0], i.e. phi_da(r_i=-1)

    phi_b = bb->vert[2].b; // Now phi_b[lmax*qb]
    // need to interpolate to grid points consistent with face 1
    for(i=0;i<lmax;i++)
      E->InterpToFace1(face, &phi_b[i*qb], &wk[i*qa]);
    // Now wk[lmax*qa] contains basis in b direction on grid points consistent with face 1

    // compute phi_b(xi_i)*F(xi_i), where tmp[0:qa-1] contains F(xi_i)
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_1, 1);
    // Now vec_wk_1[0:lmax-1] contains \sum_i [ phi_b[lmax][i] * F[i] ]

    // Now that vec_wk_1 contains phi_b[lmax] and vec_wk_1 contains phi_da[lmax]
    // compute matrix H_data[lmax][lmax]
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p

    break;

  } // switch(face)


  // Unpack data in H
  // copy field into E->hj for moment 
  H = H_data;
  E->vert[2].hj[0] = *H; ++H;
  E->vert[3].hj[0] = *H; ++H; 
  memcpy(e[2].hj,H,e[2].l*d); H += lmax-2;
  E->vert[1].hj[0] = *H; ++H;
  E->vert[0].hj[0] = *H; ++H;
  memcpy(e[0].hj,H,e[0].l*d); H += lmax-2;
   
  dcopy(e[1].l, H, lmax, e[1].hj, 1);++H;
  dcopy(e[3].l, H, lmax, e[3].hj, 1);++H;
  
  for(i = 0; i < f[0].l; ++H, ++i)
    dcopy(f[0].l, H, lmax, f[0].hj[i], 1);

  // Now E->vert->hj[Nmodes] contains integral \int_{boundary} [V\cdot grad(r)]*d(basis_p)/dr
  //    in J-space

  double *save = DMath::newD(E->Nmodes);
  dcopy(E->Nmodes, E->vert[0].hj, 1, save, 1);
  // Now save[0:Nmodes-1] contains resulting integral \int_{boundary} [V\cdot grad(r)]*d(basis_p)/dr


  //*************************************************************//
  // Compute \int_{boundary} [V\cdot grad(s)]*d(basis_p)/ds

  for(i=0;i<qa;i++)
    tmp[i] = (ui[i]*sx[i] + vi[i]*sy[i])*sjac[i]*wa[i];
  // Now tmp[0:qa-1] contains [V \cdot grad(s)]*Jacobian*Weight, or F(xi_i)

  // compute d(basis_p)/ds = phi_a(r_i)*phi_b^{prime}(s_j)
  switch(face) {

    //-----------------------------------------------------//
  case 0:
    
    // for face 0, phi_db(s_j) = phi_db(s_j=-1)
    phi_db = db->vert[2].b; // Now phi_db[lmax*qb], derivative basis
    dcopy(lmax, phi_db, qb, vec_wk_1, 1);
    // Now vec_wk_1[0:lmax-1] contains phi_db(s_j=-1)

    phi_a = bb->vert[1].a; // Now phi_a[lmax*qa]
    for(i=0;i<lmax;i++)
      E->InterpToFace1(face, &phi_a[i*qa], &wk[i*qa]);
    // Now wk[lmax*qa] contains basis phi_a in a direction on grids consistent with face 1

    // compute phi_a(xi_i)*F(xi_i), where tmp[0:qa-1] contains F(xi_i)
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains \sum_i [ phi_a[lmax][i]*F[i] ]

    // Now that vec_wk_1 contains phi_db[lmax] and vec_wk_2 contains phi_a[lmax]
    // compute the matrix H_data[b][a]_{lmax x lmax}, Note: firs index is b , second index is a
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p    

    break;

    //-----------------------------------------------------//
  case 1:

    // for face 1, phi_a(r_i) = phi_a(r_i=1)
    phi_a = bb->vert[1].a; // phi_a[lmax*qa]
    dcopy(lmax, &phi_a[qa-1], qa, vec_wk_2, 1);
    // Now vec_wk_2 [0:lmax-1] contains phi_a[0:lmax-1][qa-1], i.e. phi_a(r_i=1)

    phi_db = db->vert[2].b; // phi_db[lmax*qb]
    // need to interpolate to grids consistent with face 1
    for(i=0;i<lmax;i++)
      E->InterpToFace1(face, &phi_db[i*qb], &wk[i*qa]);
    // Now wk[lmax*qa] contains derivative basis on grids consistent with face 1

    // compute phi_db(xi_i)*F(xi_i)
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_1, 1);
    // Now vec_wk_1[lmax] contains \sum_i [phi_db[lmax][i] * F[i] ]

    // Now that vec_wk_1 contains phi_db[lmax] and vec_wk_2 contains phi_a[lmax]
    // compute matrix H_data[b][a]_{lmax x lmax}, Note: firs index is b , second index is a
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p    
    
    break;

    //------------------------------------------------------//
  case 2:

    // for face 2, phi_db(s_j) = phi_db(s_j=1)
    phi_db = db->vert[2].b; // Now phi_db[lmax*qb]
    dcopy(lmax, &phi_db[qb-1], qb, vec_wk_1, 1);
    // Now vec_wk_1[0:lmax-1] contains phi_db(s_j=1)

    phi_a = bb->vert[1].a; // Now phi_a[lmax*qa]
    for(i=0;i<lmax;i++)
      E->InterpToFace1(face, &phi_a[i*qa], &wk[i*qa]);
    // Now wk[lmax*qa] contains basis phi_a in a direction on grids consistent with face 1

    // compute phi_a(xi_i)*F(xi_i), where tmp[0:qa-1] contains F(xi_i)
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains \sum_i [ phi_a[lmax][i]*F[i] ]

    // Now that vec_wk_1 contains phi_db[lmax] and vec_wk_2 contains phi_a[lmax]
    // compute the matrix H_data[b][a]_{lmax x lmax}, Note: firs index is b , second index is a
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p    
    
    break;

    //------------------------------------------------------//
  case 3:

    // for face 3, phi_a(r_i) = phi_a(r_i=-1)
    phi_a = bb->vert[1].a; // phi_a[lmax*qa]
    dcopy(lmax, phi_a, qa, vec_wk_2, 1);
    // Now vec_wk_2[0:lmax-1] contains phi_a[0:lmax-1][0], i.e. phi_a(r_i=-1)

    phi_db = db->vert[2].b; // phi_db[lmax*qb]
    // need to interpolate to grids consistent with face 1
    for(i=0;i<lmax;i++)
      E->InterpToFace1(face, &phi_db[i*qb], &wk[i*qa]);
    // Now wk[lmax*qa] contains derivative basis on grids consistent with face 1

    // compute phi_db(xi_i)*F(xi_i)
    dgemv('T', qa, lmax, 1.0, wk, qa, tmp, 1, 0.0, vec_wk_1, 1);
    // Now vec_wk_1[lmax] contains \sum_i [phi_db[lmax][i] * F[i] ]

    // Now that vec_wk_1 contains phi_db[lmax] and vec_wk_2 contains phi_a[lmax]
    // compute matrix H_data[b][a]_{lmax x lmax}, Note: firs index is b , second index is a
    for(i=0;i<lmax;i++) {
      dcopy(lmax, vec_wk_2, 1, &H_data[i*lmax], 1);
      dscal(lmax, vec_wk_1[i], &H_data[i*lmax], 1);
    }
    // Now H_data[lmax][lmax] contains integral in J-space for all basis_p    
 
    break;

  } // switch(face)


  // Unpack data in H
  // copy field into E->hj for moment 
  H = H_data;
  E->vert[2].hj[0] = *H; ++H;
  E->vert[3].hj[0] = *H; ++H; 
  memcpy(e[2].hj,H,e[2].l*d); H += lmax-2;
  E->vert[1].hj[0] = *H; ++H;
  E->vert[0].hj[0] = *H; ++H;
  memcpy(e[0].hj,H,e[0].l*d); H += lmax-2;
   
  dcopy(e[1].l, H, lmax, e[1].hj, 1);++H;
  dcopy(e[3].l, H, lmax, e[3].hj, 1);++H;
  
  for(i = 0; i < f[0].l; ++H, ++i)
    dcopy(f[0].l, H, lmax, f[0].hj[i], 1);

  // Now E->vert->hj[Nmodes] contains integral \int_{boundary} [V\cdot grad(s)]*d(basis_p)/ds
  //    in J-space
  // Note:
  //   save[0:Nmodes-1] contains \int_{boundary} [V\cdot grad(r)]*d(basis_p)/dr

  //*************************************************************//
  // summ up two integrals

  daxpy(E->Nmodes, 1.0, save, 1, E->vert[0].hj, 1);
  // Now E->vert->hj[0:Nmodes-1] contains
  //   \int_{boundary} [V\cdot grad(r)]*d(basis_p)/dr + \int_{boundary} [V\cdot grad(s)]*d(basis_p)/ds
  // or \int_{boundary} [V \cdot grad(basis_p)] for all basis_p


  //*************************************************************//
  // clean up
  DMath::del(ui);
  DMath::del(rx);
  DMath::del(H_data);
  DMath::del(wk);
  DMath::del(save);
}






void MakeFlux_div_Tri(Tri *E, Bndry *B, double *u, double *v, int ifac)
{
  // Compute \int_{boundary} [ V \cdot grad(basis_p) ] for all basis_p
  //   where V=(u,v) is a vector field, B stores boundary information,
  //   and result will be stored in E->vert[0].hj[Nmodes]. Element E
  //   must be the element with id the same as B->elmt->id, i.e., have
  //   the geometry that B corresponds to.
  //
  // On input:
  //   B contains the boundary information
  //   (u,v) contains vector field data on this boundary
  //   E->vert->hj and E->h[qb][qa] contains nothing important, will be overwritten
  //
  // On exit:
  //   E->vert->hj[Nmodes] contains integral in J-space
  //   E->h[qb][qa] contains junk
  //   B unchanged
  //   (u,v) unchanged
  //

  error_msg(MakeFlux_div(): triangular element not supported);

}


//********************************************************************************************//

void NSE_Imag_Multiply(double alpha, int nz, int htot, double *data, double *work)
{
  // compute quantity
  //     i * alpha * data[:]
  // where i = sqrt(-1), alpha is a real constant, data[:] is complex data in Fourier-space
  //    and stored plane by plane
  //
  // On input:
  //  alpha : coefficient
  //  nz, htot : dimensions
  //  data : complex data in Fourier space, dimension [nz][htot], assumed to be stored in 
  //         following format:
  //           data[0][0:htot-1], if parid(0) == 0, then this is Fourier mode 0, only real part
  //                              otherwise, this is the real part of a normal Fourier mode
  //               [1][0:htot-1], if parid(1) == 1, then this Fourier mode Nztot/2, only real part
  //                              otherwise, this is the imaginary part of a normal Fourier mode
  //               [2][0:htot-1], real part of Fourier mode parid(2)/2
  //               [3][0:htot-1], imaginary part of Fourier mode parid(2)/2
  //               ... etc etc
  // work : worspace, dimension of at least [htot]
  //
  // On exit:
  //  data[*] updated and stored result: i*alpha*data[*]
  //  work : contains junk
  // 
  
  int k, k_glob;

  for(k=0;k<nz;k=k+2) {

    k_glob = parid(k);

    if(k_glob<2) {
      // this is Fourier mode 0 and mode Nztot/2, no imaginary part, only real part
      memset(&data[k*htot], '\0', sizeof(double)*htot*2);
      continue;
    }

    // Now this is a normal mode
    // Now data[k*htot:(k+1)*htot-1] contains real part of this Fourier mode
    //     data[(k+1)*htot:(k+2)*htot-1] contains imaginary part of this Fourier mode

    memset(work, '\0', sizeof(double)*htot);
    daxpy(htot, -alpha, &data[(k+1)*htot], 1, work, 1);
    // Now work contains the real part of i*alpha*data[:]

    memset(&data[(k+1)*htot], '\0', sizeof(double)*htot);
    daxpy(htot, alpha, &data[k*htot], 1, &data[(k+1)*htot], 1);
    // Now data[(k+1)*htot:(k+2)*htot-1] contains the imaginary part of i*alpha*data[:]

    dcopy(htot, work, 1, &data[k*htot], 1);
    // Now data[k*htot:(k+1)*htot-1] contains the real part of i*alpha*data[:]

  } // for(k=0 ...


}



void NSE_Imag_Multiply_BetaK(double alpha, int nz, int htot, double *data, double *work)
{
  // compute quantity
  //     i * alpha * Beta(k) * data[:]
  // where i = sqrt(-1), alpha is a real constant, data[:] is complex data in Fourier-space
  //    and stored plane by plane
  //
  // On input:
  //  alpha : coefficient
  //  nz, htot : dimensions
  //  data : complex data in Fourier space, dimension [nz][htot], assumed to be stored in 
  //         following format:
  //           data[0][0:htot-1], if parid(0) == 0, then this is Fourier mode 0, only real part
  //                              otherwise, this is the real part of a normal Fourier mode
  //               [1][0:htot-1], if parid(1) == 1, then this Fourier mode Nztot/2, only real part
  //                              otherwise, this is the imaginary part of a normal Fourier mode
  //               [2][0:htot-1], real part of Fourier mode parid(2)/2
  //               [3][0:htot-1], imaginary part of Fourier mode parid(2)/2
  //               ... etc etc
  // work : worspace, dimension of at least [htot]
  //
  // On exit:
  //  data[*] updated and stored result: i*alpha*data[*]
  //  work : contains junk
  // 
  
  double beta_k;
  double coeff;
  int k, k_glob;

  for(k=0;k<nz;k=k+2) {

    k_glob = parid(k);

    if(k_glob<2) {
      // this is Fourier mode 0 and mode Nztot/2, no imaginary part, only real part
      memset(&data[k*htot], '\0', sizeof(double)*htot*2);
      continue;
    }

    // Now this is a normal mode
    // Now data[k*htot:(k+1)*htot-1] contains real part of this Fourier mode
    //     data[(k+1)*htot:(k+2)*htot-1] contains imaginary part of this Fourier mode

    beta_k = Beta(k);
    coeff = alpha*beta_k;

    memset(work, '\0', sizeof(double)*htot);
    daxpy(htot, -coeff, &data[(k+1)*htot], 1, work, 1);
    // Now work contains the real part of i*alpha*data[:]

    memset(&data[(k+1)*htot], '\0', sizeof(double)*htot);
    daxpy(htot, coeff, &data[k*htot], 1, &data[(k+1)*htot], 1);
    // Now data[(k+1)*htot:(k+2)*htot-1] contains the imaginary part of i*alpha*data[:]

    dcopy(htot, work, 1, &data[k*htot], 1);
    // Now data[k*htot:(k+1)*htot-1] contains the real part of i*alpha*data[:]

  } // for(k=0 ...


}




//********************************************************************************************//

void FList_Calc_Vorticity(Element_List *U,
			  double *u, double *v, double *w,
			  double *vort_1, double *vort_2, double *vort_3,
			  double *work_1, double *work_2)
{
  // compute vorticity in (Q,F)-space
  //
  // On input:
  //   U : provides Element_list structure
  //   (u,v,w): contains velocity in (Q,F)-space, dimension: [U->nz*U->htot]
  //   (vort_1,vort_2,vort_3): to contains vorticity in (Q,F)-space on return, dimension: [U->nz*U->htot]
  //   work_1, work_2: workspace, dimension: [U->nz*U->htot]
  //
  // On exit:
  //   vort_1, vort_2, vort_3: contains vorticity in (Q,F)-space
  //   (u,v,w) unchanged
  //   Uf unchanged
  //   work_1, work_2 contains junk
  //

  int htot = U->htot,
    nz = U->nz;
  int htot_nz = htot*nz;

  U->Grad_z_h(nz, htot, u, vort_2); // vort_2 <--- du/dz in (Q,F)-space

  U->Grad_z_h(nz, htot, v, vort_1); // vort_1 <--- dv/dz in (Q,F)-space
  dscal(htot_nz, -1.0, vort_1, 1);
  // Now vort_1 contains -dv/dz in (Q,F)-space

  U->Grad_h(w, work_1, work_2, 0, 'a');
  // Now work_1 contains dw/dx in (Q,F)-space
  //     work_2 contains dw/dy in (Q,F)-space

  daxpy(htot_nz, 1.0, work_2, 1, vort_1, 1);
  // Now vort_1 contains x component of vorticity in (Q,F)-space: dw/dy - dv/dz
  
  daxpy(htot_nz, -1.0, work_1, 1, vort_2, 1);
  // Now vort_2 contains y component of vorticity in (Q,F)-space: du/dz - dw/dx

  U->Grad_h(u, 0, work_1, 0, 'y'); // work_1 <--- du/dy, in (Q,F)-space
  U->Grad_h(v, vort_3, 0, 0, 'x'); // vort_3 <--- dv/dx in (Q,F)-space
  daxpy(htot_nz, -1.0, work_1, 1, vort_3, 1);
  // Now vort_3 contains z component of vorticity in (Q,F)-space: dv/dx - du/dy


}





//********************************************************************************************//


void FourierList_Grad_Iprod_div(Element_List *EL,
				double *f,
				double *work_1, double *work_2,
				double *work_3)
{
  // compute 
  //    \int_{volume_2.5D} [\grad(f) \cdot \grad(basis_p)_k]
  //   = \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)]
  //     + Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]
  //
  // On input:
  //   EL : Element_list structure, workspace
  //    f : function data, in (Q,F)-space, dimension: [EL->nz*EL->htot]
  //   work_1, work_2 : workspace, dimension: [EL->nz*El->htot]
  //   work_3 : workspace, dimension: at least [EL->nz*EL->hjtot]
  //
  // On exit:
  //   EL->hjtot contains integral result
  //   f : unchanged from before the call
  //   work_1, work_2:  contains df/dx and df/dy in (Q,F)-space
  //   work_3 contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]
  //

  double beta_k;

  int k;
  int nz = EL->nz;
  int k_glob;
  int htot = EL->htot,
    hjtot = EL->hjtot;
  int htot_nz = htot*nz,
    hjtot_nz = hjtot*nz;

  //----------------------------------------------//
  // first Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]

  dcopy(htot_nz, f, 1, EL->base_h, 1);

  EL->Iprod(EL);
  // Now EL->base_hj contains \int_{volume_2D} f_k * basis_p, in (J,F)-space
  //     EL->base_h unchanged, still contains function f in (Q,F)-space

  for(k=0;k<nz;k++) 
    dscal(hjtot, Beta(k)*Beta(k), EL->flevels[k]->base_hj, 1);
  // Now EL->base_hj contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], in (J,F)-space
  
  // save this integral
  dcopy(hjtot_nz, EL->base_hj, 1, work_3, 1);
  // work_3 contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], in (J,F)-space
  

  //-----------------------------------------------------------//
  // then \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)]

  // Now El->base_hj still contains function f in (Q,F)-space
  EL->Grad_d(work_1, work_2, 0, 'a');
  // Now work_1 and work_2 contains df/dx and df/dy, in (Q,F)-space
  //     EL->base_h still contains function f in (Q,F)-space, free for other use
  //       ->base_hj still contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], in (J,F)-space,
  //                       free for other use

  FourierList_Iprod_div(EL, work_1, work_2);
  // Now EL->base_hj contains 
  //             \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)], in (J,F)-space
  //       ->base_h contains junk
  //     work_1, work_2 unchanged, still contain df/dx and df/dy in (Q,F)-space
 
  //------------------------------//
  // add two integrals

  daxpy(hjtot_nz, 1.0, work_3, 1, EL->base_hj, 1);
  // Now EL->base_hj contains 
  //          \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)]
  //             + Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], 
  //        in (J,F)-space
  //     EL->base_h contains junk
  //     work_1, work_2:  contains df/dx and df/dy in (Q,F)-space
  //     work_3 contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]
  //     f unchanged from before the call


}



void FourierList_Grad_Iprod_div(Element_List *EL,
				double *work_1, double *work_2,
				double *work_3)
{
  // compute 
  //    \int_{volume_2.5D} [\grad(f) \cdot \grad(basis_p)_k]
  //   = \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)]
  //     + Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]
  //
  // On input:
  //   EL->base_h : contains function f in (Q,F)-space
  //     ->base_hj : contains nothing
  //   work_1, work_2 : workspace, dimension: [EL->nz*El->htot]
  //   work_3 : workspace, dimension: at least [EL->nz*EL->hjtot]
  //
  // On exit:
  //   EL->base_hj contains integral result
  //     ->base_h contains junk
  //   work_1, work_2:  contains df/dx and df/dy in (Q,F)-space
  //   work_3 contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]
  //

  double beta_k;

  int k;
  int nz = EL->nz;
  int k_glob;
  int htot = EL->htot,
    hjtot = EL->hjtot;
  int htot_nz = htot*nz,
    hjtot_nz = hjtot*nz;

  //----------------------------------------------//
  // first Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]

  //dcopy(htot_nz, f, 1, EL->base_h, 1);

  // Now El->base_h contains function data in (Q,F)-space
  EL->Iprod(EL);
  // Now EL->base_hj contains \int_{volume_2D} f_k * basis_p, in (J,F)-space
  //     EL->base_h unchanged, still contains function f in (Q,F)-space

  for(k=0;k<nz;k++) 
    dscal(hjtot, Beta(k)*Beta(k), EL->flevels[k]->base_hj, 1);
  // Now EL->base_hj contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], in (J,F)-space
  
  // save this integral
  dcopy(hjtot_nz, EL->base_hj, 1, work_3, 1);
  // work_3 contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], in (J,F)-space
  

  //-----------------------------------------------------------//
  // then \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)]

  // Now El->base_hj still contains function f in (Q,F)-space
  EL->Grad_d(work_1, work_2, 0, 'a');
  // Now work_1 and work_2 contains df/dx and df/dy, in (Q,F)-space
  //     EL->base_h still contains function f in (Q,F)-space, free for other use
  //       ->base_hj still contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], in (J,F)-space,
  //                       free for other use

  FourierList_Iprod_div(EL, work_1, work_2);
  // Now EL->base_hj contains 
  //             \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)], in (J,F)-space
  //       ->base_h contains junk
  //     work_1, work_2 unchanged, still contain df/dx and df/dy in (Q,F)-space
 
  //------------------------------//
  // add two integrals

  daxpy(hjtot_nz, 1.0, work_3, 1, EL->base_hj, 1);
  // Now EL->base_hj contains 
  //          \int_{volume_2D} [ \grad_{2D}(f)_k \cdot \grad_{2D}(basis_p)]
  //             + Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p], 
  //        in (J,F)-space
  //     EL->base_h contains junk
  //     work_1, work_2:  contains df/dx and df/dy in (Q,F)-space
  //     work_3 contains Beta(k)*Beta(k)* \int_{volume_2D} [ (f)_k * basis_p]
  //     f unchanged from before the call


}



//**********************************************************************//


