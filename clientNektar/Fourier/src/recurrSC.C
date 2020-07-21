/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: /homedir/cvs/Nektar/Fourier/src/recurrSC.C,v $  
 * $Revision: 1.2 $
 * $Date: 2006/05/13 09:32:45 $    
 * $Author: ssherw $  
 * $State: Exp $   
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <veclib.h>
#include "nektarF.h"

typedef struct pvert {
  int gid;
  struct pvert *next;
} Pvert;

typedef struct pedge {
  int    nelmt_vert;
  int    *vgid;
  int    end[2];
  int    nelmt_edge;
  int    *gid;
  int    conpatch;
  struct pedge *next;
} Pedge;

typedef struct patchinfo{
  int id;
  int nvert;
  int nedge;
  Pvert *vert;
  Pedge *edge;
} Pinfo;

typedef struct plevel{
  int npatch;
  int nvs;     /* global elemental value of maximum vertices gid */
  Pinfo *pbdy; /* boundary information about patch               */
  Pinfo *pint; /* interior information about patch               */
  int *old2new; /* old to new patch mapping needed for matrix setup */
} Plevel;

static Plevel  *setup_patch_info     (Element *E, int nel, int nvs, int nes);
static void     free_patch_info      (Plevel *P);
static int      find_next_patch      (Plevel *P, Plevel **P1);
static void     print_plevel         (char *out, Plevel *P, Element *E,
				      Bsystem *B);
static Rsolver *setup_recur_numbering(Plevel **plev, int nrecur,
				      Element *E, Bsystem *B);

void Recursive_SC_decom(Element *E, Bsystem *B){
  register int cnt,trip;
  int nrecur;
  Plevel **plev,**Porig,**Pnew;
  Rsolver  *rslv;
  int nvs = B->nv_solve;
  int nes = B->ne_solve;

  nrecur = option("recursive"); 

  plev   = (Plevel **)malloc(nrecur*sizeof(Plevel *));
  Porig  = plev;
  Pnew   = plev+1;

  *Porig = setup_patch_info(E,B->nel,nvs,nes);
  cnt = 0;
  trip = nrecur;
      
  /* go through patches and set up patch information           */
  while(trip-- && find_next_patch(plev[cnt],Pnew + cnt))
    cnt++;
  
  if(cnt) B->rslv = rslv = setup_recur_numbering(Pnew,cnt,E,B);
 
  ROOT {
    if(!rslv)
      fprintf(stderr,"No recursions were possible\n");
    else if(nrecur > rslv->nrecur)
      fprintf(stderr,"Only %d recursions were possible\n",rslv->nrecur);
  }
    
  for(trip = 0; trip <= cnt; ++trip)
    free_patch_info(plev[trip]); 

  free(plev);
}

static Fctcon *SCFacetConnect(int nsols, Plevel *P, int *vmap, int *emap);

static Rsolver *setup_recur_numbering(Plevel **plev, int nrecur,
				      Element *E, Bsystem *B){
  Element      *F;
  register int i,j,k;
  Rsolver *rslv;
  Recur   *r;
  int     *ngid, *edge, *edglen;
  int     id,cnt,cnt1,nedge,one=1;
  int     nvs = B->nv_solve;
  int     nes = B->ne_solve;
  Pinfo   *pi,*pb;
  Pvert   *pv;
  Pedge   *pe;

  /* check to see if the last recursion has more than one patch, if
     not then reduce nrecur  */

  if(plev[nrecur-1]->npatch == 1) nrecur--;
  
  if(nrecur< 1) return (Rsolver *) NULL;
  
  rslv = (Rsolver *)calloc(1,sizeof(Rsolver));
  rslv->nrecur    = nrecur;
  r = rslv->rdata = (Recur *)calloc(nrecur,sizeof(Recur));

  for(i = 0; i < nrecur; ++i)
    r[i].pmap = plev[i]->old2new;

  ngid   = ivector(0,nvs);
  edge   = ivector(0,nes);
  edglen = ivector(0,nes);
  
  /* set up original list of edge lengths */
  for(F = E; F; F=F->next)
    for(j = 0; j < F->Nedges; ++j)
      if((id = F->edge[j].gid) < nes) 
	edglen[id] = F->edge[j].l;  
    
  /* setup inner most boundary system numbering */
  /* find unsolved edge and vertex id's */
  
  ifill(nes,-1,edge,1);
  ifill(nvs,-1,ngid,1);
  
  /* set up basic boundary numbering scheme for iterative solver */
  /* put global vertices first */
  cnt = 0;
  for(i = 0; i < plev[nrecur-1]->npatch; ++i){
    pb = plev[nrecur-1]->pbdy + i;
    
    for(pv = pb->vert; pv; pv = pv->next)
      if(ngid[pv->gid] == -1) ngid[pv->gid] = cnt++;
  }
  rslv->Ainfo.nv_solve = cnt;
  
  for(i = 0; i < plev[nrecur-1]->npatch; ++i){
    pb = plev[nrecur-1]->pbdy + i;
    
    /* local vertices  */
    for(pe = pb->edge; pe; pe = pe->next){
      for(k = 0; k < pe->nelmt_vert; ++k)
	if(ngid[pe->vgid[k]] == -1) ngid[pe->vgid[k]] = cnt++;
    }
    
    /* local edges */
    for(pe = pb->edge; pe; pe = pe->next)
      for(k = 0; k < pe->nelmt_edge; ++k)
	if(edge[pe->gid[k]] == -1){
	  edge[pe->gid[k]] = cnt; 
	  cnt += edglen[pe->gid[k]];
	}
  }

  if(B->smeth == direct){  /* reverse cuthill mckee sorting if direct method */
    Fctcon  *ptcon;
    int     nvert,nedge,nsols;
    int     *bwdm,*newmap;

    bwdm = ivector(0,nes+nvs-1);

    /* set up consequative list of unknowns (i.e. vertices and edges)
       as well as the backward mapping                                */

    nvert = 0;
    for(i = 0; i < nvs; ++i)
      if(ngid[i]+1) {
	ngid[i] = nvert;
	bwdm[nvert++] = i;
      }

    nedge = 0;
    for(i = 0; i < nes; ++i)
      if(edge[i]+1){
	edge[i] = nvert + nedge;
	bwdm[nvert + nedge++] = i;
      }
    nsols = nvert + nedge;
    
    newmap = ivector(0,nsols-1);
    ptcon  = SCFacetConnect(nsols,plev[nrecur-1],ngid,edge);
    
    MinOrdering(nsols,ptcon,newmap);

    /* make up list of sorted global gids */
    for(i = 0,cnt = 0; i < nsols; ++i)
      if(newmap[i] < nvert)
	ngid[bwdm[newmap[i]]] = cnt++;
      else{
	edge[bwdm[newmap[i]]] = cnt;
	cnt += edglen[bwdm[newmap[i]]];
      }	  
    
    free(bwdm); free(newmap);
    free_Fctcon(ptcon,nsols);
  }
  
  /* go through patches and set up patch information           *
   * do this in reverse order to make the numbering consistent *
   * with the standard ordering                                */
  
  for(i = nrecur-1; i >= 0; --i){
    r[i].id = i;
    r[i].npatch     = plev[i]->npatch;
    r[i].patchlen_a = ivector(0,r[i].npatch-1);
    r[i].patchlen_c = ivector(0,r[i].npatch-1);
    
    /* set up new ordering for interior patches -- putting the
       vertices in the middle */
    
    for(j = 0; j < r[i].npatch; ++j){
      pi = plev[i]->pint + j;
      cnt1 = cnt;
      nedge = pi->nedge/2;
      for(pe = pi->edge; pe&&nedge; pe = pe->next, --nedge){

	for(k = 0; k < pe->nelmt_edge/2; ++k){
	  edge[pe->gid[k]] = cnt;
	  cnt += edglen[pe->gid[k]];
	}    

	for(k = 0; k < pe->nelmt_vert; ++k)
	  ngid[pe->vgid[k]] = cnt++;
	
	for(k = pe->nelmt_edge/2; k < pe->nelmt_edge; ++k){
	  edge[pe->gid[k]] = cnt;
	  cnt += edglen[pe->gid[k]];
	}    
      }

      for(pv = pi->vert; pv; pv = pv->next)
	ngid[pv->gid] = cnt++;

      for(; pe; pe=pe->next){
	for(k = 0; k < pe->nelmt_edge/2; ++k){
	  edge[pe->gid[k]] = cnt;
	  cnt += edglen[pe->gid[k]];
	}    

	for(k = 0; k < pe->nelmt_vert; ++k)
	  ngid[pe->vgid[k]] = cnt++;
	
	for(k = pe->nelmt_edge/2; k < pe->nelmt_edge; ++k){
	  edge[pe->gid[k]] = cnt;
	  cnt += edglen[pe->gid[k]];
	}    
      }
      
      r[i].patchlen_c[j] = cnt - cnt1;
      
      /* calculate the boundary length */
      
      pi = plev[i]->pbdy + j;
      cnt1 = pi->nvert;
      for(pe = pi->edge; pe; pe = pe->next){
	cnt1 += pe->nelmt_vert;
	for(k = 0; k < pe->nelmt_edge; ++k)
	  cnt1 += edglen[pe->gid[k]];
      }
      r[i].patchlen_a[j] = cnt1;
      rslv->max_asize = max(rslv->max_asize,cnt1);
    }
  }

  /* copy back gid's to vertices */
  for(F = E; F; F = F->next)
    for(j = 0; j < F->Nverts; ++j)
      if((id = F->vert[j].gid) < nvs) 
	F->vert[j].gid = ngid[id];

  /* copy back edge id's */
  icopy(nes,edge,1,B->edge,1);
  
  /* set up start location interior systems */
  cnt = B->nsolve;
  for(i = 0; i < nrecur; ++i){
    r[i].cstart = cnt - isum(r[i].npatch,r[i].patchlen_c,1);
    cnt = r[i].cstart;
  }
  
  /* calculate the local to global boundary mapping */
  for(i = 0; i < nrecur; ++i){
    r[i].map = (int **)malloc(r[i].npatch*sizeof(int *));
    for(j = 0; j < r[i].npatch; ++j)
      if(r[i].patchlen_a[j]){
	r[i].map[j] = ivector(0,r[i].patchlen_a[j]-1);
	
	pi = plev[i]->pbdy + j;

	cnt = 0;
	/* put global vertices first at present */
	for(pv = pi->vert;pv;pv = pv->next)
	  r[i].map[j][cnt++] = ngid[pv->gid];

	/* then put local vertices */
	for(pe = pi->edge; pe; pe = pe->next){
	  for(k = 0; k < pe->nelmt_vert; ++k)
	    r[i].map[j][cnt++] = ngid[pe->vgid[k]];
	}	  
	
	/* finally store global vertices */
	for(pe = pi->edge; pe; pe = pe->next){
	  for(k = 0; k < pe->nelmt_edge; ++k){
	    iramp(edglen[pe->gid[k]],edge+pe->gid[k],&one,r[i].map[j]+cnt,1);
	    cnt += edglen[pe->gid[k]];
	  }
	}
      }
  }
    
  if((B->smeth == iterative)&&(B->Precon == 1)){
    Blockp *Blk;

    rslv->precon = (Precond *)calloc(1,sizeof(Precond));
    
    Blk = &rslv->precon->blk; 
    
    Blk->ngv = ivector(0,r[nrecur-1].npatch-1);
    Blk->nle = ivector(0,r[nrecur-1].npatch-1);

    Blk->edglen = (int **)malloc(r[nrecur-1].npatch*sizeof(int*));
    Blk->lgid   = (int **)malloc(r[nrecur-1].npatch*sizeof(int*));

    /* set up each block to be used along and edge with a unique global id */
    ifill(nes,-1,edge,1);
    ifill(nvs,-1,ngid,1);
    cnt = 0;

    for(i = 0; i < r[nrecur-1].npatch; ++i){
      pi = plev[nrecur-1]->pbdy + i;
      
      for(pe = pi->edge; pe; pe = pe->next){
	if(pe->nelmt_vert)
	  if(ngid[pe->vgid[0]] == -1){
	    for(k = 0; k < pe->nelmt_vert; ++k)
	      ngid[pe->vgid[k]] = cnt;
	    ++cnt;
	  }
      }
      
      for(pe = pi->edge; pe; pe = pe->next)
	for(k = 0; k < pe->nelmt_edge; ++k) 
	  if(edge[pe->gid[k]] == -1) 
	    edge[pe->gid[k]] = cnt++;
      
    }
    Blk->nlgid = cnt;

    for(i = 0; i < r[nrecur-1].npatch; ++i){
      pi = plev[nrecur-1]->pbdy + i;
    
      cnt = 0;
      /* count global vertices  */
      for(pv = pi->vert;pv;pv = pv->next) cnt++;
      Blk->ngv[i] = cnt;
      
      /* count local vertex blocks */
      for(cnt = 0, pe = pi->edge; pe; pe = pe->next){
	/* if vertices exist count them as one block */
	if(pe->nelmt_vert) cnt++;
      
	/* count edge blocks as well */
	cnt += pe->nelmt_edge;
      }
      Blk->nle[i] = cnt;

      Blk->edglen[i] = ivector(0,Blk->nle[i]-1);
      Blk->lgid[i]   = ivector(0,Blk->nle[i]-1);

      for(cnt=0,pe = pi->edge; pe; pe = pe->next)
	if(pe->nelmt_vert){
	  Blk->lgid[i][cnt] = ngid[pe->vgid[0]];
	  Blk->edglen[i][cnt++] = pe->nelmt_vert;
	}

      for(pe = pi->edge; pe; pe = pe->next)
	for(k = 0; k < pe->nelmt_edge; ++k){
	  Blk->edglen[i][cnt] = edglen[pe->gid[k]];
	  Blk->lgid[i][cnt++] = edge[pe->gid[k]];
	}
    }
  }

  free(ngid); free(edglen); free(edge);
  return rslv;
}

static Fctcon *SCFacetConnect(int nsols, Plevel *P, int *vmap, int *emap)
{
  register int i,k;
  const    int npatch = P->npatch;
  int        *pts, n;
  Fctcon     *connect;
  Pinfo      *pb = P->pbdy;
  Pvert      *pv;
  Pedge      *pe;

  connect = (Fctcon *)calloc(nsols,sizeof(Fctcon));
     
  pts = ivector(0,nsols-1);

  for(k = 0; k < npatch; ++k){
    /* find all facets that are unknowns */
    n = 0;
    for(pv = pb[k].vert; pv; pv = pv->next)
      pts[n++] = vmap[pv->gid];

    for(pe = pb[k].edge;pe;pe = pe->next){
      for(i = 0; i < pe->nelmt_vert; ++i)
	pts[n++] = vmap[pe->vgid[i]];
    
      for(i = 0; i < pe->nelmt_edge; ++i)
	pts[n++] = emap[pe->gid[i]];
    }

    addfct(connect,pts,n);
  }
    
  free(pts);

  return connect;
}


typedef struct plist{
  int id;
  struct plist *next;
} Plist;

static int  find_next_patch(Plevel *P, Plevel **P1){
  register int i,j,k,n;
  int   **new2old,*old2new,newpat,*index;
  int    *newpatlen,*p,*mult,*vorder;
  int    gid, trip, id, plen, count;
  Pvert  *pv,*pv1, *pv2;
  Pedge  *pe,*pe1, *pe2;
  Pinfo  *pi,*pb;
  Plist  *pl,*pl1,*pl2;

  if((P->npatch <= 2 )||!P->nvs) return 0; /* can't condense any lower */

  mult    = ivector(0,P->nvs-1);
  old2new = ivector(0,P->npatch-1);
  ifill(P->npatch,-1,old2new,1);
  
  pl = (Plist *)calloc(P->nvs,sizeof(Plist));

  /* find the patch id's surrounding all the vertices */
  for(i = 0; i < P->npatch; ++i)
    for(pv = P->pbdy[i].vert; pv; pv = pv->next){
      pl1     = (Plist *)malloc(sizeof(Plist));
      pl1->id = i;
      pl1->next = pl[pv->gid].next;
      pl[pv->gid].next = pl1;
    }
  
  /* find the multiplicity of each new vertex */
  for(i = 0; i < P->nvs; ++i)
    for(pl1 = pl[i].next, mult[i] = 0; pl1; pl1 = pl1->next)
      mult[i]++;
  
  /* count the number of global vertices in patch */
  for(i =0, count = 0; i < P->nvs; ++i) if(mult[i]) count++;
  
  if(!count){ /* deal with case where all patches don't have a global vert */
    free(mult);
    free(old2new);

    /* free vertex  list */
    for(i = 0; i < P->nvs; ++i){
      pl1 = pl[i].next;
      while(pl1){
	pl2 = pl1->next;
	free(pl1);
	pl1 = pl2;
      }
    }   
    free(pl);
    return 0;
  }
    

  vorder = ivector(0,count-1);

  for(i = 0, n = 0; i < P->nvs; ++i)
    if(mult[i]) vorder[n++] = i;  

  if(n != count) error_msg(did not sort vertices correctly in find_next_patch);
  
  /* Finally sort out old to new patch mapping around vertices */
  newpat = 0;
  for(i = 0; i < count; ++i){
    /* check to see if any patch around vertex is part of a new patch */
    trip = 1;
    pl1 = pl[vorder[i]].next;
    while(pl1){
      if(old2new[pl1->id]+1) trip = 0;
      pl1 = pl1->next;
    }
    if(trip){ /* if not then set up a new patch */
      pl1 = pl[vorder[i]].next;
      while(pl1){
	old2new[pl1->id] = newpat;
	pl1 = pl1->next;
      }
      newpat++;
    }
  }  

  if(!newpat){/* this will condense all the remaining patches into 1 */
    newpat = 1; 
    old2new[0] = 0;
  }

  /* measure length of present patches */
  newpatlen = ivector(0,newpat-1);
  izero(newpat,newpatlen,1);
  for(i = 0; i < P->npatch; ++i) 
    if(old2new[i]+1) newpatlen[old2new[i]]++;
  
  /* check for any loose patch and attach to neighbouring patch with
     lowest number of components. Add all remaining patches to an
     existing patch (if there is more than one choice use the one with
     the lowest number of components) This search may have to be
     repeated more than once for a general mesh */

  count = 100;
  while((isum(newpat,newpatlen,1) != P->npatch)&&--count)
    for(i = 0; i < P->npatch; ++i)
      if(old2new[i] == -1){
	for(pe = P->pbdy[i].edge; pe; pe = pe->next){
	  plen = 100000000;
	  if(pe->conpatch+1)
	    if(old2new[pe->conpatch]+1)
	      if(newpatlen[old2new[pe->conpatch]] < plen){
		old2new[i] = old2new[pe->conpatch];
		plen = newpatlen[old2new[pe->conpatch]];
	      }
	}
	if(old2new[i]+1) newpatlen[old2new[i]]++;
      }
  
  if(!count) error_msg(Failed to generate new patch in find_next_patch);
  
  /* set up new to old patch mapping */
  index = ivector(0,newpat-1);
  ifill(newpat,0,index,1);

  new2old = imatrix(0,newpat-1,0,newpatlen[iimax(newpat,newpatlen,1)]-1);

  for(i = 0; i < P->npatch; ++i){
    new2old[old2new[i]][index[old2new[i]]] = i;
    index[old2new[i]]++;
  }
  
  /* setup newpatch */
  P1[0] = (Plevel *)calloc(1,sizeof(Plevel));
  P1[0]->npatch = newpat;
  P1[0]->pbdy   = (Pinfo *)calloc(newpat,sizeof(Pinfo));
  P1[0]->pint   = (Pinfo *)calloc(newpat,sizeof(Pinfo));
  P1[0]->nvs    =  P->nvs;
  P1[0]->old2new = old2new;

  /* setup new patches with all the vertices and edge of old patches */
  for(i = 0; i < newpat; ++i){
    pb = P1[0]->pbdy+i;
    pi = P1[0]->pint+i;
    pb->id = i;
    pi->id = i;


    for(j = 0; j < newpatlen[i]; ++j){
      for(pv1 = P->pbdy[new2old[i][j]].vert; pv1; pv1 = pv1->next){
	pv = (Pvert *)malloc(sizeof(Pvert));
	memcpy(pv,pv1,sizeof(Pvert));
	pv->next = pb->vert;
	pb->vert = pv;
      }

      for(pe1 = P->pbdy[new2old[i][j]].edge; pe1; pe1 = pe1->next){
	/* this test is to stop duplicate copies of the interior edges */
	trip = 1;
	if(pe1->conpatch + 1)
	  if((pe1->conpatch < new2old[i][j])&&(old2new[pe1->conpatch] == i))
	    trip = 0;

	if(trip){
	  pe = (Pedge *)calloc(1,sizeof(Pedge));
	  memcpy(pe,pe1,sizeof(Pedge));
	  if(pe->nelmt_vert){
	    pe->vgid = ivector(0,pe->nelmt_vert-1);
	    icopy(pe->nelmt_vert,pe1->vgid,1,pe->vgid,1);
	  }
	  pe->gid  = ivector(0,pe->nelmt_edge-1);
	  icopy(pe->nelmt_edge,pe1->gid,1,pe->gid,1);
	  if(pe1->conpatch != -1)
	    pe->conpatch = old2new[pe1->conpatch];
	  else
	    pe->conpatch = -1;
	  pe->next = pb->edge;
	  pb->edge = pe;
	}
      }
    }
    
    /* eliminate multiply defined vertices */
    for(pv = pb->vert; pv; pv = pv->next){
      gid = pv->gid;
      pv1 = pv;
      while(pv1->next){ /* if  patch has duplicated vertex remove patch */
	if((pv1->next->gid == gid)){
	  pv2 = pv1->next;
	  pv1->next = pv2->next;
	  free(pv2);
	}
	else
	  pv1 = pv1->next;
      }
    }

    
    /* remove any interior vertices from pbdy and put in pint */
    pv = pb->vert;
    pb->vert = (Pvert *)NULL;
    while(pv){
      for(trip = 1,pl1 = pl[pv->gid].next; pl1; pl1 = pl1->next)
	if(old2new[pl1->id] != i) trip &= 0;

      if(trip){ /* put on pint list delete element */
	pv1 = pv->next;
	pv->next = pi->vert;
	pi->vert = pv;
	pv = pv1;
      }
      else{ /* add back onto original pbdy list */
	pv1 = pv->next;
	pv->next = pb->vert;
	pb->vert = pv;
	pv = pv1;
      }
    }
    
    /* put internally connected edges where conpatch == i and flux 
       boudaries where conpatch == -1 in pint */
    pe = pb->edge;
    pb->edge = (Pedge *)NULL;
    while(pe){
      if(pe->conpatch == i || pe->conpatch == -1){ /* put in pint */
	pe1 = pe->next;
	pe->next = pi->edge;
	pi->edge = pe;
	pe = pe1;
      }
      else{                  /* put edge back onto link list off pbdy */
	pe1 = pe->next;
	pe->next = pb->edge;
	pb->edge = pe;
	pe = pe1;
      }
    }
    
    /* finally condense edges and vertices which only 
       lie between the same two patches */
    pv = pb->vert;
    pb->vert = (Pvert *)NULL;

    while(pv){
      
      /* search through elements around vertex and find out how many 
	 patches are attached */
      id = old2new[pl[pv->gid].next->id];
      trip = 0;
      for(gid=-1,pl1 = pl[pv->gid].next; pl1; pl1 = pl1->next){
	if(old2new[pl1->id] - id){
	  if(!(gid+1)) gid = old2new[pl1->id];
	  if(old2new[pl1->id] - gid){
	    trip = 1;
	    break;
	  }
	}
      }
      
      /* if there are only two patches surrounding this vertex then check 
	 to see if it's between two edges and if so condense  */

      if(!trip){ 

	pe = pb->edge;
	while(pe){ /* find first edge adjacent to vertex */
	  if(pe->end[0] == pv->gid || pe->end[1] == pv->gid) break;
	  else pe = pe->next;
	}

	/* add this vertex to the edge and then delete after next check */
	p = ivector(0,pe->nelmt_vert);
	if(pe->nelmt_vert){
	  icopy(pe->nelmt_vert,pe->vgid,1,p,1);
	  free(pe->vgid);
	}
	pe->vgid = p;
	pe->nelmt_vert += 1;
	p[pe->nelmt_vert-1] = pv->gid;
	
	/* check for any other adjacent edges and condense if they exist */
	pe1 = pe->next;
	pe->next = (Pedge *)NULL;
	while(pe1){
	  /* if adjacent condense */
	  if(pe1->end[0] == pv->gid || pe1->end[1] == pv->gid){
	    p = ivector(0,pe->nelmt_vert + pe1->nelmt_vert-1);
	    if(pe->nelmt_vert){
	      icopy(pe ->nelmt_vert,pe->vgid,1,p,1);
	      free(pe->vgid);
	    }
	    if(pe1->nelmt_vert){
	      icopy(pe1->nelmt_vert,pe1->vgid,1,p + pe->nelmt_vert,1);
	      free(pe1->vgid);
	    }
	    pe->vgid = p;
	    pe->nelmt_vert +=  pe1->nelmt_vert;
	    
	    /* set up new end id's */
	    for(j = 0; j < 2; ++j)
	      if(pe->end[j] == pv->gid)
		pe->end[j] = (pe1->end[0] != pv->gid)? 
		  pe1->end[0]:pe1->end[1];


	    /* sort vgid's so that the lowest id is first -- This is only
	       really necessary for the iterative solve where we want
	       the vertices to be stored in a consistent block with side
	       of a patch */

	    for(j = 0; j < pe->nelmt_vert; ++j){
	      for(id = k = j; k < pe->nelmt_vert; ++k)
		id = (pe->vgid[k] < pe->vgid[id])? k: id;
	      trip = pe->vgid[id];
	      for(k = id; k > j; --k)
		pe->vgid[k] = pe->vgid[k-1];
	      pe->vgid[j] = trip;
	    }
	    
	    /* sort out local edges */
	    p = ivector(0,pe->nelmt_edge + pe1->nelmt_edge);
	    icopy(pe ->nelmt_edge,pe->gid,1,p,1);
	    icopy(pe1->nelmt_edge,pe1->gid,1,p+pe->nelmt_edge,1);
	    
	    free(pe->gid); free(pe1->gid);
	    pe->gid = p;
	    pe->nelmt_edge +=  pe1->nelmt_edge;
	      
	    pe2 = pe1->next;
	    free(pe1);
	    pe1 = pe2;


	  }
	  else{ /* if not adjacent add back to list */
	    pe2 = pe1->next;
	    pe1->next = pe->next;
	    pe->next  = pe1;
	    pe1 = pe2;	    
	  }
	}
	/* free up this vertex */
	pv1 = pv;
	pv = pv->next;
	free(pv1);
      }
      else{
	pv1 = pv->next; /* add vertex onto list */
	pv->next = pb->vert;
	pb->vert = pv;
	pv = pv1;
      }
    }
    /* count up number of vertices */
    for(pv = pb->vert; pv; pv = pv->next) pb->nvert ++;
    for(pv = pi->vert; pv; pv = pv->next) pi->nvert ++;
    /* count up number of edges */
    for(pe = pb->edge; pe; pe = pe->next) pb->nedge ++;    
    for(pe = pi->edge; pe; pe = pe->next) pi->nedge ++;    
  }

  /* free vertex  list */
  for(i = 0; i < P->nvs; ++i){
    pl1 = pl[i].next;
    while(pl1){
      pl2 = pl1->next;
      free(pl1);
      pl1 = pl2;
    }
  }

  free(newpatlen); free(index);
  free_imatrix(new2old,0,0);
  return 1;
}

static Plevel *setup_patch_info(Element *E, int nel, int nvs, int nes){
  Element      *F;
  register int i,k;
  Plevel *P;
  Pvert *pv;
  Pedge *pe;
  Pinfo *p;
  
  P = (Plevel *)calloc(1,sizeof(Plevel));
  P->npatch = nel;
  P->nvs    = nvs;
  P->pbdy   = p = (Pinfo *)calloc(nel,sizeof(Pinfo));

  for(F = E; F; F = F->next){
    k = F->id;
    p[k].id = k;

    for(i = 0; i < F->Nverts; ++i)
      if(F->vert[i].gid < nvs){
	p[k].nvert += 1;
	pv = (Pvert *)malloc(sizeof(Pvert));
	pv->gid = F->vert[i].gid;
	pv->next = p[k].vert;
	p[k].vert = pv;
      }
    
    for(i = 0; i < F->Nedges; ++i)
      if(F->edge[i].gid < nes){
	p[k].nedge += 1;
	pe = (Pedge *)calloc(1,sizeof(Pedge));
	pe->nelmt_vert = 0;
	pe->end[0] = F->vert[i].gid;
	pe->end[1] = F->vert[(i+1)%F->Nverts].gid;
	pe->nelmt_edge = 1;
	pe->gid = ivector(0,0);
	pe->gid[0] = F->edge[i].gid;
	if(F->edge[i].base)
	  pe->conpatch = (F->edge[i].link)?  
	    F->edge[i].link->eid: F->edge[i].base->eid;
	else
	  pe->conpatch = -1;
	pe->next = p[k].edge;
	p[k].edge = pe;
      }    
  }

  return P;
}

static void free_patch_info(Plevel *P){
  register int i;
  Pvert *pv1, *pv2;
  Pedge *pe1, *pe2;
  int npatch = P->npatch;
  Pinfo *p;

  p = P->pbdy;
  for(i = 0; i < npatch; ++i){
    if(p[i].nvert){
      pv1 = p[i].vert;
      while(pv1){
	pv2 = pv1->next;
	free(pv1);
	pv1 = pv2;	
      }    
    }

    if(p[i].nedge){
      pe1 = p[i].edge;
      while(pe1){
	pe2 = pe1->next;
	if(pe1->nelmt_vert) free(pe1->vgid);
	free(pe1->gid);
	free(pe1);
	pe1 = pe2;	
      }    
    }
  }
  free(p);

  if(p = P->pint){
    for(i = 0; i < npatch; ++i){
      if(p[i].nvert){
	pv1 = p[i].vert;
      while(pv1){
	pv2 = pv1->next;
	free(pv1);
	pv1 = pv2;	
      }    
      }
      
      if(p[i].nedge){
	pe1 = p[i].edge;
	while(pe1){
	  pe2 = pe1->next;
	  if(pe1->nelmt_vert) free(pe1->vgid);
	  free(pe1->gid);
	  free(pe1);
	  pe1 = pe2;	
	}    
      }
    }
    free(p);
  }  
  free(P);
}
