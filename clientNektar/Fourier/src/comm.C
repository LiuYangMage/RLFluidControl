#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <math.h>
#include <veclib.h>
#include "nektarF.h"

#ifdef PARALLEL

void unreduce (double *x, int n);
void reduce   (double *x, int n, double *work);
static int numproc;
static int my_node;

void gsync ()
{
  int info;

  info = MPI_Barrier(MPI_COMM_WORLD);

  return;
}


int numnodes ()
{
  int np;

  MPI_Comm_size(MPI_COMM_WORLD, &np);         /* Number of processors */

  return np;
}


int mynode ()
{
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);          /* my process id */
  
  return myid;
}

 
void csend (int type, void *buf, int len, int node, int pid)
{

  MPI_Send (buf, len, MPI_BYTE, node, type, MPI_COMM_WORLD);
  
  return;
}

void crecv (int typesel, void *buf, int len)
{
  MPI_Status status;

  MPI_Recv (buf, len, MPI_BYTE, MPI_ANY_SOURCE, typesel, MPI_COMM_WORLD, &status);
  
  return;
}


void msgwait (MPI_Request *request)
{
  MPI_Status status;

  MPI_Wait (request, &status);

  return;
}

double dclock(void)
{
  double time;

  time = MPI_Wtime();

  return time;

}
void gimax (int *x, int n, int *work) { 
  register int i;

  MPI_Allreduce (x, work, n, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* *x = *work; */
  icopy(n,work,1,x,1);

  return;
}

void gdmax (double *x, int n, double *work)
{
  register int i;

  MPI_Allreduce (x, work, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  /* *x = *work; */
  dcopy(n,work,1,x,1);

  return;
}


void gdsum (double *x, int n, double *work)
{
  register int i;

  MPI_Allreduce (x, work, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* *x = *work; */
  dcopy(n,work,1,x,1);

  return;
}

void gisum (int *x, int n, int *work)
{
  register int i;

  MPI_Allreduce (x, work, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* *x = *work; */
  icopy(n,work,1,x,1);

  return;
}


void unreduce (double *x, int n)
{
  int nprocs = numnodes(),
      pid    = mynode(),
      k, i;

  ROOTONLY
    for (k = 1; k < nprocs; k++)
      csend (MSGTAG + k, x, n*sizeof(double), k, 0);
  else
    crecv (MSGTAG + pid, x, n*sizeof(double));
  
  return;
}

void reduce (double *x, int n, double *work)
{
  int nprocs = numnodes(),
      pid    = mynode(),
      k, i;

  ROOTONLY {
    for (i = 0; i < n; i++) work[i] = x[i];
    for (k = 1; k < nprocs; k++) {
      crecv (MSGTAG + k, x, n*sizeof(double));
      for (i = 0; i < n; i++) work[i] += x[i];
    }
    for (i = 0; i < n; i++) x[i] = work[i];    
  } else
    csend (MSGTAG + pid, x, n*sizeof(double), 0, 0);
  
  return;
}

void ifexists(double *in, double *inout, int *n, MPI_Datatype *size){
  int i;
  
  for(i = 0; i < *n; ++i)
    inout[i] = (in[i] != 0.0)? in[i] : inout[i];
  
}



#ifdef METIS /* redefine default partitioner to be metis */


extern "C"
{
#include <metisproto.h>
}
#define pmetis(nel,xadj,adjncy,vwgt,ewgt,wflag,nparts,option,num\
	       ,edgecut,partition)\
(_vlib_ireg[0]=nel,_vlib_ireg[1]=wflag,_vlib_ireg[2]=num,_vlib_ireg[3]=nparts,\
  METIS_PartGraphRecursive(_vlib_ireg,xadj,adjncy,vwgt,ewgt,_vlib_ireg+1,\
			   _vlib_ireg+2,_vlib_ireg+3,option,edgecut,partition))

void default_partitioner(Element_List *EL){
  register int i,j;
  int eDIM  = EL->fhead->dim();
  int nel = EL->nel;
  int medg,edgecut,cnt;
  int *xadj, *adjncy, *partition;
  int opt[5];
  Element *E;

  ROOTONLY
    fprintf(stderr,"Partitioner         : using pmetis \n");
      
  /* count up number of local edges in patch */
  medg =0;
  if(eDIM == 2)
    for(E = EL->fhead; E; E= E->next){
      for(j = 0; j < E->Nedges; ++j) 
	if(E->edge[j].base) ++medg;
    }
  else
    for(E = EL->fhead; E; E= E->next){
      for(j = 0; j < E->Nfaces; ++j) 
	if(E->face[j].link) ++medg;
    }
  
  xadj      = ivector(0,nel);
  adjncy    = ivector(0,medg-1);
  partition = ivector(0,nel);
  
  izero(nel+1,xadj,1);
  cnt = 0;
  if(eDIM == 2)
    for(i = 0; i < nel; ++i){
      E = EL->flist[i];
      xadj[i+1] = xadj[i];
      for(j = 0; j < E->Nedges; ++j){
	if(E[i].edge[j].base){
	  if(E[i].edge[j].link){
	    adjncy[cnt++] = E->edge[j].link->eid;
	    xadj[i+1]++;
	  }
	  else{
	    adjncy[cnt++] = E->edge[j].base->eid;
	    xadj[i+1]++;
	  }
	}
      }
    }
  else
    for(i = 0; i < nel; ++i){
      E = EL->flist[i];
      xadj[i+1] = xadj[i];
      for(j = 0; j < E->Nfaces; ++j) 
	if(E->face[j].link){
	  adjncy[cnt++] = E->face[j].link->eid;
	  xadj[i+1]++;
	}
    }
  
  opt[0] = 0;
  pmetis(nel,xadj,adjncy,0,0,0,pllinfo.nprocs,opt,0,
	 &edgecut,partition);
  
  /* extract data on current partition */
  pllinfo.nloop = 0;
    for(i = 0; i < nel; ++i)
      if(partition[i] == pllinfo.procid) pllinfo.nloop++;
  
  pllinfo.eloop = ivector(0,pllinfo.nloop-1);
  
  cnt = 0;
  for(i = 0; i < nel; ++i)
    if(partition[i] == pllinfo.procid) pllinfo.eloop[cnt++] = i;
  
  free(xadj); free(adjncy); free(partition);
}
#endif

#endif



