/*
 * mpiprism -- MPI-F version of prism3d
 *
 * C.H. Crawford -- cait@cfm.brown.EDU
 * (from Dave Newman's SGI original pvm3 version)
 *
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "nektarF.h"

#ifdef PARALLEL

#define ROOTONLY  if (mynode() == 0)
#define MSGTAG    599

static int numproc;
static int my_node;
void do_main(int argc, char *argv[]);

main (int argc, char *argv[])
{
  int info, nprocs,                      /* Number of processors */
      mytid;                             /* My task id */

  info = MPI_Init (&argc, &argv);                 /* Initialize */
  if (info != MPI_SUCCESS) {
    fprintf (stderr, "MPI initialization error\n");
    exit(1);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);         /* Number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &mytid);          /* my process id */

  numproc = nprocs;
  my_node = mytid;

  MPI_Barrier  (MPI_COMM_WORLD);                  /* sync before work */

#ifdef DEBUG_LAM
  kpause();                                       /* LAM pause utility */
#endif
#ifdef DEBUG_WAIT
  sleep(60);
#endif  

  do_main (argc, argv);                           /* wrapped program  */
  
  MPI_Finalize ();                                /* exit MPI         */

  return 0;
}


#else

void do_main(int argc, char *argv[]);
main (int argc, char *argv[]){
  do_main (argc, argv);                           /* wrapped program  */
  return 0;
}

#endif
