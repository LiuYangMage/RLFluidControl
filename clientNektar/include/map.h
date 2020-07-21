typedef enum {               /* --------- List of flags --------- */
  OFF            = 0,        /* General-purpose flags for on/off  */
  ON             = 1,        /*                                   */
  ERROR          = 2,        /* Non-recoverable error             */
  WARNING        = 3         /* Not too serious                   */
} FLAG;


/* INITIALIZATION */

void     InitMap        (Domain *omega);
void s_InitMap (Domain *omega, int KP);
// Map     *allocate_map   (int nz);
Mapping     *allocate_map   (int nz);
// void     SaveBCs        (Domain *omega);
// double  *saveBC         (Bedge *Xbc);
void     ResetICs       (Domain *omega);
void s_ResetICs (Domain *omega, int KP);
void     ResetWlast     (Element *W, Element *Wlast);

/* UPDATING */

void UpdateMap (Domain *omega);
//zwang subiteration
void UpdateMap (Domain *omega, int bnoiterate);

void s_UpdateMap (Domain *omega, int k, int KP);
void gather_baseP (Domain *omega, double *basep);
void s_gather_baseP (Domain *omega, double *basep, int k);
void UpdateVbc (Domain *omega);
void s_UpdateVbc (Domain *omega, int k, int KP);

void MapField (Domain *omega, int dir); /*Apply the mapping & inverse mapping*/
void Mapfield (Domain *omega, int dir); /*Apply the mapping & inverse mapping*/
#ifdef AXIAL
//void Mapfield_w (Domain *omega, int dir); /*Apply the mapping & inverse mapping*/
#endif
// void gradz_map (Map *map);
void gradz_map (Mapping *map);
void gradz_map_new (Mapping *map);

// void update_forc (Map *map, double amp, double freq, double phiz, double phit);
// void update_forc (Map *map, double amp, double freq, double phiz, double phit, double BETA, int forc);
// void update_free (Map *map, double dt, double mass, double wn, double wnc, double wnb, double zeta);
// void update_free (Map *map, double dt, double mass, double wn, double wnc, double wnb, double zeta, int forc);

void update_forc (Mapping *map, double amp, double freq, double phiz, double phit);
void update_forc (Mapping *map, double amp, double freq, double phiz, double phit,
    double BETA, int forc);
#ifdef FLOW_CONTROL
void update_forc (Mapping *map, double amp, double freq, double phiz, double phit,
		  double BETA, int forc,  double *hilbert_coeff, double *dy, double *ty);
#endif

void update_free (Mapping *map, double dt, double mass, double wn, double wnc,
    double wnb, double zeta);
#ifdef AXIAL 
void update_free (Mapping *map, double dt, double mass, double wn, double wnc,
    double wnb, double wn_ea,  double zeta, int forc);
#endif
void update_free (Mapping *map, Mapping *mapz, double dt, double mass, double wn, double wnc,
    double wnb, double wn_ea,  double zeta, int forc);
void update_free (Mapping *map, Mapping *mapx, Mapping *mapy, double dt, double mass, double wn, double wnc,
    double wnb, double wn_ea,  double zeta, int forc);
void update_free (Mapping *map, Mapping *mapz, double dt, double mass, double wn, double wnc,
    double wnb, double wn_ea,  double zeta, int forc, double *Tension, double *Stiffness);
void update_free (Mapping *map, double dt, double mass, double wn, double wnc,
    double wnb, double zeta, int forc);
void update_free (Mapping *map, double dt, double mass, double wn, double wnc,
    double wnb, double zeta, int forc, double *Tension, double *Stiffness);
void s_MapField (Domain *omega, int dir, int KP);
// void filter (Map *map, char type); /* n-point linear filter (2-d only) */
void filter (Mapping *map, char type);

void AddAccelPBCs(Domain *omega);
void run_obj_vars(Domain *omega);
void s_run_obj_vars(Domain *omega, int KP);

/* TEST */
void Compare (Domain *omega);
void Compare (Domain *omega, ACTION space);
#ifdef TEST
void velocity (double, double, double, double, 
	       double*, double*, double*, double*);
void drive_force (double, double, double, double, 
		  double*, double*, double*);
void mapping (double, double, 
	      double*, double*, double*, double*, double*, double*);
void UpdateFFs (Domain *omega);
void compare_errors (Element *X, double *exact);
#endif
/* IO */

int  ReadMap  (Domain *omega); /* read the map file */
int s_ReadMap (Domain *omega, int k);

// void readMap  (FILE *fp, Map *map); /* read a single map */
void readMap  (FILE *fp, Mapping *map); /* read a single map */

void WriteMap (Domain *omega); /* write the map file (for restarts etc.) */
void s_WriteMap (Domain *omega, int KP);

// void writeMap (FILE *fp, Map *map); /* write a single map */
void writeMap (FILE *fp, Mapping *map);
// void readmap  (char *name, int nz, Map *mapx, Map *mapy);
void readmap  (char *name, int nz, Mapping *mapx, Mapping *mapy);


/* EXTRA */

void ErrorHandler(char *section, char *message, FLAG type);
double *Zmesh (int nz);

#if DIM == 3
void     dcopy_0all     (int np, double *x0, double *y);
void     dcopy_all0     (int np, double *y, double *x0);
void     dbroadcast     (int n, double *x);
void     dallgather     (int n, double *x, double *work);
void     dgather        (int n, double *x);
#endif




