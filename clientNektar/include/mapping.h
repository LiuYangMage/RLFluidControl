#ifndef _MAPPING_H
#define _MAPPING_H

typedef struct mppng {            /* ......... Mapping ............... */
  int       NZ                  ; /* Number of z-planes                */
  double    time                ; /* Current time                      */
  double   *d                   ; /* Displacement                      */
  double   *z                   ; /*   (z-deriv)                       */
  double   *zz                  ; /*   (zz-deriv)                      */
#ifdef NONLINEAR
  double   *zzz                  ; /*   (zzz-deriv)                      */
#endif
  double   *t                   ; /* Velocity                          */
  double   *tt                  ; /* Acceleration                      */
  double   *tz                  ; /*   (tz-deriv)                      */
  double   *tzz                 ; /*   (tzz-deriv)                     */
  double   *f                   ; /* Force                             */
//} Map;  be very careful on certain architecture (Trapezoid, XD1-Cray) there is a name conflict and there is already a class Map defined on the system. Therefore the name should be changed to another name.
} Mapping;

#endif
