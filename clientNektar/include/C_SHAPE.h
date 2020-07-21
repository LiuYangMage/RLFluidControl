#ifndef C_SHAPE_H
#define C_SHAPE_H


class C_SHAPE{

public:
  char type[32];
  double *data;
  double *bound_box;
  double ***BCKGRND;
  void read_geometry(char *file_name);
  void compute_indicator(int npt, int plane, double *x, double *y, double *indicator_value);
  void update_position(double dt);
  void update_velocity(double f_x, double f_y, double f_z, double dt);
  void get_velocity(double *u, double *v, double *w);
  void get_position(double *x, double *y, double *z);
  void update_position_and_velocity(double f_x, double f_y, double f_z, double dt);
};

struct C_VERT{
  double x;
  double y;
};
#endif
