#pragma once

double msun = 1.989e+33;
double rsun = 695700e+5;
double G = 6.67259e-8;

//////// you need to set ////////
double m1 = 10*msun;
double r1 = 4*rsun;
double pori_n1 = 1.5;

double m2 = 10*msun;
double r2 = 4*rsun;
double pori_n2 = 1.5;

double r_mul = 1.0;
double v_rel_inf = 1e+6;
double r_rel = 10 * r1;
/////////////////////////////////

double r_rel_p = r_mul * (r1 + r2);

class Particle {
public:
  long long int id;
  long long int istar;
  double mass;
  double pos[3];
  double vel[3];
  double acc[3];
  double uene;
  double alph;
  double alphu;
  double dens;
  double ksr;
  double np;
  double vsnd;
  double pres;
  double emp;
  double divv;
  double rotv;
  double bswt;
  double pot;
  double abar;
  double zbar;
  double enuc;
  double vsmx;
  double udot;
  double dnuc;
  double cmps[18];
};



