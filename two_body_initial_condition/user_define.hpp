#pragma once

double msun = 1.989e+33;
double rsun = 695700e+5;
double G = 6.67259e-8;

double m1 = msun;
double m2 = msun;
double r1 = rsun;
double r2 = rsun;

double r_mul = 1.5;
double r_rel_p = r_mul * (r1 + r2);
double v_rel_inf = 1e+6;
//double r_rel = 80 * 695700e+5;
//double r_rel = 20 * 695700e+5;
double r_rel = 10 * r1;
double pori_n = 1.5;
// double pori_n = 2.5;

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



