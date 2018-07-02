#pragma once

//////// you need to set ////////
double m_mul = 10.0;
double r_mul = 4.0;
/////////////////////////////////

double v_mul = sqrt(m_mul/r_mul);

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