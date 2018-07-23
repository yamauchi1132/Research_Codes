#pragma once

class Particle
{
public:
  int id;
  double x;
  double y;
  double z;
  double r;
  double m_total;

  int lane_id;
  double x_pori;
  double y_pori;
  double z_pori;
  double r_pori;
  double m_pori_total;
  
  long long int istar;
  double mass;
  double x_real;
  double y_real;
  double z_real;
  double r_real;
  double v[3];
  double pres;
  double dens;
  double vsnd;
  double uene;
  double alph;
  double alphu;
  double ksr;
  double cmps[13];
};

class LaneEmden
{
public:
  int id;
  double r; 
  double rho; 
  double v; //drow/dr
  double m_total; 
};
