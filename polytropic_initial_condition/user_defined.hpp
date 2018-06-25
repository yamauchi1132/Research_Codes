#pragma once

double dt = 0.05; //interval of particles in uniform sphere 
double radius = 1.0; //radius of uniform sphere
double pori_num = 1.5; //poritorope number
//double pori_num = 2.5; //poritorope number

double msun = 1.989e+33;
double rsun =  695700e+5;

double MASS = msun; //mass of real star
double RADIUS = rsun; //* pow(2.0, 1./3.); //radius of real star
double G = 6.67259e-8; //gravittational constant
/*
double MASS = 1.0; //mass of real star
double RADIUS = 1.0; //radius of real star
double G = 6.67259e-8; //gravittational constant
*/
double eps = 1e-6; 

//define class 
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
