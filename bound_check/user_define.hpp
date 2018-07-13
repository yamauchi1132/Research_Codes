#pragma once

double rsun = 695700e+5;
double G = 6.67259e-8;

//////// you need to set ////////
double R1 = 10*rsun;
double R2 = 10*rsun;
/////////////////////////////////

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
