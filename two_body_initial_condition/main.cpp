#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sysexits.h>
#include <sys/stat.h>
#include "user_define.hpp"

//Read file and input same infomation into array p[1:N] and p[N:2*N] in order to make two body 
void readfile(char *filename, Particle *p, long long int N) {
  FILE *fp = fopen(filename, "r");
  for(int i = 0; i < N; i++) {
    fscanf(fp, "%lld %lld %lf", &p[i].id, &p[i].istar, &p[i].mass); //3
    fscanf(fp, "%lf %lf %lf", &p[i].pos[0], &p[i].pos[1], &p[i].pos[2]); //6
    fscanf(fp, "%lf %lf %lf", &p[i].vel[0], &p[i].vel[1], &p[i].vel[2]);  //9
    fscanf(fp, "%lf %lf %lf", &p[i].acc[0], &p[i].acc[1], &p[i].acc[2]);  //12
    fscanf(fp, "%lf %lf %lf", &p[i].uene, &p[i].alph, &p[i].alphu);  //15
    fscanf(fp, "%lf %lf %lf", &p[i].dens, &p[i].ksr, &p[i].np);   //18
    fscanf(fp, "%lf %lf %lf", &p[i].vsnd, &p[i].pres, &p[i].emp);  //21
    fscanf(fp, "%lf %lf %lf", &p[i].divv, &p[i].rotv, &p[i].bswt);  //24
    fscanf(fp, "%lf %lf %lf", &p[i].pot, &p[i].abar, &p[i].zbar);  //27
    fscanf(fp, "%lf %lf %lf", &p[i].enuc, &p[i].vsmx, &p[i].udot); //30
    fscanf(fp, "%lf", &p[i].dnuc);  //31
    for(int i = 0; i < 18; i++) {  //32-49
      fscanf(fp, "%lf", &p[i].cmps[i]);
    }
  }

  FILE *fp2 = fopen(filename, "r");
  for(int i = N; i < 2*N; i++) {
    fscanf(fp2, "%lld %lld %lf", &p[i].id, &p[i].istar, &p[i].mass); //3
    fscanf(fp2, "%lf %lf %lf", &p[i].pos[0], &p[i].pos[1], &p[i].pos[2]); //6
    fscanf(fp2, "%lf %lf %lf", &p[i].vel[0], &p[i].vel[1], &p[i].vel[2]);  //9
    fscanf(fp2, "%lf %lf %lf", &p[i].acc[0], &p[i].acc[1], &p[i].acc[2]);  //12
    fscanf(fp2, "%lf %lf %lf", &p[i].uene, &p[i].alph, &p[i].alphu);  //15
    fscanf(fp2, "%lf %lf %lf", &p[i].dens, &p[i].ksr, &p[i].np);   //18
    fscanf(fp2, "%lf %lf %lf", &p[i].vsnd, &p[i].pres, &p[i].emp);  //21
    fscanf(fp2, "%lf %lf %lf", &p[i].divv, &p[i].rotv, &p[i].bswt);  //24
    fscanf(fp2, "%lf %lf %lf", &p[i].pot, &p[i].abar, &p[i].zbar);  //27
    fscanf(fp2, "%lf %lf %lf", &p[i].enuc, &p[i].vsmx, &p[i].udot); //30
    fscanf(fp2, "%lf", &p[i].dnuc);  //31
    for(int i = 0; i < 18; i++) {  //32-49
      fscanf(fp2, "%lf", &p[i].cmps[i]);
    }
  }
}

//Calculation of impact parameter b
double calc_impact_parameter(double k_e_inf) {
  //double mass_rel = m1 * m2 / (m1 + m2);
  double p_e_p = -(G * (m1 + m2)) / r_rel_p;
  double v_rel_p = sqrt((k_e_inf - p_e_p) * 2.0);
  double b = (r_rel_p * v_rel_p) / v_rel_inf;
  return b;
}

//Calculation of eccentricity e
double calc_eccentricity(double b, double k_e_inf) {
  double ang_mom2 = b * b * v_rel_inf * v_rel_inf;
  double term1 = (2 * k_e_inf * ang_mom2) / ((G*(m1+m2)) * (G*(m1+m2)));
  double e = sqrt(1 + term1);
  return e;
}

//Caluculate relative position and velocity of particle1
void calc_relative_position_and_velocity(double b, double e, double *pos1_rel, double *vel1_rel, double k_e_inf) {
  double a = (G * (m1 + m2)) / (2 * k_e_inf);
  double term1 = (a * (e*e - 1)) / r_rel;
  double term2 = term1 - 1;
  double cos_f = term2 / e;
  double sin_f = -sqrt(1 - cos_f*cos_f);
  
  pos1_rel[0] = r_rel * cos_f;
  pos1_rel[1] = r_rel * sin_f;

  double term3 = (G * (m1 + m2)) / (a * a * a);
  double n = sqrt(term3);
  double term4 = (a * n) / sqrt(e*e - 1);

  vel1_rel[0] = -term4 * sin_f;
  vel1_rel[1] = term4 * (cos_f + e);
}

//Caluculate origingal position and velocity of particle1 and particle2
void calc_original_position_and_velocity(double *pos1_rel, double *vel1_rel, double *pos1, double *pos2, double *vel1, double *vel2) {
  double term1 = m2 / (m1 + m2);
  double term2 = -m1 / (m1 + m2);
  pos1[0] = term1 * pos1_rel[0];
  pos1[1] = term1 * pos1_rel[1];
  pos2[0] = term2 * pos1_rel[0];
  pos2[1] = term2 * pos1_rel[1];

  vel1[0] = term1 * vel1_rel[0];
  vel1[1] = term1 * vel1_rel[1];
  vel2[0] = term2 * vel1_rel[0];
  vel2[1] = term2 * vel1_rel[1];
}

//Output data file and header file
void output_file(Particle *p, long long int N, double *pos1, double *pos2, double *vel1, double *vel2) {
  for(int i = 0; i < N; i++) {
    p[i].pos[0] += pos1[0];
    p[i].pos[1] += pos1[1];
    p[i].vel[0] += vel1[0];
    p[i].vel[1] += vel1[1];

    //p[i].pos[2] = 0.;
    //p[i].vel[2] = 0.;
    // fprintf(stdout, "%e %e\n", p[i].pos[0], p[i].pos[1]);
  }
  
  for(int i = N; i < 2*N; i++) {
    p[i].id = p[i].id + N;
    p[i].pos[0] += pos2[0];
    p[i].pos[1] += pos2[1];
    p[i].vel[0] += vel2[0];
    p[i].vel[1] += vel2[1];

    //p[i].pos[2] = 0.;
    //p[i].vel[2] = 0.;
    // fprintf(stdout, "%e %e\n", p[i].pos[0], p[i].pos[1]);
  }

  char dir[256];
  sprintf(dir, "data_n%dk_m%.2e_m%.2e_rp%.2lf*(%.2e+%.2e)_vinf%.2e_pori%.1lf", (int)N / 10000, m1, m2, r_rel_p/(r1+r2), r1, r2,  v_rel_inf, pori_n);
  mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);

  FILE *data;
  FILE *head;
  char data_name[256];
  char head_name[256];

  sprintf(data_name, "%s/two_body.00_t000.data", dir);
  data = fopen(data_name, "w");

  sprintf(head_name, "%s/two_body.00_t000.head", dir);
  head = fopen(head_name, "w");

  for(int i = 0; i < 2*N; i++) {
    fprintf(data, "%lld %lld %e %e %e %e %e %e %e %e %e %e %e",
	    p[i].id, p[i].istar, p[i].mass, p[i].pos[0], //4
	    p[i].pos[1], p[i].pos[2], p[i].vel[0], p[i].vel[1], //8
	    p[i].vel[2], p[i].uene, p[i].alph, p[i].alphu, p[i].ksr); //13

    for(int k = 0; k < 13; k++) {
      p[i].cmps[k] = (k==1 || k==2) ? 5.0e-1 : 0;
      fprintf(data, " %e", p[i].cmps[k]);  //14~26
    }
    fprintf(data, "\n");
  }


  //make header firle
  fprintf(head, "+0.000000e+00 +1.0e+4 1.0e+3 1.0e+3\n");
  fprintf(head, "0 0\n");
  fprintf(head, "+2.000000e+00 +1.000000e-01\n");
  fprintf(head, "+0.000000e+00 +0.000000e+00\n");
  fprintf(head, "+1.0000000e-01\n");
  fprintf(head, "%lld\n", 2*N);
}

//Make two body itnitial condition.
void make_initial_condition(Particle *p, long long int N)
{
  double k_e_inf = 0.5 * v_rel_inf * v_rel_inf;

  double b = calc_impact_parameter(k_e_inf);  
  double e = calc_eccentricity(b, k_e_inf);
  // fprintf(stderr, "b,e,k_e_inf : %e %e %e\n\n", b, e, k_e_inf);

  double pos1_rel[2], vel1_rel[2];
  calc_relative_position_and_velocity(b, e, pos1_rel, vel1_rel, k_e_inf);
  fprintf(stderr, "relative pos and vel : %e %e %e %e\n\n", pos1_rel[0], pos1_rel[1], vel1_rel[0], vel1_rel[1]);

  double pos1[2], pos2[2], vel1[2], vel2[2];
  calc_original_position_and_velocity(pos1_rel, vel1_rel, pos1, pos2, vel1, vel2);
  fprintf(stderr, "pos and vel : %e %e %e %e %e %e %e %e\n", pos1[0], pos1[1], pos2[0], pos2[1], vel1[0], vel1[1], vel2[0], vel2[1]);
  // fprintf(stderr, "%e %e\n", 2*sqrt(pos1[0]*pos1[0] + pos2[1]*pos2[1]), r_rel);

  output_file(p, N, pos1, pos2, vel1, vel2);
}

int main(int argc, char **argv)
{
  if(argc != 2) {
    fprintf(stderr, "Error : Input file\n");
    exit(1);
  }

  //File processing and count lines to decide particle number N.
  FILE *fp = fopen(argv[1], "r");
  long long int N = 0;
  if(fp == false) {
    fprintf(stderr, "Error : Don't open the file\n");
  } else {
    char buff[1024];
    while(fgets(buff, 1024, fp) != NULL) {
      N++;
    }
    fprintf(stderr, "particle number : %lld\n", N);
  }

  Particle *p;
  p = new Particle[2*N];

  //Read file into array p
  readfile(argv[1], p, N);

  //Make initial condition
  make_initial_condition(p, N);

  delete[] p;
  return 0;
}
