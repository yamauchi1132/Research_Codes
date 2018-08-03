#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "class.hpp"

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
}

void quick_sort(Particle *p, int left, int right)
{
  int i = left;  
  int j = right;
  double pivot = p[(left+right)/2].id;

  while(1) {
    while (p[i].id < pivot) {
      i++;
    }

    while (pivot < p[j].id) {
      j--;
    }

    if (i >= j) {
      break;
    }

    Particle tmp = p[i];
    p[i] = p[j];
    p[j] = tmp;

    i++;
    j--;
  }

  if (left < i-1) {
    quick_sort(p, left, i-1);
  }
  if (j+1 < right) {
    quick_sort(p, j+1, right);
  }
}

double calc_potential_energy(Particle *p, long long int s, long long int e) {
  double p_e = 0.;
  for(i = s; i < e-1; i++) {
    for(j = i+1; j < e; j++){
      double x_ij = p[i].pos[0] - p[j].pos[0];
      double y_ij = p[i].pos[1] - p[j].pos[1];
      double z_ij = p[i].pos[2] - p[j].pos[2];
      double r2_ij = x_ij*x_ij + y_ij*y_ij + z_ij*z_ij;
      double r_ij = sqrt(r2_ij);
      p_e += (p[i].m * p[j].m) / (r_ij + 10e-6);
    }
  }
  return p_e;
}

void energy_check(Particle *p_ini, Particle *p_fin, long long int N) {
  calc_potential_energy(p_ini, 0, N/2);
}

int main(int argc, char **argv) {
  if(argc != 3) {
    fprintf(stderr, "Error : no input file\n");
  }

  //File processing and count lines to decide particle number N.
  FILE *fp1 = fopen(argv[1], "r");
  FILE *fp2 = fopen(argv[2], "r");
  long long int N = 0;
  long long int N2 = 0;
  if(fp1 == false || fp1 == false) {
    fprintf(stderr, "Error : Don't open the file\n");
  } else {
    char buff[1024];
    while(fgets(buff, 1024, fp1) != NULL) {
      N++;
    }
    while(fgets(buff, 1024, fp2) != NULL) {
      N2++;
    }
    if(N == N2) {
      fprintf(stderr, "particle number : %lld\n", N);
    } else {
      fprintf(stderr, "particle number is different\n");
    }
  }

  Particle *p_ini, *p_fin;
  p_ini = new Particle[N];
  p_fin = new Particle[N];

  readfile(argv[1], p_ini, N);
  readfile(argv[2], p_fin, N);

  quick_sort(p_ini, 0, N-1);
  quick_sort(p_ini, 0, N-1);
  /*
  for(int i = 0; i < N/2; i++) {
    fprintf(stdout, "%e %e\n", p1[i].pos[0], p2[i].pos[1]);
    fprintf(stdout, "%e %e\n", p1[i].pos[0], p2[i].pos[1]);
  }
  */
  energy_check(p_ini, p_fin, N);

}