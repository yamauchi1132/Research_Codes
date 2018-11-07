#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sysexits.h>
#include <sys/stat.h>
#include "class.hpp"
#include "user_define.hpp"

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

void making_data_simulation(Particle *p, long long int N) {
  char dir[256];
  sprintf(dir, "n%dk_%.1lfMsun_%.1lfRsun_NewPori%.1lf_PolyChange", (int)N / 10000, m_mul, r_mul, new_pori);
  mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);

  FILE *data;
  FILE *head;
  char data_name[256];
  char head_name[256];

  sprintf(data_name, "%s/final.data", dir);
  sprintf(head_name, "%s/final.head", dir);

  data = fopen(data_name, "w");
  head = fopen(head_name, "w");

  for(int i = 0; i < N; i++) {
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
  fprintf(head, "+0.000000e+00 +1.0e+6 1.0e+4 1.0e+4\n");
  fprintf(head, "0 0\n");
  fprintf(head, "+2.000000e+00 +1.000000e-01\n");
  fprintf(head, "+0.000000e+00 +0.000000e+00\n");
  fprintf(head, "+1.0000000e-01\n");
  fprintf(head, "%lld\n", N);
}

void making_data(Particle *p, long long int N) {
  char dir[256];
  sprintf(dir, "snap_n%dk_%.1lfMsun_%.1lfRsun_NewPori%.1lf_PolyChange", (int)N / 10000, m_mul, r_mul, new_pori);
  mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);

  FILE *data;
  char data_name[256];
  sprintf(data_name, "%s/final.dat", dir);
  data = fopen(data_name, "w");

  for(int i = 0; i < N; i++) {
    fprintf(data, "%lld %lld %lf", p[i].id, p[i].istar, p[i].mass); //3
    fprintf(data, " %lf %lf %lf", p[i].pos[0], p[i].pos[1], p[i].pos[2]); //6
    fprintf(data, " %lf %lf %lf", p[i].vel[0], p[i].vel[1], p[i].vel[2]);  //9
    fprintf(data, " %lf %lf %lf", p[i].acc[0], p[i].acc[1], p[i].acc[2]);  //12
    fprintf(data, " %lf %lf %lf", p[i].uene, p[i].alph, p[i].alphu);  //15
    fprintf(data, " %lf %lf %lf", p[i].dens, p[i].ksr, p[i].np);   //18
    fprintf(data, " %lf %lf %lf", p[i].vsnd, p[i].pres, p[i].emp);  //21
    fprintf(data, " %lf %lf %lf", p[i].divv, p[i].rotv, p[i].bswt);  //24
    fprintf(data, " %lf %lf %lf", p[i].pot, p[i].abar, p[i].zbar);  //27
    fprintf(data, " %lf %lf %lf", p[i].enuc, p[i].vsmx, p[i].udot); //30
    fprintf(data, " %lf", p[i].dnuc);  //31
    for(int i = 0; i < 18; i++) {  //32-49
      fprintf(data, " %lf", p[i].cmps[i]);
    }
    fprintf(data, "\n");
  }
}

void change_poly_num(Particle *p, long long int N) {
    for(int i = 0; i < N; i++) {
    // double new_uene = p[i].pres / ((new_ganma - 1) * p[i].dens);
    double ratio = (old_ganma - 1) / (new_ganma - 1);
    double new_uene = ratio * p[i].uene;
    //printf("%e %e %e %e %e\n", p[i].pos[0], p[i].pos[1], p[i].pos[2], p[i].uene, new_uene);
    p[i].uene = new_uene;
  }
}
/*
void change_poly_num(Particle *p, long long int N) {
  for(int i = 0; i < N; i++) {
    double new_k = p[i].pres / pow(p[i].dens, new_ganma);
    double new_uene = (new_k / new_ganma-1) * pow(p[i].dens, new_ganma-1);
    //printf("%e %e %e %e %e\n", p[i].pos[0], p[i].pos[1], p[i].pos[2], p[i].uene, new_uene);
  }
}
*/
int main(int argc, char **argv) {
  if(argc != 2) {
    fprintf(stderr, "Error : No input file\n");
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
  p = new Particle[N];
  readfile(argv[1], p, N);

  change_poly_num(p, N);

  //making_data_simulation(p, N);
  making_data(p, N);

  return 0;
}