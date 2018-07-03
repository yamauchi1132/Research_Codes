#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sysexits.h>
#include <sys/stat.h>
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

void change_scale(Particle *p, long long int N) {
  for(int i = 0; i < N; i++) {
    p[i].mass *= m_mul;

    p[i].pos[0] *= r_mul;
    p[i].pos[1] *= r_mul;
    p[i].pos[2] *= r_mul;

    p[i].vel[0] *= v_mul;
    p[i].vel[1] *= v_mul;
    p[i].vel[2] *= v_mul;

    p[i].uene *= uene_mul;

    p[i].ksr *= r_mul;
  }

  char dir[256];
  sprintf(dir, "n%dk_%.1lfMsun_%.1lfRsun_pori%.1lf_ScaleChange", (int)N / 10000, m_mul, r_mul, pori_num);
  mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);

  FILE *data;
  char data_name[256];
  sprintf(data_name, "%s/final.dat", dir);
  data = fopen(data_name, "w");

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
}

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

  change_scale(p, N);
  
  return 0;
}