#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

void calc_point_physics(Particle *p, long long int N, double *point_m, double (*point_pos)[2], double (*point_vel)[2]) {
	double row1_sum = 0;
	double row2_sum = 0;
	double pos_sum[2][2];
	double vel_sum[2][2];
	
	for(int i = 0; i < N/2; i++) {
		row1_sum += p[i].dens;
		pos_sum[0][0] += p[i].dens * p[i].pos[0];
		pos_sum[0][1] += p[i].dens * p[i].pos[1];
		vel_sum[0][0] += p[i].dens * p[i].vel[0];
		vel_sum[0][1] += p[i].dens * p[i].vel[1];
	}
	for(int i = N/2; i < N; i++) {
		row2_sum += p[i].dens;
		pos_sum[1][0] += p[i].dens * p[i].pos[0];
		pos_sum[1][1] += p[i].dens * p[i].pos[1];
		vel_sum[1][0] += p[i].dens * p[i].vel[0];
		vel_sum[1][1] += p[i].dens * p[i].vel[1];
	}
	
	point_pos[0][0] = pos_sum[0][0] / row1_sum;
	point_pos[0][1] = pos_sum[0][1] / row1_sum;
	point_pos[1][0] = pos_sum[1][0] / row2_sum;
	point_pos[1][1] = pos_sum[1][1] / row2_sum;

	point_vel[0][0] = vel_sum[0][0] / row1_sum;
	point_vel[0][1] = vel_sum[0][1] / row1_sum;
	point_vel[1][0] = vel_sum[1][0] / row2_sum;
	point_vel[1][1] = vel_sum[1][1] / row2_sum;

	for(int i = 0; i < N/2; i++) {
		double dx = p[i].pos[0]-point_pos[0][0];
		double dy = p[i].pos[1]-point_pos[0][1];
		double r_2 = dx*dx + dy*dy;
		double r = sqrt(r_2);
		if(r < 2*R1) {
			point_m[0] += p[i].mass;
		}
	}
	
	for(int i = N/2; i < N; i++) {
		double dx = p[i].pos[0]-point_pos[1][0];
		double dy = p[i].pos[1]-point_pos[1][1];
		double r_2 = dx*dx + dy*dy;
		double r = sqrt(r_2);
		if(r < 2*R2) {
			point_m[1] += p[i].mass;
		}
	}
	//fprintf(stdout , "%e %e %e %e %e %e\n", point_pos[0][0], point_pos[0][1], point_pos[1][0], point_pos[1][1], point_m[0], point_m[1]);
}

int main(int argc, char **argv) {

	if(argc != 2) {
		fprintf(stderr, "Error : no input file\n");
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
	quick_sort(p, 0, N-1);

	double point_m[2];
	double point_pos[2][2];
	double point_vel[2][2];
	calc_point_physics(p, N, point_m, point_pos, point_vel);

	for(int i = 0; i < N/2; i++) {
		fprintf(stdout , "%lld %lf %lf\n", p[i].id, p[i].pos[0], p[i].pos[1]);
	}

	return 0;
}