#pragma once

double C = 10;

void decide_particle_info(Particle *p, int i, double mass)
{
  p[i].istar = 0;
  p[i].mass = mass;
  p[i].v[0] = 0.;
  p[i].v[1] = 0.;
  p[i].v[2] = 0.;
  p[i].alph = 1.;
  p[i].alphu = 0.;
  p[i].ksr = C * pow((mass / p[i].dens), 1./3.);
  for(int k = 0; k < 13; k++) {
    p[i].cmps[k] = (k==1 || k==2) ? 5.0e-1 : 0;
  }
}

void make_data(Particle *p, int p_num, char *dir)
{
  FILE *data;
  FILE *head;
  char data_name[256];
  char head_name[256];

  sprintf(data_name, "%s/poly.00_t000.data", dir);
  data = fopen(data_name, "w");

  sprintf(head_name, "%s/poly.00_t000.head", dir);
  head = fopen(head_name, "w");

//make data file
  double mass = MASS / p_num;
  for(int i = 0; i < p_num; i++) {
    decide_particle_info(p, i, mass);
    fprintf(data, "%d %lld %e %e %e %e %e %e %e %e %e %e %e",
     p[i].id, p[i].istar, p[i].mass, p[i].x_real, //4
     p[i].y_real, p[i].z_real, p[i].v[0], p[i].v[1], //8
     p[i].v[2], p[i].uene, p[i].alph, p[i].alphu, p[i].ksr); //13

    for(int k = 0; k < 13; k++) {
      fprintf(data, " %e", p[i].cmps[k]); //14~26
    }
    fprintf(data, "\n");
  }

//make head file
  fprintf(head, "+0.000000e+00 +50000.000000e+00 1000.000000e+00 1000.000000e+00\n");
  fprintf(head, "0 0\n");
  fprintf(head, "+2.000000e+00 +1.000000e-01\n");
  fprintf(head, "+0.000000e+00 +0.000000e+00\n");
  fprintf(head, "+1.0000000e-01\n");
  fprintf(head, "%d\n", p_num-1);
}

void output_sphere(Particle *p, int p_num, char *dir_name)
{
  FILE *sphere;
  char filename[256];
  sprintf(filename, "%s/sphere.data", dir_name);
  sphere = fopen(filename, "w");

  for(int i = 0; i < p_num; i++) {
    fprintf(sphere, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
     p[i].id, p[i].x, p[i].y, p[i].z, p[i].r, //5 
     p[i].m_total, p[i].x_pori, p[i].y_pori,  //8
     p[i].z_pori, p[i].r_pori, p[i].m_pori_total, p[i].lane_id); //12
  }
}

void output_lane_emden(LaneEmden *le, int lane_num, char *dir_name)
{
  FILE *lane_emden;
  char filename[256];
  sprintf(filename, "%s/lane_emden.data", dir_name);
  lane_emden = fopen(filename, "w");
  
  for(int i = 0; i < lane_num; i++) {
    fprintf(lane_emden, "%d %lf %lf %lf %lf\n", 
    le[i].id, le[i].r, le[i].rho, le[i].v, le[i].m_total); //5
  }
}

void output_real(Particle *p, int p_num, char *dir_name)
{
  FILE *real;
  char filename[256];
  sprintf(filename, "%s/real.data", dir_name);
  real = fopen(filename, "w");

  for(int i = 0; i < p_num; i++) {
   fprintf(real, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
   p[i].id, p[i].x_real, p[i].y_real, p[i].z_real, p[i].r_real, //5
   p[i].dens, p[i].pres, p[i].vsnd, p[i].uene);//9
 }
}
