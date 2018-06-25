#include <stdio.h>
#include <sys/stat.h>
#include <math.h>
#include <sysexits.h>
#include <stdlib.h>
#include "user_defined.hpp"
#include "lane_emden.hpp"
#include "make_and_output_data.hpp"

//make uniform sphere 
int make_uniform_sphere(Particle *p, double dt) 
{
  int num = 0;
  
  for(double i = -1; i < 1; i += dt) {
    for(double j = -1; j < 1; j += dt) {
      for(double k = -1; k < 1; k += dt) {
        //double r_p = (i-RADIUS)*(i-RADIUS) + (j-RADIUS)*(j-RADIUS) + (k-RADIUS)*(k-RADIUS);
        double r_p = i*i + j*j + k*k;
        if(r_p < radius*radius) {
          p[num].x = i;
          p[num].y = j;
          p[num].z = k;
          p[num].r = sqrt(r_p);
          num++;
        }
      }
    }
  }

  return num;
}

void quick_sort(Particle *p, int left, int right)
{
  int i = left;  
  int j = right;
  double pivot = p[(left+right)/2].r;

  while(1) {
    while (p[i].r < pivot) {
      i++;
    }

    while (pivot < p[j].r) {
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

void decide_mass_distribution(Particle *p, int p_num, LaneEmden *le, int lane_num)
{
//uniform mass distribution(cumulative mass distribution)
  double m = 1. / p_num;
  double sum_m = 0.;
  int i = 0;
  int j = 0;
  
  while(1) {
    if(i+1 == p_num) {
      sum_m += m;
      for(int k = j; k < i+1; k++) {
        p[k].m_total = sum_m;
      }
      break;
    }
    if(p[i].r == p[i+1].r) {
      sum_m += m;
    }
    if(p[i].r != p[i+1].r) {
      sum_m += m;
      for(int k = j; k < i+1; k++) {
        p[k].m_total = sum_m;
      }
      j = i + 1;
    }

    i++;
  }

//poritorope mass distribution(integraton over r)
  le[0].m_total = (4 * M_PI * le[0].r * le[0].r * le[0].r * le[0].rho) / 3;

  for(int i = 1; i < lane_num; i++) {
    //double dm = 4 * M_PI * le[i].r * le[i].r * (le[i].r - le[i-1].r) * le[i].rho;
    double dv = le[i].r * le[i].r * le[i].r - le[i-1].r * le[i-1].r * le[i-1].r;
    double dm = (4 * M_PI * dv * le[i].rho) / 3.;
    le[i].m_total = le[i-1].m_total + dm;
  }

//normalization
  double max_m = le[lane_num-1].m_total;
  //double max_r = le[lane_num-1].t;
  for(int i = 0; i < lane_num; i++) {
    le[i].m_total = le[i].m_total / max_m;
    //printf("%e %e\n", le[i].r, le[i].m_total);
    //le[i].t = le[i].t / max_r;
  }
}

void binary_search(Particle *p, int p_num, LaneEmden *le, int lane_num, int i)
{
  int left = 0;
  int right = lane_num;

  while(1) {
    int mid = (left + right) / 2;
    if(le[mid].m_total == p[i].m_total) {
      p[i].r_pori = le[mid].r;
      p[i].lane_id = le[mid].id;
      //printf("%lf %lf\n", p[i].m, le[mid].m);
    }
    if(le[mid].m_total < p[i].m_total) {
      left = mid + 1;
    } 
    if(le[mid].m_total > p[i].m_total) {
      right = mid - 1;
    }
    if(left >= right) {
      p[i].r_pori = le[mid].r;
      p[i].lane_id = le[mid].id;
      //printf("%lf %lf %lf %lf %lf\n", p[i].m, le[mid].m, p[i].r, le[mid].r, p[i].r_pori);
      break;
    }
  }
}

void expand_pori_position(Particle *p, int i)
{
  double ratio = p[i].r_pori / (p[i].r + eps);
  //printf("%lf\n", ratio);
  p[i].x_pori = p[i].x * ratio;
  p[i].y_pori = p[i].y * ratio;
  p[i].z_pori = p[i].z * ratio;
}

void convert_pori_quantities(Particle *p, int p_num, LaneEmden *le, int lane_num)
{
  for(int i = 0; i < p_num; i++) {
    binary_search(p, p_num, le, lane_num, i);
    expand_pori_position(p, i);
  }
}

void check_mass_distribution(Particle *p, int p_num)
{
  double m = 1. / p_num;
  double sum_m = 0.;
  int i = 0;
  int j = 0;

  while(1) {
    if(i+1 == p_num) {
      sum_m += m;
      for(int k = j; k < i+1; k++) {
        p[k].m_pori_total = sum_m;
      }
      break;
    }
    if(p[i].r_pori == p[i+1].r_pori) {
      sum_m += m;
    }
    if(p[i].r_pori != p[i+1].r_pori) {
      sum_m += m;
      for(int k = j; k < i+1; k++) {
        p[k].m_pori_total = sum_m;
      }
      j = i + 1;
    }

    i++;
  }
}

void expand_real_position(Particle *p, int i, double a)
{
  double r_real = p[i].r_pori * a;
  double ratio = r_real / (p[i].r_pori + eps);
  p[i].x_real = p[i].x_pori * ratio; 
  p[i].y_real = p[i].y_pori * ratio; 
  p[i].z_real = p[i].z_pori * ratio; 
  p[i].r_real = r_real;
}

void convert_real_quantities(Particle *p, int p_num, LaneEmden *le, int lane_num)
{
  double rho_ave = (3 * MASS) / (4 * M_PI * RADIUS * RADIUS * RADIUS);
  double rho_c = -(rho_ave) / ((3 * le[lane_num-1].v) / le[lane_num-1].r);
  double k = (4 * M_PI * G * RADIUS * RADIUS * pow(rho_c, (-1+pori_num)/pori_num))
  / ((pori_num+1) * le[lane_num-1].r * le[lane_num-1].r);
  double a = sqrt(((pori_num+1) * k * pow(rho_c, (1-pori_num)/pori_num)) / (4 * M_PI * G));

  fprintf(stderr, "k = %e\n", k);
  
  for(int i = 0; i < p_num; i++) {
    double term1 = -MASS * le[lane_num-1].r * pow(le[p[i].lane_id].rho, pori_num);
    double term2 = 4 * M_PI * RADIUS * RADIUS * RADIUS * le[lane_num-1].v;
    double rho = term1 / term2;
    double pres = k * pow(rho, 1+1/pori_num);

    expand_real_position(p, i, a);
    p[i].dens = rho;
    p[i].pres = pres;
    p[i].vsnd = sqrt(((1+(1/pori_num)) * pres) / rho);
    p[i].uene = pori_num * pres;
  }
}

int main(void)
{
  int n = 10e+6;
  Particle *p = new Particle[n];
  
//make uniform_sphere
  int p_num = make_uniform_sphere(p, dt);
  fprintf(stderr, "p_num = %d\n", p_num);
  quick_sort(p, 0, p_num-1);
  //add ID for particles
  for(int i = 0; i < p_num; i++) {
    p[i].id = i;
  }
  
//solve lane_emden equation
  LaneEmden *le = new LaneEmden[n];
  int lane_num = lane_emden(le, pori_num);
//make mass distribution of uniform sphere and poritirope  
  decide_mass_distribution(p, p_num, le, lane_num);

//mapping particle of uniform mass coordinate r to poritorope coordinate r'
  convert_pori_quantities(p, p_num, le, lane_num);
  check_mass_distribution(p, p_num);
 
//convert dimensionless quantities to dimensonal quantities
  convert_real_quantities(p, p_num, le, lane_num);

//making data file
  char dir[256];
  sprintf(dir, "n%dk_%.1lfMsun_%.1lfRsun_pori%.2lf", (int)p_num/10000, MASS/msun, RADIUS/rsun, pori_num);
  mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);
  make_data(p, p_num, dir);

  
//output file 
  char dir_name[256];
  sprintf(dir_name, "output_n_%.1f", pori_num);
  mkdir(dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
  output_sphere(p, p_num, dir_name);
  output_lane_emden(le, lane_num, dir_name);
  output_real(p, p_num, dir_name);


  delete[] p;
  delete[] le;

  return 0;
}
