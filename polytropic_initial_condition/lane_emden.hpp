#pragma once

double f1(double r,double rho,double v)
{
  return v;
}

double f2(double r,double rho,double v)
{
  //printf("%lf %lf %lf %lf\n", -pow(rho, pori_num),rho, 2.0/r*v,-pow(rho,pori_num)-2.0/r*v);
  return -pow(rho,pori_num) - 2.0/r*v;
}

int lane_emden(LaneEmden *le, double pori_num)
{
  double k1[2],k2[2],k3[2],k4[2];

  double r0=0.0001;      //初期値
  double rho = 1.0 - 1.0/6.0*r0*r0 + pori_num/120.0*r0*r0*r0*r0;        //初期値
  double v = -1.0/3.0*r0 + pori_num/30.0*r0*r0*r0;  //初期値
  double dr = 0.001;       //刻み幅
  double rmax = 10000.0;   //最大値

  int lane_num = 0;
  for(double r = r0; r < rmax; r += dr, lane_num++) {
    k1[0] = dr * f1(r,rho,v);
    k1[1] = dr * f2(r,rho,v);

    if(rho+k1[0]/2.0 < 0.) break;
    k2[0] = dr * f1(r+dr/2.0, rho+k1[0]/2.0, v+k1[1]/2.0);
    k2[1] = dr * f2(r+dr/2.0, rho+k1[0]/2.0, v+k1[1]/2.0);

    if(rho+k2[0]/2.0 < 0.) break;
    k3[0] = dr * f1(r+dr/2.0, rho+k2[0]/2.0, v+k2[1]/2.0);
    k3[1] = dr * f2(r+dr/2.0, rho+k2[0]/2.0, v+k2[1]/2.0);

    if(rho+k3[0]/2.0 < 0.) break;
    k4[0] = dr * f1(r+dr, rho+k3[0], v+k3[1]);
    k4[1] = dr * f2(r+dr, rho+k3[0], v+k3[1]);

    //printf("%e %e %e %e %e %e\n", rho, (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]) / 6.0, k1[1], k2[1], k3[1], k4[1]);
    rho = rho + (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]) / 6.0;
    v = v + (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]) / 6.0;
    le[lane_num].id = lane_num;
    le[lane_num].r = r;
    le[lane_num].rho = rho;
    le[lane_num].v = v;
    //fprintf(stdout,"%f %f %f %d\n",le[lane_num].r,le[lane_num].rho,le[lane_num].v, lane_num);
    if(le[lane_num].rho < 0.) break;
  }

  return lane_num;
}