/**************************************************
 ****************************
template.c: template file for practice 1
***************************************************
***************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "cpgplot.h"
#define NMAX 16384
#define SQR(x) ((x)*(x))
#define WSIZE 6.0 /* size of window in inch */
#define VMAX 1.25 /* limit of variables in window */
#define XPLOT 0   /* coordinate for x */
#define YPLOT 1   /* coordinate for y */
#define TX  0.9   /* x coordinate of time */
#define TY -1.1   /* y coordinate of time */
#define CFG 5     /* foreground color */
#define CBG 0     /* background color */


float make_spherical_df(int n,  float *m, float **x, float **v,float r_v, float eps2);
float gaussian(void);
void calc_force(int n, float *m, float **x, float **a, float eps2);
void leap_frog(int n, float *m, float **x, float **v, float **a,float dt, float eps2,float t);
void open_window();
void clse_window();

float make_spherical_df(int n, float *m, float **x, float **v,float r_v, float eps2)
{
  int i,j;
  float a, b, c, r;
  float w = 0.0;
  float k_e_0 = 0.0;
  float v_rms;
  float mass = 1.0;
  float e_ini = 0.0;
  //粒子分布の作成  
  //printf( "test\n");
  // printf("%e\n", r_v);

  for(i=0; i<n; i++){
    do{
      a=2.0*drand48()-1.0;
      b=2.0*drand48()-1.0;
      c=2.0*drand48()-1.0;
      r=sqrt(a*a+b*b+c*c);
    } while(r >= 1.0); //半径1の球の内部に存在しないときループする
    x[i][0] = a;
    x[i][1] = b;
    x[i][2] = c;
  }
  
  //printf( "test\n");
  
  //粒子の質量
  for(i=0; i<n; i++){
    m[i] = mass/n;
  }

  //fprintf(stderr,"test %e\n",m[1]);
  //ポテンシャルエネルギーの計算
  for(i=0; i<n-1; i++){
    for(j=i+1; j<n; j++){
      w += (m[i]*m[j]) / sqrt(pow(x[i][0]-x[j][0],2) + pow(x[i][1]-x[j][1],2) + pow(x[i][2]-x[j][2],2) + eps2*eps2);
    }
  }
  
  //速度分散の計算
  v_rms = sqrt((2*r_v*w) / (3*mass));


  //速度分布の作成
  for(i=0; i<n; i++){
    for(j=0; j<3; j++){
      v[i][j] = v_rms*gaussian();
    }
  }

  //運動エネルギーの計算
  for(i=0;i<n;i++){
    k_e_0 += 0.5*m[i]*(pow(v[i][0],2)+pow(v[i][1],2)+pow(v[i][2],2));
  } 

  //初期エネルギーの計算
  e_ini = k_e_0 - w;

  return e_ini;
  
}


//粒子分布を作るのに必要な関数
float gaussian(void)
/* Gaussian with mean = 0.0 and dispersion = 1.0 by Box-Muller   method */
{
  float x, y, r2;
  float z;
  do{
    x = 2.0*drand48() - 1.0;
    y = 2.0*drand48() - 1.0;
    r2 = x*x + y*y;
  }while(r2 >= 1.0 || r2 == 0.0);
  z = sqrt(-2.0*log(r2)/r2)*x; /* discard another Gaussian (
				  ....*y) */
  return z ;
}




void calc_force(int n, float *m, float **x, float **a, float eps2)
{
  int i, j, k;
  float r[3];
  float r3inv;
  //加速度の初期化
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      a[i][k] = 0.0;
    }
  }
  //相互重力の計算
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      //相対距離r[3]、rの3乗の逆数の計算
      r[0] = x[j][0]-x[i][0];
      r[1] = x[j][1]-x[i][1];
      r[2] = x[j][2]-x[i][2];
      r3inv = 1.0 / (pow(sqrt(r[0]*r[0] +r[1]*r[1]+ r[2]*r[2]+eps2*eps2),3));
      for(k=0;k<3;k++){
	a[i][k] += m[j]*r[k]*r3inv;
	a[j][k] += -m[i]*r[k]*r3inv;
      }
    }
  }
}


void leap_frog(int n, float *m, float **x, float **v, float **a,float dt, float eps2, float t)
{
  int i, k;
  float v_half[NMAX][3], x_start[NMAX][3];
  //t=0のときx0からa0を計算する
  if(t==0.0){
    calc_force(n, m, x, a, eps2);
  }

  /*
    for(i=0;i<n;i++){
    fprintf(stdout,"%e %e %e\n", a[i][0], a[i][1], a[i][2]);
    }
    exit(0);
  */

  //配列の初期化
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      x_start[i][k] = x[i][k];
    }
  }
  
  //時間積分
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      v_half[i][k] = v[i][k] + a[i][k]*dt*0.5;
      x[i][k] = x_start[i][k] + v_half[i][k]*dt;
    }
  }
  calc_force(n, m, x, a, eps2);
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      v[i][k] = v_half[i][k] + a[i][k]*dt*0.5;
    }
  }
}   


void movie(FILE *gnuplot, float **x, int n)
{
  int i;
  fprintf(gnuplot,"set term x11\n");
  fprintf(gnuplot,"set size square\n");
  fprintf(gnuplot,"set xrange [-1:1]\n");
  fprintf(gnuplot,"set yrange [-1:1]\n");
  fprintf(gnuplot,"set zrange [-1:1]\n");
  fprintf(gnuplot, "splot '-'\n");
  for(i=0;i<n;i++){
    fprintf(gnuplot,"%e %e %e\n", x[i][0], x[i][1], x[i][2]);
  }
  fprintf(gnuplot, "e\n");

}

int main(int argc, char** argv)
{

  /* number of particles */
  int n;
  /* particle data */
  static float *m, **x, **v, **a;
  /* system energy, Virial ratio */
  float r_v;
  /* squared softening parameter */
  float eps2;
  //time, timestep, end time, data interval
  float t, dt, t_end, t_out;
  float k_e, p_e, e, d_e, r_e, e_ini;
  //float x2[3][NMAX];
  t = 0.0;
  e_ini = 0;
  
  int i,j;

  //粒子数の入力
  fprintf(stderr, "n = ");
  scanf("%d", &n);
  //ビリアル比の入力
  fprintf(stderr, "r_v = ");
  scanf("%e", &r_v);
  //ソフトニングパラメーターの入力
  fprintf(stderr, "eps2 = ");
  scanf("%e", &eps2);
  //タイムステップの入力
  fprintf(stderr, "dt = ");
  scanf("%e", &dt);
  //シュミレーション終了時刻の入力
  fprintf(stderr, "t_end = ");
  scanf("%e", &t_end);
  //データ解析/出力の時間間隔の入力
  fprintf(stderr, "t_out = ");
  scanf("%e", &t_out);
  //確認
  fprintf(stderr, "n = %d, r_v = %e, eps2 = %e, dt = %e, t_end = %e, t_out = %e\n",n, r_v,eps2, dt, t_end, t_out);
  


  m=(float*)malloc(sizeof(int)*n);


  x=(float**)malloc(sizeof(float *)*n);
  v=(float**)malloc(sizeof(float *)*n);
  a=(float**)malloc(sizeof(float *)*n);
  for(i=0;i<n;i++){
    x[i]=(float*)malloc(sizeof(float)*3);
    v[i]=(float*)malloc(sizeof(float)*3);
    a[i]=(float*)malloc(sizeof(float)*3);
  }

  //初期分布の設定
  e_ini = make_spherical_df(n, m, x, v, r_v, eps2); 
  // 確認のための粒子分布の出力
  // n1024.snpファイルに出力
  /*
    for(i=0; i<n; i++){
    fprintf(stdout,"%e %e %e\n", x[i][0], x[i][1], x[i][2]);
    }
  */

  //ウィンドウを開く
  FILE* gnuplot=popen("gnuplot", "w");

  //open_window();


  clock_t start = clock();
  while(t<t_end){
    // realtime analysis 
    if(fmod(t, t_out) == 0){
      //運動エネルギーの計算  
      for(i=0;i<n;i++){
	k_e += 0.5*m[i]*(v[i][0]*v[i][0] + v[i][1]*v[i][1]+ v[i][2]*v[i][2]);
      } 
      //ポテンシャルエネルギーの計算
      for(i=0;i<n-1;i++){
	for(j=i+1;j<n;j++){
	  p_e += (m[i]*m[j]) / sqrt((x[i][0]-x[j][0])*(x[i][0]-x[j][0]) + (x[i][1]-x[j][1])*(x[i][1]-x[j][1]) + (x[i][2]-x[j][2])*(x[i][2]-x[j][2]) + eps2*eps2);
	}
      }
     
      r_v = k_e/p_e;
    }

    //time integration 
    leap_frog(n, m, x, v, a, dt, eps2,t);

    movie(gnuplot, x, n);

    t += dt;
    k_e = 0.0;
    p_e = 0.0;
  }
  clock_t end = clock();
  fprintf(stderr, "all-time = %lf seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
  pclose(gnuplot);

  /* integration error check */
  //  エネルギー誤差の確認
  //運動エネルギーの計算  
  for(i=0;i<n;i++){
    k_e += 0.5*m[i]*(pow(v[i][0],2)+pow(v[i][1],2)+pow(v[i][2],2));
  } 

  // fprintf(stderr, "%e", k_e);
  //ポテンシャルエネルギーの計算
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      p_e += (m[i]*m[j]) / sqrt(pow(x[i][0]-x[j][0],2) + pow(x[i][1]-x[j][1],2) + pow(x[i][2]-x[j][2],2) + eps2*eps2);
    }
  }

  //系のエネルギーEの計算
  e = k_e - p_e;
  //ΔEの計算
  d_e = e - e_ini;
  //相対誤差の計算
  r_e = fabs(d_e/e_ini);

  //  fprintf(stderr, "%f\t%f\n", dt,r_e);

  //t=10の時の粒子分布の表示
  /*
    for(i=0; i<n; i++){
    fprintf(stdout,"%e %e %e\n", x[i][0], x[i][1], x[i][2]);
    }
  */

  free(m);

  for(i=0;i<n;i++){
    free(x[i]);
    free(v[i]);
    free(a[i]);
  }
  free(x);
  free(v);
  free(a);

  return(0);
}


