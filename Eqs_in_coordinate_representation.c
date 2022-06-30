#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <omp.h>
#include <pthread.h>

/* Defining global constants*/
float PI=3.14159265358979323846;
double at=.1;
double ax=.5;
int Nt=200;
int Nf=15;
int Nx=21; /*always keep odd number of lattice point to have symmetric lattice around zero */
float A=3;
int r=0;
double u=0.5;

/* Define sigmaF*/
double sigmaF11(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((pow(f11[i][j][k][l],2) + f12[i][j][k][l]*f12[i][j][k][l] + f21[i][j][k][l]*f21[i][j][k][l]+pow(f22[i][j][k][l],2)- 1/4*(pow(rho11[i][j][k][l],2) + rho12[i][j][k][l]*rho12[i][j][k][l]+ rho21[i][j][k][l]*rho21[i][j][k][l]+pow(rho22[i][j][k][l],2)))*phi1[i][k]*phi1[j][l] +
           2*(f11[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                           f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l])-
              1/4*(rho11[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                              rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l]) )+

           f11[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
           f21[i][j][k][l]*(f21[i][j][k][l]*phi1[i][k]*phi1[j][l] + f22[i][j][k][l]*phi1[i][k]*phi2[j][l])-
      1/4*(rho11[i][j][k][l]* (rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      rho21[i][j][k][l]*(rho21[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi1[i][k]*phi2[j][l]))+

      f11[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      f12[i][j][k][l]*(f12[i][j][k][l]*phi1[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi1[j][l])-
 1/4*(rho11[i][j][k][l]* (rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
 rho12[i][j][k][l]*(rho12[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi1[j][l]))+

     f11[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l]+ f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     f12[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi2[j][l]+ f21[i][j][k][l]*phi2[i][k]*phi2[j][l])-
  1/4*(rho11[i][j][k][l]* (rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     rho12[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi2[j][l]))
   )
 );
  return r;
}

double sigmaF22(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((pow(f11[i][j][k][l],2) + f12[i][j][k][l]*f12[i][j][k][l] + f21[i][j][k][l]*f21[i][j][k][l]+pow(f22[i][j][k][l],2)- 1/4*(pow(rho11[i][j][k][l],2) + rho12[i][j][k][l]*rho12[i][j][k][l]+ rho21[i][j][k][l]*rho21[i][j][k][l]+pow(rho22[i][j][k][l],2)))*phi2[i][k]*phi2[j][l] +
           2*(f22[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l]+
                           f12[i][j][k][l]*phi1[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi1[j][l])-
              1/4*(rho22[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l]+
                              rho12[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi1[j][l]) )+

           f22[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
           f12[i][j][k][l]*(f12[i][j][k][l]*phi2[i][k]*phi2[j][l] + f11[i][j][k][l]*phi2[i][k]*phi1[j][l])-
      1/4*(rho22[i][j][k][l]* (rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      rho12[i][j][k][l]*(rho12[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi2[i][k]*phi1[j][l]))+

      f22[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      f21[i][j][k][l]*(f21[i][j][k][l]*phi2[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi2[j][l])-
 1/4*(rho22[i][j][k][l]* (rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
 rho21[i][j][k][l]*(rho21[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi2[j][l]))+

     f22[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l]+ f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     f21[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi1[j][l]+ f12[i][j][k][l]*phi1[i][k]*phi1[j][l])-
  1/4*(rho22[i][j][k][l]* (rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     rho21[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi1[j][l]))
   )
 );
  return r;
}

double sigmaF21(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((pow(f11[i][j][k][l],2) + f12[i][j][k][l]*f12[i][j][k][l] + f21[i][j][k][l]*f21[i][j][k][l]+pow(f22[i][j][k][l],2)- 1/4*(pow(rho11[i][j][k][l],2) + rho12[i][j][k][l]*rho12[i][j][k][l]+ rho21[i][j][k][l]*rho21[i][j][k][l]+pow(rho22[i][j][k][l],2)))*phi2[i][k]*phi1[j][l] +
           2*(f21[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                           f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l])-
              1/4*(rho21[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                              rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l]) )+

           f11[i][j][k][l]*(f11[i][j][k][l]*phi2[i][k]*phi1[j][l] + f12[i][j][k][l]*phi2[i][k]*phi2[j][l])+
           f21[i][j][k][l]*(f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l])-
      1/4*(rho11[i][j][k][l]* (rho11[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi2[i][k]*phi2[j][l])+
      rho21[i][j][k][l]*(rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l]))+

      f21[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      f22[i][j][k][l]*(f12[i][j][k][l]*phi1[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi1[j][l])-
 1/4*(rho21[i][j][k][l]* (rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
 rho22[i][j][k][l]*(rho12[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi1[j][l]))+

     f21[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l]+ f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     f22[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi2[j][l]+ f21[i][j][k][l]*phi2[i][k]*phi2[j][l])-
  1/4*(rho21[i][j][k][l]* (rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     rho22[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi2[j][l]))
   )
 );
  return r;
}

double sigmaF12(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((pow(f11[i][j][k][l],2) + f12[i][j][k][l]*f12[i][j][k][l] + f21[i][j][k][l]*f21[i][j][k][l]+pow(f22[i][j][k][l],2)- 1/4*(pow(rho11[i][j][k][l],2) + rho12[i][j][k][l]*rho12[i][j][k][l]+ rho21[i][j][k][l]*rho21[i][j][k][l]+pow(rho22[i][j][k][l],2)))*phi1[i][k]*phi2[j][l] +
           2*(f12[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                           f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l])-
              1/4*(rho12[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                              rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l]) )+

           f22[i][j][k][l]*(f22[i][j][k][l]*phi1[i][k]*phi2[j][l] + f21[i][j][k][l]*phi1[i][k]*phi1[j][l])+
           f12[i][j][k][l]*(f12[i][j][k][l]*phi1[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi1[j][l])-
      1/4*(rho22[i][j][k][l]* (rho22[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi1[i][k]*phi1[j][l])+
      rho12[i][j][k][l]*(rho12[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi1[j][l]))+

      f12[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      f11[i][j][k][l]*(f21[i][j][k][l]*phi2[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi2[j][l])-
 1/4*(rho12[i][j][k][l]* (rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
 rho11[i][j][k][l]*(rho21[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi2[j][l]))+

     f12[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l]+ f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     f11[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi1[j][l]+ f12[i][j][k][l]*phi1[i][k]*phi1[j][l])-
  1/4*(rho12[i][j][k][l]* (rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     rho11[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi1[j][l]))
   )
 );
  return r;
}

double sigmarho11(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((2*f11[i][j][k][l]*rho11[i][j][k][l] + 2*f12[i][j][k][l]*rho12[i][j][k][l]+2*f21[i][j][k][l]*rho21[i][j][k][l]+2*f22[i][j][k][l]*rho22[i][j][k][l])*phi1[i][k]*phi1[j][l] +

              2*(f11[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                           rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l])+
              1*(rho11[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                              f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l]) )+

           f11[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
           f21[i][j][k][l]*(rho21[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      1*(rho11[i][j][k][l]* (f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      rho21[i][j][k][l]*(f21[i][j][k][l]*phi1[i][k]*phi1[j][l] + f22[i][j][k][l]*phi1[i][k]*phi2[j][l]))+

      f11[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      f12[i][j][k][l]*(rho12[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi1[j][l])+
 1*(rho11[i][j][k][l]* (f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
 rho12[i][j][k][l]*(f12[i][j][k][l]*phi1[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi1[j][l]))+

     f11[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l]+ rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     f12[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi2[j][l]+ rho21[i][j][k][l]*phi2[i][k]*phi2[j][l])+
  1*(rho11[i][j][k][l]* (f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     rho12[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi2[j][l] + f21[i][j][k][l]*phi2[i][k]*phi2[j][l]))
   )
 );
  return r;
}

double sigmarho22(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((2*f11[i][j][k][l]*rho11[i][j][k][l] + 2*f12[i][j][k][l]*rho12[i][j][k][l]+2*f21[i][j][k][l]*rho21[i][j][k][l]+2*f22[i][j][k][l]*rho22[i][j][k][l])*phi2[i][k]*phi2[j][l] +

           2*(f22[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l]+
                           rho12[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi1[j][l])+
              1*(rho22[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l]+
                              f12[i][j][k][l]*phi1[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi1[j][l]) )+

           f22[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
           f12[i][j][k][l]*(rho12[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      1*(rho22[i][j][k][l]* (f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      rho12[i][j][k][l]*(f12[i][j][k][l]*phi2[i][k]*phi2[j][l] + f11[i][j][k][l]*phi2[i][k]*phi1[j][l]))+

      f22[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      f21[i][j][k][l]*(rho21[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi2[j][l])+
 1*(rho22[i][j][k][l]* (f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
 rho21[i][j][k][l]*(f21[i][j][k][l]*phi2[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi2[j][l]))+

     f22[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l]+ rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     f21[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi1[j][l]+ rho12[i][j][k][l]*phi1[i][k]*phi1[j][l])+
  1*(rho22[i][j][k][l]* (f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     rho21[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi1[j][l]))
   )
 );
  return r;
}

double sigmarho21(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((2*f11[i][j][k][l]*rho11[i][j][k][l] + 2*f12[i][j][k][l]*rho12[i][j][k][l]+2*f21[i][j][k][l]*rho21[i][j][k][l]+2*f22[i][j][k][l]*rho22[i][j][k][l])*phi2[i][k]*phi1[j][l] +

           2*(f21[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                           rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l])+
              1*(rho21[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                              f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l]) )+

           f11[i][j][k][l]*(rho11[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi2[i][k]*phi2[j][l])+
           f21[i][j][k][l]*(rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l])+
      1*(rho11[i][j][k][l]* (f11[i][j][k][l]*phi2[i][k]*phi1[j][l] + f12[i][j][k][l]*phi2[i][k]*phi2[j][l])+
      rho21[i][j][k][l]*(f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l]))+

      f21[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
      f22[i][j][k][l]*(rho12[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi1[j][l])+
 1*(rho21[i][j][k][l]* (f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
 rho22[i][j][k][l]*(f12[i][j][k][l]*phi1[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi1[j][l]))+

     f21[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l]+ rho21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     f22[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi2[j][l]+ rho21[i][j][k][l]*phi2[i][k]*phi2[j][l])+
  1*(rho21[i][j][k][l]* (f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f21[i][j][k][l]*phi2[i][k]*phi1[j][l])+
     rho22[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi2[j][l] + f21[i][j][k][l]*phi2[i][k]*phi2[j][l]))
   )
 );
  return r;
}

double sigmarho12(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r= -2* pow((u/4),2)*
          ((2*f11[i][j][k][l]*rho11[i][j][k][l] + 2*f12[i][j][k][l]*rho12[i][j][k][l]+2*f21[i][j][k][l]*rho21[i][j][k][l]+2*f22[i][j][k][l]*rho22[i][j][k][l])*phi1[i][k]*phi2[j][l] +

           2*(f12[i][j][k][l]*(rho11[i][j][k][l]*phi1[i][k]*phi1[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                           rho21[i][j][k][l]*phi2[i][k]*phi1[j][l] + rho22[i][j][k][l]*phi2[i][k]*phi2[j][l])+
              1*(rho12[i][j][k][l]*(f11[i][j][k][l]*phi1[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l]+
                              f21[i][j][k][l]*phi2[i][k]*phi1[j][l] + f22[i][j][k][l]*phi2[i][k]*phi2[j][l]) )+

           f22[i][j][k][l]*(rho22[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho21[i][j][k][l]*phi1[i][k]*phi1[j][l])+
           f12[i][j][k][l]*(rho12[i][j][k][l]*phi1[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi1[j][l])+
      1*(rho22[i][j][k][l]* (f22[i][j][k][l]*phi1[i][k]*phi2[j][l] + f21[i][j][k][l]*phi1[i][k]*phi1[j][l])+
      rho12[i][j][k][l]*(f12[i][j][k][l]*phi1[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi1[j][l]))+

      f12[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
      f11[i][j][k][l]*(rho21[i][j][k][l]*phi2[i][k]*phi2[j][l] + rho11[i][j][k][l]*phi1[i][k]*phi2[j][l])+
 1*(rho12[i][j][k][l]* (f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
 rho11[i][j][k][l]*(f21[i][j][k][l]*phi2[i][k]*phi2[j][l] + f11[i][j][k][l]*phi1[i][k]*phi2[j][l]))+

     f12[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi2[j][l]+ rho12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     f11[i][j][k][l]*(rho22[i][j][k][l]*phi2[i][k]*phi1[j][l]+ rho12[i][j][k][l]*phi1[i][k]*phi1[j][l])+
  1*(rho12[i][j][k][l]* (f22[i][j][k][l]*phi2[i][k]*phi2[j][l] + f12[i][j][k][l]*phi1[i][k]*phi2[j][l])+
     rho11[i][j][k][l]*(f22[i][j][k][l]*phi2[i][k]*phi1[j][l] + f12[i][j][k][l]*phi1[i][k]*phi1[j][l]))
   )
 );
  return r;
}

/* Define integrals ijkl i is x0*/
double i11a(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  for (a=0;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f11[a][j][b][l] + sigmarho12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f21[a][j][b][l] ) +r;
    }
  }
  for (a=0;a<j;a++){
    for (b=0; b<Nx;b++){
      r= -at * ax * (sigmaF11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho11[a][j][b][l] + sigmaF12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho21[a][j][b][l] ) +r;
    }
  }
  return r;
}

double i12a(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  for (a=0;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f12[a][j][b][l] + sigmarho12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f22[a][j][b][l] ) +r;
    }
  }
  for (a=0;a<j;a++){
    for (b=0; b<Nx;b++){
      r= -at * ax * (sigmaF11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho12[a][j][b][l] + sigmaF12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho22[a][j][b][l] ) +r;
    }
  }
  return r;
}

double i21a(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  for (a=0;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho22(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f21[a][j][b][l] + sigmarho21(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f11[a][j][b][l] ) +r;
    }
  }
  for (a=0;a<j;a++){
    for (b=0; b<Nx;b++){
      r= -at * ax * (sigmaF22(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho21[a][j][b][l] + sigmaF21(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho11[a][j][b][l] ) +r;
    }
  }
  return r;
}

double i22a(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  for (a=0;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho22(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f22[a][j][b][l] + sigmarho21(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*f12[a][j][b][l] ) +r;
    }
  }
  for (a=0;a<j;a++){
    for (b=0; b<Nx;b++){
      r= -at * ax * (sigmaF22(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho22[a][j][b][l] + sigmaF21(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho12[a][j][b][l] ) +r;
    }
  }
  return r;
}

double i11rhoa(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  if (i > j ){
  for (a=j;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho11[a][j][b][l] + sigmarho12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho21[a][j][b][l] ) +r;
    }
  }
}
if (j>i) {
  for (a=i;a<j;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho11[a][j][b][l] + sigmarho12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho21[a][j][b][l] ) +r;
    }
  }
}
  return r;
}

double i12rhoa(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  if (i >j){
  for (a=j;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho12[a][j][b][l] + sigmarho12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho22[a][j][b][l] ) +r;
    }
  }
}
if (j>i){
  for (a=i;a<j;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho11(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho12[a][j][b][l] + sigmarho12(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho22[a][j][b][l] ) +r;
    }
  }

}
  return r;
}

double i21rhoa(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  r=0;
  int a,b;
  if (i> j){
  for (a=j;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho21(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho11[a][j][b][l] + sigmarho22(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho21[a][j][b][l] ) +r;
    }
  }
}
if (j >i){
  for (a=i;a<j;a++){
    for (b=0; b<Nx;b++){
      r=at * ax * (sigmarho21(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho11[a][j][b][l] + sigmarho22(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho21[a][j][b][l] ) +r;
    }
  }
}
  return r;
}

double i22rhoa(double at, double ax,int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i, int j,int k,int l, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  if (i>j){
  for (a=j;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho21(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho12[a][j][b][l] + sigmarho22(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho22[a][j][b][l] ) +r;
    }
  }
}
if (j>i){
  for (a=i;a<j;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * (sigmarho21(Nt, Nx, u, f11, f12, f21, f22,
        rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho12[a][j][b][l] + sigmarho22(Nt, Nx, u, f11, f12, f21, f22,
          rho11,rho12,rho21,rho22, i, a,k, b,  phi1, phi2)*rho22[a][j][b][l] ) +r;
    }
  }
}
  return r;
}

double DGamma1(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i,int k, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  for (a=0;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * 2*pow(u/4,2)*(
        (pow(f11[i][a][k][b],2) + f12[i][a][k][b]*f12[i][a][k][b] + f21[i][a][k][b]*f21[i][a][k][b]+pow(f22[i][a][k][b],2)- 1/4*(pow(rho11[i][a][k][b],2) + rho12[i][a][k][b]*rho12[i][a][k][b]+ rho21[i][a][k][b]*rho21[i][a][k][b]+ pow(rho22[i][a][k][b],2)))*(rho11[i][a][k][b]*phi1[a][b] + rho12[i][a][k][b]*phi2[a][b])+

           2*(f11[i][a][k][b]*rho11[i][a][k][b]+ f12[i][a][k][b]*rho12[i][a][k][b]+ f21[i][a][k][b]*rho21[i][a][k][b]+f22[i][a][k][b]*rho22[i][a][k][b])*(f11[i][a][k][b]*phi1[a][b] + f12[i][a][k][b]*phi2[a][b])
      )
      + at * ax * 4*pow(u/4,2)*(
             f11[i][a][k][b]*( rho11[i][a][k][b]*(f12[i][a][k][b]*phi2[a][b] + f11[i][a][k][b]*phi1[a][b]) + rho21[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b]))+
             f12[i][a][k][b]*( rho12[i][a][k][b]*(f11[i][a][k][b]*phi1[a][b]+ f12[i][a][k][b]*phi2[a][b]) + rho22[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b]))+

             f11[i][a][k][b]*( f11[i][a][k][b]*(rho12[i][a][k][b]*phi2[a][b] + rho11[i][a][k][b]*phi1[a][b]) + f21[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b]+rho22[i][a][k][b]*phi2[a][b]))+
             f12[i][a][k][b]*( f12[i][a][k][b]*(rho11[i][a][k][b]*phi1[a][b] + rho12[i][a][k][b]*phi2[a][b])+ f22[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b] +rho22[i][a][k][b]*phi2[a][b]))+
             rho11[i][a][k][b]*(f11[i][a][k][b]*(f11[i][a][k][b]*phi1[a][b]+f12[i][a][k][b]*phi2[a][b]) + f21[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b]+f22[i][a][k][b]*phi2[a][b]))+
             rho12[i][a][k][b]*(f12[i][a][k][b]*(f11[i][a][k][b]*phi1[a][b]+f12[i][a][k][b]*phi2[a][b])+ f22[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b]))+
             -1/4*rho11[i][a][k][b]*(rho11[i][a][k][b]*(rho11[i][a][k][b]*phi1[a][b]+ rho12[i][a][k][b]*phi2[a][b])+ rho21[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b] +rho22[i][a][k][b]*phi2[a][b]))+
             -1/4*rho12[i][a][k][b]*(rho12[i][a][k][b]*(rho11[i][a][k][b]*phi1[a][b]+ rho12[i][a][k][b]*phi2[a][b])+ rho22[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b] + rho22[i][a][k][b]*phi2[a][b]))
       )
      +r;
    }
  }
  return r;
}

double DGamma2(int Nt,int Nx, double u, double ****f11,double ****f12,double ****f21,double ****f22,
  double ****rho11,double ****rho12,double ****rho21,double ****rho22, int i,int k, double phi1[Nt][Nx],double phi2[Nt][Nx]){
  double r;
  int a,b;
  r=0;
  for (a=0;a<i;a++){
    for (b=0; b<Nx;b++){
      r= at * ax * 2*pow(u/4,2)*(
        (pow(f11[i][a][k][b],2) + f12[i][a][k][b]*f12[i][a][k][b]+  f21[i][a][k][b]*f21[i][a][k][b] +pow(f22[i][a][k][b],2)- 1/4*(pow(rho11[i][a][k][b],2) + rho12[i][a][k][b]*rho12[i][a][k][b]+ rho21[i][a][k][b]*rho21[i][a][k][b] +pow(rho22[i][a][k][b],2)))*(rho21[i][a][k][b]*phi1[a][b] + rho22[i][a][k][b]*phi2[a][b])+

           2*(f11[i][a][k][b]*rho11[i][a][k][b]+ f12[i][a][k][b]*rho12[i][a][k][b]+ f21[i][a][k][b]*rho21[i][a][k][b]+f22[i][a][k][b]*rho22[i][a][k][b])*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b])
      )
      + at * ax * 4*pow(u/4,2)*(
             f21[i][a][k][b]*( rho11[i][a][k][b]*(f12[i][a][k][b]*phi2[a][b] + f11[i][a][k][b]*phi1[a][b]) + rho21[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b]))+
             f22[i][a][k][b]*( rho12[i][a][k][b]*(f11[i][a][k][b]*phi1[a][b]+ f12[i][a][k][b]*phi2[a][b]) + rho22[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b]))+

             f21[i][a][k][b]*( f11[i][a][k][b]*(rho12[i][a][k][b]*phi2[a][b] + rho11[i][a][k][b]*phi1[a][b]) + f21[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b]+rho22[i][a][k][b]*phi2[a][b]))+
             f22[i][a][k][b]*( f12[i][a][k][b]*(rho11[i][a][k][b]*phi1[a][b] + rho12[i][a][k][b]*phi2[a][b])+ f22[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b] +rho22[i][a][k][b]*phi2[a][b]))+
             rho21[i][a][k][b]*(f11[i][a][k][b]*(f11[i][a][k][b]*phi1[a][b]+f12[i][a][k][b]*phi2[a][b]) + f21[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b]+f22[i][a][k][b]*phi2[a][b]))+
             rho22[i][a][k][b]*(f12[i][a][k][b]*(f11[i][a][k][b]*phi1[a][b]+f12[i][a][k][b]*phi2[a][b])+ f22[i][a][k][b]*(f21[i][a][k][b]*phi1[a][b] + f22[i][a][k][b]*phi2[a][b]))+
             -1/4*rho21[i][a][k][b]*(rho11[i][a][k][b]*(rho11[i][a][k][b]*phi1[a][b]+ rho12[i][a][k][b]*phi2[a][b])+ rho21[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b] +rho22[i][a][k][b]*phi2[a][b]))+
             -1/4*rho22[i][a][k][b]*(rho12[i][a][k][b]*(rho11[i][a][k][b]*phi1[a][b]+ rho12[i][a][k][b]*phi2[a][b])+ rho22[i][a][k][b]*(rho21[i][a][k][b]*phi1[a][b] + rho22[i][a][k][b]*phi2[a][b]))
       )
      +r;
    }
  }
  return r;
}

/* Define central right and left derivative*/
double D(int Nt,int Nx,int i, int j,int k, int l,double ****a,double b){
      double z;
      if (k==0){
        z= 1/pow(b,2)*(a[i][j][1][l]+a[i][j][Nx-1][l]-2*a[i][j][0][l]);
      }
      else if (k== Nx-1){
          z= 1/pow(b,2)*(a[i][j][0][l]+a[i][j][Nx-2][l]-2*a[i][j][Nx-1][l]);
       }
       else {
         z= 1/pow(b,2)*(a[i][j][k+1][l]+a[i][j][k-1][l]-2*a[i][j][k][l]);
       }
      return z;
    }

double D0(int Nt,int Nx,int i, int k, double a[Nt][Nx],double b){
      double l;
      if (k==0){
        l= 1/pow(b,2)*(a[i][k+1]+a[i][Nx-1]-2*a[i][k]);
      }
      else if (k== Nx-1){
          l= 1/pow(b,2)*(a[i][0]+a[i][k-1]-2*a[i][k]);
        }
        else {
          l= 1/pow(b,2)*(a[i][k+1]+a[i][k-1]-2*a[i][k]);
        }
        return l;
    }

double np0(int Nx, double ax, double sigma, int i,int j,int s){
  double l;
  l= 2*cos( 2*PI/(Nx)* s *(i-j))* (1*(sqrt(2*PI/pow(sigma,2)))/pow((ax*Nx),2)*exp(-1/(2*pow(sigma,2))*pow((2*PI*s/(Nx*ax)),2))+.5/(Nx*ax));
return l;
}

double f00(int Nx, double ax, double sigma, int i,int j){
      double l;
      int n;
      double s;
      l=0;
      for (n=1;n<=(21-1)/2;n++) {
        l=l + np0(Nx,ax,sigma,i,j,n)/(sqrt(pow((2*PI*n/(Nx*ax)),2)+1));
      }
      l=l+ np0(Nx,ax,sigma,i,j,0)/2;
      /*return l ;*/
      return l;
    }

double ren1(int Nx,double ax, double sigma){
  double l;
  int n;
  l= f00(Nx, ax,sigma, 0,0);
  /*return -l - (1*(sqrt(2*PI/pow(sigma,2)))*9/2/(Nx*ax)-0/(ax*Nx));*/
  /*return l;*/
  return -l;
}

double f00_22(int Nx, double ax, double sigma, int i,int j){
          double l;
          int n;
          double s;
          l=0;
          for (n=1;n<=(21-1)/2;n++) {
            l= l + 2*cos( 2*PI/(Nx*ax)* n *(i-j))* (.5/(Nx*ax));
          }
          return l + .5/(Nx*ax) ;
        }

/*        double f00_22(int Nx, double ax, double sigma, int i,int j){
              double l;
              if (i==j){l=.5;}
              else {l=0;}
              return l;
            } */

            double ren2(int Nx,double ax, double sigma){
              double l;
              int n;
              l= f00(Nx,ax,sigma,0,0);

              /*return -l;*/
              return -l;
            }

            double f00_11(int Nx, double ax, double sigma, int i,int j){
                  double l;
                  int n;
                  double s;
                  l=0;
                  for (n=1;n<=(21-1)/2;n++) {

                    l=l + np0(Nx,ax,sigma,i,j,n)*(sqrt(pow((2*PI*n/(Nx*ax)),2)+1));
                  }
                  l=l+ np0(Nx,ax,sigma,i,j,0)/2;
                  /*return l ;*/
                  return l;
                }

int main () {
  /*Allocating memory in the heap for the arrays, this will be painful for the brain*/
  int r1=Nt,r2=Nt,r3=Nx,r4=Nx;
  double a;
  double sigma;
  sigma =2*PI/(Nx*ax);
  int i1,j1,k1;
  double *allElements11 = malloc(r1*r2*r3*r4*sizeof(double));
  double *allElements12 = malloc(r1*r2*r3*r4*sizeof(double));
  double *allElements21 = malloc(r1*r2*r3*r4*sizeof(double));
  double *allElements22 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rallElements11 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rallElements12 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rallElements21 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rallElements22 = malloc(r1*r2*r3*r4*sizeof(double));

  double *rrallElements11 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrallElements12 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrallElements21 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrallElements22 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrallElements11 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrallElements12 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrallElements21 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrallElements22 = malloc(r1*r2*r3*r4*sizeof(double));

  double *rrrrallElements11 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrallElements12 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrallElements21 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrallElements22 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrrallElements11 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrrallElements12 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrrallElements21 = malloc(r1*r2*r3*r4*sizeof(double));
  double *rrrrrallElements22 = malloc(r1*r2*r3*r4*sizeof(double));

  double ****f11= malloc(r1 * sizeof(double ***));
  double ****f12= malloc(r1 * sizeof(double ***));
  double ****f21= malloc(r1 * sizeof(double ***));
  double ****f22= malloc(r1 * sizeof(double ***));
  double ****rho11= malloc(r1 * sizeof(double ***));
  double ****rho12= malloc(r1 * sizeof(double ***));
  double ****rho21= malloc(r1 * sizeof(double ***));
  double ****rho22= malloc(r1 * sizeof(double ***));

  double ****i11= malloc(r1 * sizeof(double ***));
  double ****i12= malloc(r1 * sizeof(double ***));
  double ****i21= malloc(r1 * sizeof(double ***));
  double ****i22= malloc(r1 * sizeof(double ***));
  double ****i11rho= malloc(r1 * sizeof(double ***));
  double ****i12rho= malloc(r1 * sizeof(double ***));
  double ****i21rho= malloc(r1 * sizeof(double ***));
  double ****i22rho= malloc(r1 * sizeof(double ***));

  double ****ii11= malloc(r1 * sizeof(double ***));
  double ****ii12= malloc(r1 * sizeof(double ***));
  double ****ii21= malloc(r1 * sizeof(double ***));
  double ****ii22= malloc(r1 * sizeof(double ***));
  double ****ii11rho= malloc(r1 * sizeof(double ***));
  double ****ii12rho= malloc(r1 * sizeof(double ***));
  double ****ii21rho= malloc(r1 * sizeof(double ***));
  double ****ii22rho= malloc(r1 * sizeof(double ***));

  for(i1=0; i1 < r1;i1++){
    f11[i1] = malloc(r2*sizeof(double **));
    f12[i1] = malloc(r2*sizeof(double **));
    f21[i1] = malloc(r2*sizeof(double **));
    f22[i1] = malloc(r2*sizeof(double **));
    rho11[i1] = malloc(r2*sizeof(double **));
    rho12[i1] = malloc(r2*sizeof(double **));
    rho21[i1] = malloc(r2*sizeof(double **));
    rho22[i1] = malloc(r2*sizeof(double **));

    i11[i1] = malloc(r2*sizeof(double **));
    i12[i1] = malloc(r2*sizeof(double **));
    i21[i1] = malloc(r2*sizeof(double **));
    i22[i1] = malloc(r2*sizeof(double **));
    i11rho[i1] = malloc(r2*sizeof(double **));
    i12rho[i1] = malloc(r2*sizeof(double **));
    i21rho[i1] = malloc(r2*sizeof(double **));
    i22rho[i1] = malloc(r2*sizeof(double **));

    ii11[i1] = malloc(r2*sizeof(double **));
    ii12[i1] = malloc(r2*sizeof(double **));
    ii21[i1] = malloc(r2*sizeof(double **));
    ii22[i1] = malloc(r2*sizeof(double **));
    ii11rho[i1] = malloc(r2*sizeof(double **));
    ii12rho[i1] = malloc(r2*sizeof(double **));
    ii21rho[i1] = malloc(r2*sizeof(double **));
    ii22rho[i1] = malloc(r2*sizeof(double **));

    for(j1=0;j1<r2;j1++){
      f11[i1][j1] = malloc(r3*sizeof(double *));
      f12[i1][j1] = malloc(r3*sizeof(double *));
      f21[i1][j1] = malloc(r3*sizeof(double *));
      f22[i1][j1] = malloc(r3*sizeof(double *));
      rho11[i1][j1] = malloc(r3*sizeof(double *));
      rho12[i1][j1] = malloc(r3*sizeof(double *));
      rho21[i1][j1] = malloc(r3*sizeof(double *));
      rho22[i1][j1] = malloc(r3*sizeof(double *));
      i11[i1][j1] = malloc(r3*sizeof(double *));
      i12[i1][j1] = malloc(r3*sizeof(double *));
      i21[i1][j1] = malloc(r3*sizeof(double *));
      i22[i1][j1] = malloc(r3*sizeof(double *));
      i11rho[i1][j1] = malloc(r3*sizeof(double *));
      i12rho[i1][j1] = malloc(r3*sizeof(double *));
      i21rho[i1][j1] = malloc(r3*sizeof(double *));
      i22rho[i1][j1] = malloc(r3*sizeof(double *));
      ii11[i1][j1] = malloc(r3*sizeof(double *));
      ii12[i1][j1] = malloc(r3*sizeof(double *));
      ii21[i1][j1] = malloc(r3*sizeof(double *));
      ii22[i1][j1] = malloc(r3*sizeof(double *));
      ii11rho[i1][j1] = malloc(r3*sizeof(double *));
      ii12rho[i1][j1] = malloc(r3*sizeof(double *));
      ii21rho[i1][j1] = malloc(r3*sizeof(double *));
      ii22rho[i1][j1] = malloc(r3*sizeof(double *));

      for(k1=0;k1<r3;k1++){
        f11[i1][j1][k1]=allElements11 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        f12[i1][j1][k1]=allElements12 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        f21[i1][j1][k1]=allElements21 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        f22[i1][j1][k1]=allElements22 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        rho11[i1][j1][k1]=rallElements11 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        rho12[i1][j1][k1]=rallElements12 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        rho21[i1][j1][k1]=rallElements21 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        rho22[i1][j1][k1]=rallElements22 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);

        i11[i1][j1][k1]=rrallElements11 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i12[i1][j1][k1]=rrallElements12 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i21[i1][j1][k1]=rrallElements21 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i22[i1][j1][k1]=rrallElements22 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i11rho[i1][j1][k1]=rrrallElements11 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i12rho[i1][j1][k1]=rrrallElements12 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i21rho[i1][j1][k1]=rrrallElements21 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        i22rho[i1][j1][k1]=rrrallElements22 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);

        ii11[i1][j1][k1]=rrrrallElements11 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii12[i1][j1][k1]=rrrrallElements12 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii21[i1][j1][k1]=rrrrallElements21 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii22[i1][j1][k1]=rrrrallElements22 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii11rho[i1][j1][k1]=rrrrrallElements11 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii12rho[i1][j1][k1]=rrrrrallElements12 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii21rho[i1][j1][k1]=rrrrrallElements21 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
        ii22rho[i1][j1][k1]=rrrrrallElements22 + (i1 * r2*r3*r4) + (j1*r3*r4)+ (k1*r4);
      }
    }
  }


  /* Commenting this out (this would be for the variables in the stack)
   double f11[Nt][Nt][Nx][Nx];
    double f12[Nt][Nt][Nx][Nx];
    double f21[Nt][Nt][Nx][Nx];
    double f22[Nt][Nt][Nx][Nx];
    */
    double phi1[Nt][Nx];
    double phi2[Nt][Nx];

    int i,j,l,m,n,k,o,p,q,s,t,n0,m0,n1,z;
    /* Initial conditions of the cycle*/
    for (s=0; s< Nx;s++){
        phi1[0][s]=0.5;
        phi2[0][s]=0;
        phi1[1][s]= phi1[0][s];
        phi2[1][s]= phi1[0][s]* at;
        }
    for (t=0; t<Nx;t++){
       for (z=0; z<Nx;z++){
          f11[0][0][t][z]=f00(Nx,ax,sigma,t,z);
         f12[0][0][t][z]=0;
         f21[0][0][t][z]=0;
         f22[0][0][t][z]=f00(Nx,ax,sigma,t,z);
         rho11[0][0][t][z]=0;
         rho22[0][0][t][z]=0;
         rho12[0][0][t][z]=0;
         rho21[0][0][t][z]=0;
       }}
       for (t=0; t<Nx;t++){
          for (z=0; z<Nx;z++){
         f11[1][0][t][z]=f00(Nx,ax,sigma,t,z) + 0.5*pow(at,2)*(D(Nt,Nx,0,0,t,z,f11,ax)-f11[0][0][t][z]);
         f12[1][0][t][z]=0;
         f21[1][0][t][z]=0;
         f22[1][0][t][z]=f00(Nx,ax,sigma,t,z) + 0.5*pow(at,2)*(D(Nt,Nx,0,0,t,z,f22,ax)-f22[0][0][t][z]);
          f11[0][1][z][t]=f00(Nx,ax,sigma,z,t) + 0.5*pow(at,2)*(D(Nt,Nx,0,0,z,t,f11,ax)-f11[0][0][z][t]);
         f12[0][1][t][z]=0;
         f21[0][1][t][z]=0;
         f22[0][1][z][t]=f00(Nx,ax,sigma,z,t) + 0.5*pow(at,2)*(D(Nt,Nx,0,0,z,t,f22,ax)-f22[0][0][z][t]);
         rho11[1][0][t][z]=0;
         rho22[1][0][t][z]=0;
         rho12[1][0][t][z]=0;
         rho21[1][0][t][z]=0;
         rho11[0][1][t][z]=0;
         rho22[0][1][t][z]=0;
         rho12[0][1][t][z]=0;
         rho21[0][1][t][z]=0;
       }}
       for (t=0; t<Nx;t++){
          for (z=0; z<Nx;z++){
         f11[1][1][t][z]=pow(at,2)*f00_11(Nx,ax,sigma,t,z)+f11[1][0][t][z] + .5*pow(at,2)*(D(Nt,Nx,0,1,t,z,f11,ax)- f11[0][0][t][z]);
         f12[1][1][t][z]=0;
         f21[1][1][t][z]=0;
         f22[1][1][t][z]=pow(at,2)*f00_11(Nx,ax,sigma,t,z)+f22[1][0][t][z] + .5*pow(at,2)*(D(Nt,Nx,0,1,t,z,f22,ax)- f22[0][0][t][z]);
         if(t==z){rho11[1][0][t][z]=-at/ax;
                  rho22[1][0][t][z]=-at/ax;
                  rho11[0][1][z][t]=at/ax;
                  rho22[0][1][z][t]=at/ax;}
        else if (t!=z){rho11[1][0][t][z]=0;
                 rho22[1][0][t][z]=0;
                 rho11[0][1][t][z]=0;
                 rho22[0][1][t][z]=0;}
              }
        }

        for (t=0; t<Nx;t++){
           for (z=0; z<Nx;z++){
             f11[1][1][t][z]= .5*(f11[1][1][t][z]+ f11[1][1][z][t]);
             f22[1][1][t][z]= .5*(f22[1][1][t][z]+ f22[1][1][z][t]);
           }}

  for (j=2;j<Nt;j++){


    /*Equation of motion for phi*/

      for (n1=0;n1<Nx;n1++){
        phi1[j][n1]=2*phi1[j-1][n1]- phi1[j-2][n1]- pow(at,2)*(-D0(Nt,Nx,j-1,n1,phi1,ax) + phi1[j-1][n1]
                                    + u/4*(pow(phi1[j-1][n1],2)+pow(phi2[j-1][n1],2))*phi1[j-1][n1]
                                    + u/4*phi1[j-1][n1]*(3*(f11[j-1][j-1][n1][n1]+ ren1(Nx,ax,sigma))+ f12[j-1][j-1][n1][n1] + f21[j-1][j-1][n1][n1] + (f22[j-1][j-1][n1][n1]+ren2(Nx,ax,sigma))) +
                                    DGamma1(Nt,Nx,u,f11,f12,f21,f22,
                                      rho11,rho12,rho21,rho22, j-1, n1, phi1, phi2)
                                  );

                                  phi2[j][n1]=2*phi2[j-1][n1]-phi2[j-2][n1]- pow(at,2)*( -D0(Nt,Nx,j-1,n1,phi2,ax) + phi2[j-1][n1]
                                                              + u/4*(pow(phi1[j-1][n1],2)+pow(phi2[j-1][n1],2))*phi2[j-1][n1]
                                                              + u/4*phi2[j-1][n1]*(3*(f22[j-1][j-1][n1][n1]+ren2(Nx,ax,sigma))+ f21[j-1][j-1][n1][n1] + f12[j-1][j-1][n1][n1] + (f11[j-1][j-1][n1][n1]+ ren1(Nx,ax,sigma))) +
                                                              DGamma2(Nt,Nx,u,f11,f12,f21,f22,
                                                                rho11,rho12,rho21,rho22, j-1, n1, phi1, phi2)
                                                            );

  }

              for (k=0; k<j; k++){
              #pragma omp parallel shared(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j,k, n,m,phi1,phi2,i11,i12,i21,i22,i11rho,i12rho,i21rho,i22rho,ii11,ii12,ii21,ii22)
                {
                  #pragma omp for collapse(2)
                  for (n=0; n< Nx;n++){
                   for (m=0; m<Nx;m++){
                     ii11[j-1][k][n][m]= i11a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                     f11[j][k][n][m]= 2*f11[j-1][k][n][m]-f11[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f21[j-1][k][n][m]
                                       - D(Nt,Nx,j-1,k,n,m,f11,ax) + f11[j-1][k][n][m]
                                       + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*f11[j-1][k][n][m]
                                       + u/4*(3* (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f11[j-1][k][n][m]+ (f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f11[j-1][k][n][m]+
                                               2*f12[j-1][j-1][n][n]*f21[j-1][k][n][m])+

                                         ii11[j-1][k][n][m]
                                   );
                                   ii11rho[j-1][k][n][m]= i11rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                                   rho11[j][k][n][m]= 2*rho11[j-1][k][n][m]-rho11[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho21[j-1][k][n][m]
                                                     - D(Nt,Nx,j-1,k,n,m,rho11,ax) + rho11[j-1][k][n][m]
                                                     + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*rho11[j-1][k][n][m]
                                                     + u/4*(3* (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho11[j-1][k][n][m]+ (f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho11[j-1][k][n][m]+
                                                           2*  f12[j-1][j-1][n][n]*rho21[j-1][k][n][m]) +

                                                             ii11rho[j-1][k][n][m]

                                                     );

                                                     ii12[j-1][k][n][m]= i12a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                                                     f12[j][k][n][m]= 2*f12[j-1][k][n][m]-f12[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f22[j-1][k][n][m]
                                                                       - D(Nt,Nx,j-1,k,n,m,f12,ax) + f12[j-1][k][n][m]
                                                                       + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*f12[j-1][k][n][m]
                                                                       + u/4*(3*(f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f12[j-1][k][n][m]+ (f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f12[j-1][k][n][m]+ 2* f12[j-1][j-1][n][n]*f22[j-1][k][n][m])+

                                                                       ii12[j-1][k][n][m]
                                                                     );

                                   ii12rho[j-1][k][n][m]= i12rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                                   rho12[j][k][n][m]= 2*rho12[j-1][k][n][m]-rho12[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho22[j-1][k][n][m]
                                                     -D(Nt,Nx,j-1,k,n,m,rho12,ax) + rho12[j-1][k][n][m]
                                                     + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*rho12[j-1][k][n][m]
                                                     + u/4*(3*(f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho12[j-1][k][n][m]+(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho12[j-1][k][n][m]+ 2* f12[j-1][j-1][n][n]*rho22[j-1][k][n][m])+

                                                     ii12rho[j-1][k][n][m]
                                                   );

                                                   ii21[j-1][k][n][m]=i21a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2) ;
                                                   f21[j][k][n][m]= 2*f21[j-1][k][n][m]- f21[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f11[j-1][k][n][m]
                                                                     - D(Nt,Nx,j-1,k,n,m,f21,ax) + f21[j-1][k][n][m]
                                                                     + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*f21[j-1][k][n][m]
                                                                     + u/4*((f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f21[j-1][k][n][m]+ 3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f21[j-1][k][n][m]+
                                                                             2*f21[j-1][j-1][n][n]*f11[j-1][k][n][m])+

                                                                             ii21[j-1][k][n][m]
                                                                           );

                                   ii21rho[j-1][k][n][m]= i21rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                                   rho21[j][k][n][m]= 2*rho21[j-1][k][n][m]-rho21[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho11[j-1][k][n][m]
                                                     -D(Nt,Nx,j-1,k,n,m,rho21,ax) + rho21[j-1][k][n][m]
                                                     + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*rho21[j-1][k][n][m]
                                                     + u/4*((f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho21[j-1][k][n][m]+ 3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho21[j-1][k][n][m]+
                                                           2*  f21[j-1][j-1][n][n]*rho11[j-1][k][n][m])+

                                                       ii21rho[j-1][k][n][m]
                                                     );

                                                     ii22[j-1][k][n][m]= i22a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                                                     f22[j][k][n][m]= 2*f22[j-1][k][n][m]-f22[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f12[j-1][k][n][m]
                                                                       - D(Nt,Nx,j-1,k,n,m,f22,ax) + f22[j-1][k][n][m]
                                                                       + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*f22[j-1][k][n][m]
                                                                       + u/4*(2*f21[j-1][j-1][n][n]*f12[j-1][k][n][m]+ (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f22[j-1][k][n][m]+
                                                                               3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f22[j-1][k][n][m])+

                                                                               ii22[j-1][k][n][m]
                                                                             );

                                                     ii22rho[j-1][k][n][m]= i22rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,k, n,m,phi1,phi2);
                                                     rho22[j][k][n][m]= 2*rho22[j-1][k][n][m]- rho22[j-2][k][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho12[j-1][k][n][m]
                                                                       - D(Nt,Nx,j-1,k,n,m,rho22,ax) + rho22[j-1][k][n][m]
                                                                       + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*rho22[j-1][k][n][m]
                                                                       + u/4*(2*f21[j-1][j-1][n][n]*rho12[j-1][k][n][m]+ (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho22[j-1][k][n][m]+
                                                                               3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho22[j-1][k][n][m])+

                                                                               ii22rho[j-1][k][n][m]
                                                                             );
                                                                           }}}}
                                                      for (k=0; k<j; k++){
                                                        for (n=0; n< Nx;n++){
                                                         for (m=0; m<Nx;m++){

                                                    f11[k][j][m][n]=f11[j][k][n][m];
                                                    f12[k][j][m][n]=f21[j][k][n][m];
                                                    f21[k][j][m][n]=f12[j][k][n][m];
                                                    f22[k][j][m][n]=f22[j][k][n][m];
                                                    rho11[k][j][m][n]=-rho11[j][k][n][m];
                                                    rho12[k][j][m][n]=-rho21[j][k][n][m];
                                                    rho21[k][j][m][n]=-rho12[j][k][n][m];
                                                    rho22[k][j][m][n]=-rho22[j][k][n][m];
               }
             }
           }


        #pragma omp parallel shared(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,i,j, n,m,phi1,phi2,i11,i12,i21,i22,i11rho,i12rho,i21rho,i22rho)
          {
            #pragma omp for collapse(2)
          for (n=0; n<Nx;n++){
            for (m=0;m<Nx;m++){
              i11[j-1][j][n][m]= i11a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
              f11[j][j][n][m]= 2*f11[j-1][j][n][m]- f11[j-2][j][n][m] - pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f21[j-1][j][n][m]
                                - D(Nt,Nx,j-1,j,n,m,f11,ax) + f11[j-1][j][n][m]
                                + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*f11[j-1][j][n][m]
                                + u/4*(3* (f11[j-1][j-1][n][n]+ ren1(Nx,ax,sigma))*f11[j-1][j][n][m]+ (f22[j-1][j-1][n][n]+ren2(Nx,ax,sigma))*f11[j-1][j][n][m]+
                                        2*f12[j-1][j-1][n][n]*f21[j-1][j][n][m]) +

                                        i11[j-1][j][n][m]

                                );
                            i11rho[j-1][j][n][m]=i11rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
                            rho11[j][j][n][m]= 2*rho11[j-1][j][n][m]- rho11[j-2][j][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho21[j-1][j][n][m]
                                              - D(Nt,Nx,j-1,j,n,m,rho11,ax) + rho11[j-1][j][n][m]
                                              + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*rho11[j-1][j][n][m]
                                              + u/4*(3* (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho11[j-1][j][n][m]+ (f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho11[j-1][j][n][m]+
                                                      2*f12[j-1][j-1][n][n]*rho21[j-1][j][n][m]) +

                                                    i11rho[j-1][j][n][m]

                                              );

                              i12[j-1][j][n][m]= i12a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
                              f12[j][j][n][m]= 2*f12[j-1][j][n][m]- f12[j-2][j][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f22[j-1][j][n][m]
                            - D(Nt,Nx,j-1,j,n,m,f12,ax) + f12[j-1][j][n][m]
                            + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*f12[j-1][j][n][m]
                            + u/4*(3*(f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f12[j-1][j][n][m]+ (f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f12[j-1][j][n][m]+ 2* f12[j-1][j-1][n][n]*f22[j-1][j][n][m])+

                            i12[j-1][j][n][m]
                          );
            i12rho[j-1][j][n][m]=i12rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
            rho12[j][j][n][m]= 2*rho12[j-1][j][n][m]- rho12[j-2][j][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho22[j-1][j][n][m]
                              - D(Nt,Nx,j-1,j,n,m,rho12,ax) + rho12[j-1][j][n][m]
                              + u/4*(3*pow(phi1[j-1][n],2)+pow(phi2[j-1][n],2))*rho12[j-1][j][n][m]
                              + u/4*(3*(f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho12[j-1][j][n][m]+ (f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho12[j-1][j][n][m]+ 2* f12[j-1][j-1][n][n]*rho22[j-1][j][n][m])+

                              i12rho[j-1][j][n][m]
                            );

i21[j-1][j][n][m]= i21a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
f21[j][j][n][m]= 2*f21[j-1][j][n][m] - f21[j-2][j][n][m] - pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f11[j-1][j][n][m]
                  - D(Nt,Nx,j-1,j,n,m,f21,ax) + f21[j-1][j][n][m]
                  + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*f21[j-1][j][n][m]
                  + u/4*((f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f21[j-1][j][n][m]+ 3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f21[j-1][j][n][m]+
                          2*f21[j-1][j-1][n][n]*f11[j-1][j][n][m])+

                    i21[j-1][j][n][m]
                  );
i21rho[j-1][j][n][m]=i21rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
rho21[j][j][n][m]= 2*rho21[j-1][j][n][m]-rho21[j-2][j][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho11[j-1][j][n][m]
                  - D(Nt,Nx,j-1,j,n,m,rho21,ax) + rho21[j-1][j][n][m]
                  + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*rho21[j-1][j][n][m]
                  + u/4*((f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho21[j-1][j][n][m]+ 3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho21[j-1][j][n][m]+
                          2*f21[j-1][j-1][n][n]*rho11[j-1][j][n][m])+

                    i21rho[j-1][j][n][m]
                  );

i22[j-1][j][n][m]= i22a(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
f22[j][j][n][m]= 2*f22[j-1][j][n][m]- f22[j-2][j][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*f12[j-1][j][n][m]
                  - D(Nt,Nx,j-1,j,n,m,f22,ax) + f22[j-1][j][n][m]
                  + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*f22[j-1][j][n][m]
                  + u/4*(2*f21[j-1][j-1][n][n]*f12[j-1][j][n][m]+ (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*f22[j-1][j][n][m]+
                          3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*f22[j-1][j][n][m])+

                          i22[j-1][j][n][m]
                        );
i22rho[j-1][j][n][m]=i22rhoa(at,ax,Nt,Nx,u,f11,f12,f21,f22,rho11,rho12,rho21,rho22,j-1,j, n,m,phi1,phi2);
rho22[j][j][n][m]= 2*rho22[j-1][j][n][m] - rho22[j-2][j][n][m]- pow(at,2)* (u/2*phi1[j-1][n]*phi2[j-1][n]*rho12[j-1][j][n][m]
                  - D(Nt,Nx,j-1,j,n,m,rho22,ax)+ rho22[j-1][j][n][m]
                  + u/4*(pow(phi1[j-1][n],2)+3*pow(phi2[j-1][n],2))*rho22[j-1][j][n][m]
                  + u/4*(2*f21[j-1][j-1][n][n]*rho12[j-1][j][n][m]+ (f11[j-1][j-1][n][n] + ren1(Nx,ax,sigma))*rho22[j-1][j][n][m]+
                          3*(f22[j-1][j-1][n][n] + ren2(Nx,ax,sigma))*rho22[j-1][j][n][m])+

                          i22rho[j-1][j][n][m]
                        );
      /*end of paralleliz*/    }}
        }

        for (n=0; n<Nx;n++){
          for (m=0;m<Nx;m++){
            f11[j][j][n][m]= (f11[j][j][n][m]+f11[j][j][m][n])*.5;
            f22[j][j][n][m]= (f22[j][j][n][m]+f22[j][j][m][n])*.5;
            f12[j][j][n][m]=(f12[j][j][n][m] + f21[j][j][m][n])*.5;
            f21[j][j][n][m]=(f21[j][j][n][m]+f12[j][j][m][n])*.5;
            rho11[j][j][n][m]=0;
            rho22[j][j][n][m]=0;
            rho12[j][j][n][m]=0;
            rho21[j][j][n][m]=0;
       }}

      /* Recovering the missing elements in the previous cycle(the ones requiring the diagonal term)*/

    }
/*Modulus behaviour printed above, f11(t,x;0,0), f12(t,x,0,0), f22(t,x,0,0), f11(t,x:t,x), f12(t,x;t,x),f22(t,x,t,x),rho11(t,x,t,x),rho12(t,x,t,x),rho12(t,xn,t,xn+1), rho12(t,x,0,0) */



printf("mod={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=sqrt(pow(phi1[p][m],2)+pow(phi2[p][m],2))/(sqrt(2));
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f1100={");
    for (p = 0; p <Nx; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f11[0][0][p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nx-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho1101xy={");
    for (p = 0; p <Nx; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho11[0][1][p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nx-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f2200={");
    for (p = 0; p <Nx; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f22[0][0][p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nx-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f1101xy={");
    for (p = 0; p <Nx; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f11[0][1][p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nx-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");


printf("f1111xy={");
    for (p = 0; p <Nx; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f11[1][1][p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nx-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f1122xy={");
    for (p = 0; p <Nx; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f11[2][2][p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nx-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");


printf("f11tx00={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f11[p][0][m][0];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f12tx00={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f12[p][0][m][Nx/2];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho11tt1xx={");
    for (p = 0; p <Nt-1; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho11[p][p+1][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-2){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f21tx00={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f21[p][0][m][Nx/2];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f22tx00={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f22[p][0][m][Nx/2];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f11txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f11[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f12txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f12[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f21txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f21[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("f22txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=f22[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho11txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho11[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho12tx00={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho12[p][0][m][Nx/2];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho11tx00={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho11[p][0][m][Nx/2];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho12txtx1={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=1;m<Nx-1;m++) {
        a=rho12[p][p][m][m+1];
        if (m!= Nx-2){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("phi1={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=phi1[p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("phi2={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=phi2[p][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho12txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho12[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("rho21txtx={");
    for (p = 0; p <Nt; p++ ) {
      printf("{");
      for (m=0;m<Nx;m++) {
        a=rho21[p][p][m][m];
        if (m!= Nx-1){printf("%f," , a );}
        else {printf("%f" , a );}
        }
        if (p< Nt-1){
        printf("},");}
        else {printf("}");}
      }
printf("};\n");

printf("(*u=%f, Nt=%i , at=%f , ax=%f, Nx=%i, f110=tsunami, sigma=%f, f220=delta ,phi1=%f , phi2=%f*)\n", u, Nt, at, ax, Nx, sigma,phi1[0][1],phi2[0][1]);

 /*Clearing the heap before exiting */
 free(allElements11);
 free(allElements12);
 free(allElements21);
 free(allElements22);
 free(rallElements11);
 free(rallElements12);
 free(rallElements21);
 free(rallElements22);
 free(rrallElements11);
 free(rrallElements12);
 free(rrallElements21);
 free(rrallElements22);
 free(rrrallElements11);
 free(rrrallElements12);
 free(rrrallElements21);
 free(rrrallElements22);

 free(rrrrallElements11);
 free(rrrrallElements12);
 free(rrrrallElements21);
 free(rrrrallElements22);
 free(rrrrrallElements11);
 free(rrrrrallElements12);
 free(rrrrrallElements21);
 free(rrrrrallElements22);
 for(i1=0; i1 < r1;i1++){
   for(j1=0;j1<r2;j1++){
     free(f11[i1][j1]);
     free(f12[i1][j1]);
     free(f21[i1][j1]);
     free(f22[i1][j1]);
     free(rho11[i1][j1]);
     free(rho12[i1][j1]);
     free(rho21[i1][j1]);
     free(rho22[i1][j1]);

     free(i11[i1][j1]);
     free(i12[i1][j1]);
     free(i21[i1][j1]);
     free(i22[i1][j1]);
     free(i11rho[i1][j1]);
     free(i12rho[i1][j1]);
     free(i21rho[i1][j1]);
     free(i22rho[i1][j1]);

     free(ii11[i1][j1]);
     free(ii12[i1][j1]);
     free(ii21[i1][j1]);
     free(ii22[i1][j1]);
     free(ii11rho[i1][j1]);
     free(ii12rho[i1][j1]);
     free(ii21rho[i1][j1]);
     free(ii22rho[i1][j1]);
 }}
 for(i1=0; i1 < r1;i1++){
   free(f11[i1]);
   free(f12[i1]);
   free(f21[i1]);
   free(f22[i1]);
   free(rho11[i1]);
   free(rho12[i1]);
   free(rho21[i1]);
   free(rho22[i1]);

   free(i11[i1]);
   free(i12[i1]);
   free(i21[i1]);
   free(i22[i1]);
   free(i11rho[i1]);
   free(i12rho[i1]);
   free(i21rho[i1]);
   free(i22rho[i1]);
   free(ii11[i1]);
   free(ii12[i1]);
   free(ii21[i1]);
   free(ii22[i1]);
   free(ii11rho[i1]);
   free(ii12rho[i1]);
   free(ii21rho[i1]);
   free(ii22rho[i1]);
 }
 free(f11);
 free(f12);
 free(f21);
 free(f22);
 free(rho11);
 free(rho12);
 free(rho21);
 free(rho22);

 free(i11);
 free(i12);
 free(i21);
 free(i22);
 free(i11rho);
 free(i12rho);
 free(i21rho);
 free(i22rho);

 free(ii11);
 free(ii12);
 free(ii21);
 free(ii22);
 free(ii11rho);
 free(ii12rho);
 free(ii21rho);
 free(ii22rho);


      return 0;
    }
