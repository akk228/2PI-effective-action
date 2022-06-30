//
//  main.cpp
//  BEC2l
//
//  Created by Andrei Kovtun on 12.11.2019.
//  Copyright © 2019 anonymousfaggots. All rights reserved.
//
//
//  main.cpp
//  1loopselfener
//
//  Created by Andrei Kovtun on 11.11.2019.
//  Copyright © 2019 anonymousfaggots. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

#define Pi 3.14159265358979
int N = 2;

double F0 (int a, int b, double m, double omega, double lambda, double phi0, double L, double p){
    
    double gp, gn, fp, fn, f, z, sgn1(1), sgn2(1);
    
    z = 0.5*pow(phi0,2);
    
    gp = sqrt( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda + pow(omega,2) + 0.5*sqrt( pow(lambda*z,2) + 16.*pow(omega,2)*( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda ) ) );
    gn = sqrt( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda + pow(omega,2) - 0.5*sqrt( pow(lambda*z,2) + 16.*pow(omega,2)*( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda ) ) );
    
    if (a == 0 && b == a) {
        sgn1 = 1;
        sgn2 = 1;
    }
    else{
        if (b == a){
            sgn1 = -1;
            sgn2 = 1;
        }
        else{
            sgn1 = 0;
            sgn2 = 0;
        }
    };
    
    fp = - sgn2*( pow(m,2) + pow(2*Pi*p/L,2) - pow(gp,2) + z*lambda - pow(omega,2) ) + sgn1*0.5*z*lambda;
    fn = - sgn2*( pow(m,2) + pow(2*Pi*p/L,2) - pow(gn,2) + z*lambda - pow(omega,2) ) + sgn1*0.5*z*lambda;
    
    f = ( 0.5*fp/gp - 0.5*fn/gn )/( pow(gp,2) - pow(gn,2) );
    
    return f;
}

double dF0 (int a, int b, double m, double omega, double lambda, double phi0, double L, double p){
    
    double gp, gn, fp, fn, f, z;
    
    z = 0.5*pow(phi0,2);
    
    gp = sqrt( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda + pow(omega,2) + 0.5*sqrt( pow(lambda*z,2) + 16.*pow(omega,2)*( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda ) ) );
    gn = sqrt( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda + pow(omega,2) - 0.5*sqrt( pow(lambda*z,2) + 16.*pow(omega,2)*( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda ) ) );
    
    if ( b == a ) {
        fp = 0;
        fn = 0;
    }
    else{
        if (a == 0){
         
            fp = - 2.*pow(gp,2)*omega - omega*( pow(m,2) + pow(2*Pi*p/L,2) - pow(gp,2) + z*lambda - pow(omega,2) ) - 0.5*omega*z*lambda;
            fn = - 2.*pow(gn,2)*omega - omega*( pow(m,2) + pow(2*Pi*p/L,2) - pow(gn,2) + z*lambda - pow(omega,2) ) - 0.5*omega*z*lambda;
            
        }
        else{
            
            fp = 2.*pow(gp,2)*omega + omega*( pow(m,2) + pow(2*Pi*p/L,2) - pow(gp,2) + z*lambda - pow(omega,2) ) - 0.5*omega*z*lambda;
            fn = 2.*pow(gn,2)*omega + omega*( pow(m,2) + pow(2*Pi*p/L,2) - pow(gn,2) + z*lambda - pow(omega,2) ) - 0.5*omega*z*lambda;
            
        }
    }
    
    f = ( 0.5*fp/gp - 0.5*fn/gn )/( pow(gp,2) - pow(gn,2) );
    
    return f;
}

double ddF0 (int a, int b, double m, double omega, double lambda, double phi0, double L, double p, double dt){
    
    double gp, gn, fp, fn, f, z;
    
    z = 0.5*pow(phi0,2);
    
    gp = sqrt( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda + pow(omega,2) + 0.5*sqrt( pow(lambda*z,2) + 16.*pow(omega,2)*( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda ) ) );
    gn = sqrt( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda + pow(omega,2) - 0.5*sqrt( pow(lambda*z,2) + 16.*pow(omega,2)*( pow(m,2) + pow(2*Pi*p/L,2) + z*lambda ) ) );
    
    if ( a == b ) {
        if ( a == 0 ){
            fp = 0.5*pow(dt,2)*( -2*lambda*z*pow(omega,2) );
            fn = 0.5*pow(dt,2)*( -2*lambda*z*pow(omega,2) );
        }
        else{
            fp = - 0.5*pow(dt,2)*( -2*lambda*z*pow(omega,2) );
            fn = - 0.5*pow(dt,2)*( -2*lambda*z*pow(omega,2) );
        }
        
    }
    else{
        fp = - dt*lambda*z*omega;
        fn = - dt*lambda*z*omega;
    }
    
    f = ( 0.5*fp/gp - 0.5*fn/gn )/( pow(gp,2) - pow(gn,2) );
    
    return f;
}

double Gphi( int a, int x0, int Np, double dt, double L, double lambda, double *****F, double *****rho,double **phi){
    double kek(0), coef;
    
    for ( int y0 = 0; y0 <= x0; y0++){
        for ( int p = - ( Np - 1 ); p < Np; p++){
            for ( int q = - ( Np - 1 ); q < Np; q++){
                for ( int b = 0; b < N; b++){
                    for ( int c = 0; c < N; c++){
                        for ( int d = 0; d < N; d++){
                            
                            if (y0 == 0 or y0 == x0) {
                                coef = 0.5;
                            }else{
                                coef = 1.0;
                            };
                            
                            kek = kek +
                                    coef*(
                                          
                                    ( F[c][d][x0][y0][abs(p)]*F[c][d][x0][y0][ abs((p-q)) % Np ] - 0.25*rho[c][d][x0][y0][abs(p)]*rho[c][d][x0][y0][ abs((p-q)) % Np ] )*rho[a][b][x0][y0][abs(q)]*phi[b][y0] +
                            
                                    2.*F[c][d][x0][y0][abs(p)]*rho[c][d][x0][y0][ abs((p-q)) % Np ]*F[a][b][x0][y0][abs(q)]*phi[b][y0] +
                            
                                    2.*( rho[a][b][x0][y0][abs(p)]*F[c][b][x0][y0][ abs((p-q)) % Np ] + F[a][b][x0][y0][abs(p)]*rho[c][b][x0][y0] [abs((p-q)) % Np ] )*F[c][d][x0][y0][abs(q)]*phi[d][y0] -
                            
                                    0.5*rho[a][b][x0][y0][abs(p)]*rho[c][b][x0][y0][ abs((p-q)) % Np ]*rho[c][d][x0][y0][abs(q)]*phi[d][y0]
                                    
                                          );
                            
                        };
                    };
                };
            };
        };
    };
    
    return 0.125*pow(lambda,2)*kek*dt/pow(L,2);
    
};

double memF( int a, int b, int x0, int y0, int p, int Np, double dt, double L, double lambda, double *****F, double *****rho,double **phi){
    double sumFrho(0), sumrhoF(0), coef;
    
    //Sigma^{F}_{a,c})*rho_{c,b}
    
    for ( int z0 = 0; z0 <= y0; z0++){
        for ( int q = - ( Np - 1 ); q < Np; q++){
            for ( int m = 0; m < N; m++){
                for ( int n = 0; n < N; n++){
                    for ( int f = 0; f < N; f++){
                            
                        if (z0 == 0 or z0 == y0) {
                            coef = 0.5;
                        }else{
                            coef = 1.0;
                        };
                            
                        sumFrho = sumFrho -
                                sumFrho*(
                                         ( F[m][n][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] - 0.25*rho[m][n][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] )*phi[a][x0]*phi[f][z0]*rho[f][b][z0][y0][p] +
                                      2*(
                                         ( F[m][n][x0][z0][abs(p-q) % Np]*F[a][f][x0][z0][abs(q)] - 0.25*rho[m][n][x0][z0][abs(p-q) % Np]*rho[a][f][x0][z0][abs(q)] )*phi[m][x0]*phi[n][z0]*rho[f][b][z0][y0][p] +
                                         ( F[m][f][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] - 0.25*rho[m][f][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] )*phi[a][x0]*phi[n][z0]*rho[f][b][z0][y0][p] +
                                         ( F[a][n][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] - 0.25*rho[a][n][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] )*phi[m][x0]*phi[f][z0]*rho[f][b][z0][y0][p] +
          /*Check!!!!!!!!!!!!!!!*/       ( F[a][m][x0][z0][abs(p-q) % Np]*F[n][f][x0][z0][abs(q)] - 0.25*rho[a][m][x0][z0][abs(p-q) % Np]*rho[n][f][x0][z0][abs(q)] )*phi[n][x0]*phi[m][z0]*rho[f][b][z0][y0][p]
                                        )
                                     );
                            
                    };
                };
            };
        };
    };
    
    //Sigma^{rho}_{a,c})*F_{c,b}
    
    for ( int z0 = 0; z0 <= x0; z0++){
        for ( int q = - ( Np - 1 ); q < Np; q++){
            for ( int m = 0; m < N; m++){
                for ( int n = 0; n < N; n++){
                    for ( int f = 0; f < N; f++){
                            
                        if (z0 == 0 or z0 == x0) {
                            coef = 0.5;
                        }else{
                            coef = 1.0;
                        };
                            
                        sumrhoF = sumrhoF -
                                coef*(
                                         ( F[m][n][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] + rho[m][n][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] )*phi[a][x0]*phi[f][z0]*F[f][b][z0][y0][p] +
                                      2*(
                                         ( F[m][n][x0][z0][abs(p-q) % Np]*rho[a][f][x0][z0][abs(q)] + rho[m][n][x0][z0][abs(p-q) % Np]*F[a][f][x0][z0][abs(q)] )*phi[m][x0]*phi[n][z0]*F[f][b][z0][y0][p] +
                                         ( F[m][f][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] + rho[m][f][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] )*phi[a][x0]*phi[n][z0]*F[f][b][z0][y0][p] +
                                         ( F[a][n][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] + rho[a][n][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] )*phi[m][x0]*phi[f][z0]*F[f][b][z0][y0][p] +
           /*Check!!!!!!!!!!!!!!!*/      ( F[a][m][x0][z0][abs(p-q) % Np]*rho[n][f][x0][z0][abs(q)] + rho[a][m][x0][z0][abs(p-q) % Np]*F[n][f][x0][z0][abs(q)] )*phi[n][x0]*phi[m][z0]*F[f][b][z0][y0][p]
                                        )
                                     );
                            
                    };
                };
            };
        };
    };
    
    return 0.125*pow(lambda,2)*( - sumrhoF + sumFrho )*dt/L;
    
};

double memrho( int a, int b, int x0, int y0, int p, int Np, double dt, double L, double lambda, double *****F, double *****rho,double **phi){
    double sumrhorho(0), coef, sgn;
    int z1, z2;
    
    if (x0 > y0){
        sgn = 1;
        z2 = x0;
        z1 = y0;
    }else{
        sgn = -1;
        z2 = y0;
        z1 = x0;
    }
    for ( int z0 = z1; z0 <= z2; z0++){
        for ( int q = - ( Np - 1 ); q < Np; q++){
            for ( int m = 0; m < N; m++){
                for ( int n = 0; n < N; n++){
                    for ( int f = 0; f < N; f++){
                        
                        if (z0 == z1 or z0 == z2) {
                            coef = 0.5;
                        }else{
                            coef = 1.0;
                        };
                        
                        sumrhorho = sumrhorho -
                            sgn*coef*(
                                     ( F[m][n][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] + rho[m][n][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] )*phi[a][x0]*phi[f][z0]*rho[f][b][z0][y0][p] +
                                  2*(
                                     ( F[m][n][x0][z0][abs(p-q) % Np]*rho[a][f][x0][z0][abs(q)] + rho[m][n][x0][z0][abs(p-q) % Np]*F[a][f][x0][z0][abs(q)] )*phi[m][x0]*phi[n][z0]*rho[f][b][z0][y0][p] +
                                     ( F[m][f][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] + rho[m][f][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] )*phi[a][x0]*phi[n][z0]*rho[f][b][z0][y0][p] +
                                     ( F[a][n][x0][z0][abs(p-q) % Np]*rho[m][n][x0][z0][abs(q)] + rho[a][n][x0][z0][abs(p-q) % Np]*F[m][n][x0][z0][abs(q)] )*phi[m][x0]*phi[f][z0]*rho[f][b][z0][y0][p] +
        /*Check!!!!!!!!!!!!!!!*/     ( F[a][m][x0][z0][abs(p-q) % Np]*rho[n][f][x0][z0][abs(q)] + rho[a][m][x0][z0][abs(p-q) % Np]*F[n][f][x0][z0][abs(q)] )*phi[n][x0]*phi[m][z0]*rho[f][b][z0][y0][p]
                                     )
                                     );
                        
                };
            };
        };
    };
};
        return - 0.125*pow(lambda,2)*sumrhorho*dt/L;
};

int main() {
    
    double **phi, *****F, *****rho;
    double dx(0.5), dt(0.1), L, m(1), dm2, phi0, lambda, omega, SE[2][2];
    int Nt(80), Np(40);
   
    
    lambda = 0.1;
    omega  = 1.5;
    phi0   = sqrt(2*25.1417807081529893765);//2.*sqrt( ( pow(omega,2) - pow(m,2) )/lambda );;
    L = Np*dx;
    //cout << 0.5*pow(phi0,2) << endl;
    
    //Initzialization of functions
    
    phi = new double* [2];
    F   = new double****[2];
    rho = new double****[2];
    
    for (int a = 0; a < 2; a++){
        
        phi[a] = new double [Nt];
        F[a]   = new double***[2];
        rho[a] = new double***[2];
        
        for (int b = 0; b < 2; b++){
        
            F[a][b]   = new double**[Nt];
            rho[a][b] = new double**[Nt];
    
            for (int t = 0; t<Nt; t++) {
        
                F[a][b][t]   = new double*[Nt];
                rho[a][b][t] = new double*[Nt];
            
                for (int tau = 0; tau < Nt; tau++) {
                
                    F[a][b][t][tau]   = new double [Np];
                    rho[a][b][t][tau] = new double [Np];
                
                };
            };
        };
    };
    
    //Initial conditions
    
    //Renormalization
    
    dm2 = - lambda*0.5/m/L;
    
    for (int p = 1; p < Np; p++){
        dm2 = dm2 - 2.*(lambda/L)*0.5/sqrt( pow(2*Pi*p/L,2) + pow(m,2) );
    };
    
    //Self Energies
    
    for ( int a = 0; a < 2; a++ ){
        for ( int b = 0; b < 2; b++ ){
            for ( int p = 0; p < Np; p++){
                if ( p == 0){
                    SE[a][b] = (1./L)*F0 (a, b, m, omega, lambda, phi0, L, 0);
                }else{
                    SE[a][b] = SE[a][b] + 2*(1./L)*F0(a, b, m, omega, lambda, phi0, L, p);
                };
            };
        };
    };
    
    cout << dm2 + 0.25*lambda*( 3*SE[0][0] + SE[1][1] ) << endl;
    
    //phi
    
    phi[0][0] = phi0;
    phi[1][0] = 0;
    
    phi[0][1] = phi[0][0] + 0*dt -          0.5*pow(dt,2)*( ( pow(m,2) + dm2 + 0.25*lambda*( pow(phi0,2) + SE[0][0] + SE[1][1] ) )*phi[0][0] + 0.5*lambda*( SE[0][0]*phi[0][0] + SE[0][1]*phi[1][1] ) );
    phi[1][1] = phi[1][0] - omega*phi0*dt - 0.5*pow(dt,2)*( ( pow(m,2) + dm2 + 0.25*lambda*( pow(phi0,2) + SE[0][0] + SE[1][1] ) )*phi[1][0] + 0.5*lambda*( SE[1][0]*phi[0][0] + SE[1][1]*phi[1][1] ) );
    
    //Green's functions
    
    for (int p = 0; p < Np; p++){
        
        for ( int a = 0; a < 2; a++ ){
            for ( int b = 0; b < 2; b++ ){
                //F
                /*if ( p == 0 ){
                    F[a][b][0][0][p] = 0;
                    F[a][b][1][0][p] = 0;
                    F[a][b][1][1][p] = 0;
                }else{*/
                    F[a][b][0][0][p] = F0 (a, b, m, omega, lambda, phi0, L, p);
                    F[a][b][1][0][p] = F[a][b][0][0][p] + dt*dF0 (a, b, m, omega, lambda, phi0, L, p) -
                                                        0.5*pow(dt,2)*( ( pow(2*Pi*p/L,2) + pow(m,2) + 0.25*lambda*pow(phi0,2) )*F[a][b][0][0][p] + 0.5*lambda* ( phi[a][0]*phi[0][0]*F[0][b][0][0][p] + phi[a][0]*phi[1][0]*F[1][b][0][0][p] ) );
                    F[a][b][1][1][p] = F[a][b][0][0][p] + ddF0 (a, b, m, omega, lambda, phi0, L, p, dt);
                /*};*/
                
                F[b][a][0][1][p] = F[a][b][1][0][p];
                
                //rho
                rho[a][b][0][0][p] = 0;
                
                if ( a == b ){
                    rho[a][b][1][0][p] = rho[a][b][0][0][p] + dt;
                }else{
                    rho[a][b][1][0][p] = rho[a][b][0][0][p];
                };
                
                rho[b][a][0][1][p] = - rho[a][b][1][0][p];
                
                rho[a][b][1][1][p] = 0;
                
            };
        };
    };
    
    //SOOOOLVINGGGG BOOYAKA BOOYAKA
    
    for (int t = 2; t < Nt; t++){
        
        //Self Energies
        
        for ( int a = 0; a < 2; a++ ){
            for ( int b = 0; b < 2; b++ ){
                for ( int p = 0; p < Np; p++){
                    if ( p == 0){
                        SE[a][b] = (1./L)*F[a][b][t-1][t-1][0];
                    }else{
                        SE[a][b] = SE[a][b] + 2.*(1./L)*F[a][b][t-1][t-1][p];
                    };
                };
            };
        };
        
        for (int a = 0; a < 2; a++){
            //phi
            phi[a][t] = 2*phi[a][t-1] - phi[a][t-2]
                        - pow(dt,2)*(
                                     (pow(m,2) + dm2 + 0.25*lambda*( pow(phi[0][t-1],2) + pow(phi[1][t-1],2) + SE[0][0] + SE[1][1] ) )*phi[a][t-1] + 0.5*lambda*( SE[a][0]*phi[0][t-1] + SE[a][1]*phi[1][t-1] )
                                    ) +
                          pow(dt,2)*Gphi(a, t-1, Np, dt, L, lambda, F, rho, phi);
        };
       
        for (int tau = 0; tau < t; tau++){
            for (int p = 0; p < Np; p++){
                //F
                for ( int a = 0; a < 2; a++ ){
                    for ( int b = 0; b < 2; b++ ){
                        F[a][b][t][tau][p] = 2.*F[a][b][t-1][tau][p] - F[a][b][t-2][tau][p] -
                                pow(dt,2)*(
                                            ( pow(2*Pi*p/L,2) + pow(m,2) + dm2 + 0.25*lambda*( pow(phi[0][t-1],2) + pow(phi[1][t-1],2) ) )*F[a][b][t-1][tau][p] + 0.5*lambda*( phi[a][t-1]*phi[0][t-1]*F[0][b][t-1][tau][p] + phi[a][t-1]*phi[1][t-1]*F[1][b][t-1][tau][p] ) +
                                                    0.25*lambda*( SE[0][0] + SE[1][1] )*F[a][b][t-1][tau][p] + 0.5*lambda*( SE[a][0]*F[0][b][t-1][tau][p] + SE[a][1]*F[1][b][t-1][tau][p] )
                                           ) +
                                pow(dt,2)*memF( a, b, t-1, tau, p, Np, dt, L, lambda, F, rho, phi);
                    };
                };
                //rho
                for ( int a = 0; a < 2; a++ ){
                    for ( int b = 0; b < 2; b++ ){
                
                        rho[a][b][t][tau][p] = 2.*rho[a][b][t-1][tau][p] - rho[a][b][t-2][tau][p] -
                            pow(dt,2)*(
                                        ( pow(2*Pi*p/L,2) + pow(m,2) + dm2 + 0.25*lambda*( pow(phi[0][t-1],2) + pow(phi[1][t-1],2) ) )*rho[a][b][t-1][tau][p] + 0.5*lambda*( phi[a][t-1]*phi[0][t-1]*rho[0][b][t-1][tau][p] + phi[a][t-1]*phi[1][t-1]*rho[1][b][t-1][tau][p] ) +
                                                    0.25*lambda*( SE[0][0] + SE[1][1] )*rho[a][b][t-1][tau][p] + 0.5*lambda*( SE[a][0]*rho[0][b][t-1][tau][p] + SE[a][1]*rho[1][b][t-1][tau][p] )
                                       ) +
                            pow(dt,2)*memrho( a, b, t-1, tau, p, Np, dt, L, lambda, F, rho, phi);
                    };
                };
            
            };
        };
        
        for ( int tau = 0; tau < t; tau++ ){
            for ( int p = 0; p < Np; p++ ){
                for ( int a = 0; a < 2; a++ ){
                    for ( int b = 0; b < 2; b++ ){
                        F[b][a][tau][t][p] = F[a][b][t][tau][p];
                        rho[b][a][tau][t][p] = - rho[a][b][t][tau][p];
                    }
                }
            }
        }
        
        //Diagonal elements
        for (int p = 0; p < Np; p++){
                
            for ( int a = 0; a < 2; a++ ){
                for ( int b = 0; b < 2; b++ ){
            
                F[a][b][t][t][p] = 2.*F[a][b][t-1][t][p] - F[a][b][t-2][t][p] -
                                    pow(dt,2)*(
                                                 ( pow(2*Pi*p/L,2) + pow(m,2) + dm2 + 0.25*lambda*( pow(phi[0][t-1],2) + pow(phi[1][t-1],2) ) )*F[a][b][t-1][t][p] + 0.5*lambda* ( phi[a][t-1]*phi[0][t-1]*F[0][b][t-1][t][p] + phi[a][t-1]*phi[1][t-1]*F[1][b][t-1][t][p] ) +
                                                         0.25*lambda*( SE[0][0] + SE[1][1] )*F[a][b][t-1][t][p] + 0.5*lambda*( SE[a][0]*F[0][b][t-1][t][p] + SE[a][1]*F[1][b][t-1][t][p] )
                                               ) +
                                    pow(dt,2)*memF( a, b, t-1, t, p, Np, dt, L, lambda, F, rho, phi);
                    
                rho[a][b][t][t][p] = 0;
                    
                };
            };
            
        };
        
        for (int p = 0; p < Np; p++){
            
            for ( int a = 0; a < 2; a++ ){
                for ( int b = 0; b < 2; b++ ){
                          F[a][b][t][t][p] = 0.5*( F[a][b][t][t][p] + F[b][a][t][t][p] );
                }
            }
        }
        
    };
    
    //Output
    ofstream FG;
    
    FG.open ("2l.txt");
    
    FG << "dt=" << dt << "\n";
    FG << "dx=" << dx << "\n";
    FG << "Nt=" << Nt << "\n";
    FG << "Np=" << Np << "\n";
    FG << "phi0=" << phi0 << "\n";
    FG << "omega=" << omega << "\n";
    FG << "Len=" << L << "\n";
    FG << "lambda=" << lambda << "\n";
    
    //phi1
    
    FG << "phi1=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        
        FG << phi[0][t];
        if ( t != Nt - 1 ){
            FG << ",";
        }
        
    };
    
    FG << "};\n";
    
    //phi2
    
    FG << "phi2=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        
        FG << phi[1][t];
        if ( t != Nt - 1 ){
            FG << ",";
        }
        
    };
    
    FG << "};\n";
    
    //F11
    
    FG << "f11t0p=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        
         FG << "{";
        
        for (int p = 0; p < Np; p++){
            
            FG << F[0][0][t][0][p];
            if ( p != Np - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 1 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    
    //F12
    
    FG << "f12t0p=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        
         FG << "{";
        
        for (int p = 0; p < Np; p++){
            
            FG << F[0][1][t][0][p];
            if ( p != Np - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 1 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    
    //F12
    
    FG << "f21t0p=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        
         FG << "{";
        
        for (int p = 0; p < Np; p++){
            
            FG << F[1][0][t][0][p];
            if ( p != Np - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 1 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    
    //rho11
    
    FG << "rho11t0p=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        
         FG << "{";
        
        for (int p = 0; p < Np; p++){
            
            FG << rho[0][0][t][0][p];
            if ( p != Np - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 1 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    
    //drho11
    
    FG << "drho11t0p=" << "\n" << "{";
    
    for (int t = 1; t < Nt-1; t++){
        
         FG << "{";
        
        for (int p = 0; p < Np; p++){
            
            FG << ( rho[0][0][t+1][t][p] - rho[0][0][t-1][t][p])/2./dt;
            if ( p != Np - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 2 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    
    
    double q;
    //f12tt
    
    FG << "f12ttau=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        FG << "{";
        for (int tau = 0; tau < Nt; tau++){
            q = F[0][1][t][tau][0]/L;
            for (int p = 1; p < Np; p++){
                q = q + 2*F[0][1][t][tau][p]/L;
            };
            FG << q;
            if ( tau != Nt - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 1 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    //f21tt
    FG << "f21ttau=" << "\n" << "{";
    
    for (int t = 0; t < Nt; t++){
        FG << "{";
        for (int tau = 0; tau < Nt; tau++){
            q = F[1][0][t][tau][0]/L;
            for (int p = 1; p < Np; p++){
                q = q + 2*F[1][0][t][tau][p]/L;
            };
            FG << q;
            if ( tau != Nt - 1 ){
                FG << ",";
            };
            
        };
        
        FG << "}\n";
        if ( t != Nt - 1 ){
            FG << ",";
        };
    };
        
    FG << "};\n";
    
    FG.close();
    
    return 0;
}


