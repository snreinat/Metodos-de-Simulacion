
// SIMULACIÓN DEL MOVIMIENTO DE N PARTÍCULAS BAJO FUERZAS DE REPULSIÓN ELECTRÓSTÁTICA Y FUERZAS CENTRALES

#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double K=1000, Kcentral=1.0, gam=0.5;
const double Lx=200, Ly=200;
const int N=6;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R,Q; 
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); 
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  // cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
}


//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Particula);
  void CalculeFuerzaEntre(Cuerpo & Particula1, Cuerpo & Particula2);
  void FuerzaCentral(Cuerpo & Particula1);
  void FuerzaViscosa(Cuerpo & Particula1);
};

void Colisionador::CalculeFuerzas(Cuerpo * Particula){
  int i,j;
  
   //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++)
    Particula[i].BorreFuerza();
  
  //--- Sumar la fuerza central ---
  for(i=0;i<N;i++){
    FuerzaCentral(Particula[i]);// Calcula la fuerza central
    FuerzaViscosa(Particula[i]);//Calcula la fuerza viscosa
  }
  //--- Calcular Fuerzas entre pares de particulas ---
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++){
      CalculeFuerzaEntre(Particula[i], Particula[j]);
    }
}


//Fuerza de repulsión electrostática

void Colisionador::CalculeFuerzaEntre(Cuerpo & Particula1,Cuerpo & Particula2){
  vector3D r21=Particula2.r-Particula1.r;
  double d=r21.norm(); //Norma del vector r_ij
  vector3D n=r21*(1.0/d);//vector r_ij unitario
  
  vector3D F12=((K*Particula1.Q*Particula2.Q)/(d*d))*n; //Fuerza que la partícula i le hace a la j
  Particula2.AdicioneFuerza(F12);   Particula1.AdicioneFuerza(F12*(-1));
  }



void Colisionador::FuerzaCentral(Cuerpo & Particula1){
  vector3D Fi=(-1)*Kcentral*Particula1.r;
  Particula1.AdicioneFuerza(Fi);
}

void Colisionador::FuerzaViscosa(Cuerpo & Particula1){
  vector3D Fviscosai=(-1)*gam*Particula1.V;
  Particula1.AdicioneFuerza(Fviscosai);
}
//----------------- Funciones de Animacion ----------

void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'Gas2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-110:110]"<<endl;
  cout<<"set yrange[-110:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}

void InicieCuadro(void){
    cout<<"plot 0,0 ";
    // cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    //cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    //cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    //cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
    }

void TermineCuadro(void){
    cout<<endl;
    }

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Particula[N];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=3.0, Q=1;
  int i,ix,iy;
  double T=2*M_PI*sqrt(m0/Kcentral);
  double t,tdibujo,tmax=5*T, tcuadro=tmax/1000,dt=0.001;
  double Theta;
  double Ranillo=40.0;

  
  
  //Inicializar Molécula:
      
 //-------------------------(       x0,   y0,           Vx0,           Vy0, m0, R0)
  for(i=0; i<N; i++){
    Particula[i].Inicie(40*cos(i*2*M_PI/N), 40*sin(i*2*M_PI/N), 0, 0, m0, R0); // La partícula debe estar a una distancia de Ranillo del centro, en el eje x
  }
  
  
 InicieAnimacion(); //Dibujar
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    
    //Dibujar animacion
    if(tdibujo>tcuadro){
  
      InicieCuadro();
      for(i=0; i<N; i++){
	Particula[i].Dibujese();
      }
      TermineCuadro();
      
       tdibujo=0;
    }
    
    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,epsilon);  

  }   

 
  
  return 0;
}
