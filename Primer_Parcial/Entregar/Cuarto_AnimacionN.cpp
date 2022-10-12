// SIMULACIÓN DEL MOVIMIENTO DE N PARTÍCULAS  BAJO UNA FUERZA CENTRAL, UNA FUERZA DISIPATIVA Y UNA FUERZA DE REPULSIÓN

#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double K=1000, Kcentral=1.0, gam=0.5;
const int N=6; //Ingresar el número de partículas que se desean

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
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double Q0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double Q0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0; Q=Q0;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); 
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Particula);
  void FuerzaCentral(Cuerpo & Particula1);
  void FuerzaViscosa(Cuerpo & Particula1);
  void FuerzaRepulsiva(Cuerpo &Particula1, Cuerpo &Particula2);
};


void Colisionador::CalculeFuerzas(Cuerpo * Particula){
  int i,j;
  
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++)
    Particula[i].BorreFuerza();
  
  //--- Sumar la fuerza central ---
  for(i=0;i<N;i++){
    FuerzaCentral(Particula[i]);
    FuerzaViscosa(Particula[i]);
  }

  //---Sumar la fuerza de Repulsión---
  for(i=0;i<N;i++){
    for(j=i+1; j<N; j++){
      FuerzaRepulsiva(Particula[i],Particula[j]);
    }
  }
}

//--- Calcular la fuerza central---
void Colisionador::FuerzaCentral(Cuerpo & Particula1){
  vector3D Fi=(-1)*Kcentral*Particula1.r;
  Particula1.AdicioneFuerza(Fi);
}

//---Calcular la fuerza disipativa---
void Colisionador::FuerzaViscosa(Cuerpo & Particula1){
  vector3D Fviscosai=(-1)*gam*Particula1.V;
  Particula1.AdicioneFuerza(Fviscosai);
}

//---Calcular la fuerza repulsiva---
void Colisionador::FuerzaRepulsiva(Cuerpo &Particula1, Cuerpo &Particula2){
  vector3D r12=Particula2.r-Particula1.r; //Vector que va desde la partícula 1 a la 2
  double d=r12.norm(); //Norma del vector r12
  vector3D n=r12*(1.0/d);//Vector unitario de r12

  vector3D Fij=n*K*Particula1.Q*Particula2.Q/pow(d,2);
  Particula1.AdicioneFuerza((-1)*Fij); Particula2.AdicioneFuerza(Fij);
}

//----------------- Funciones de Animacion ----------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Cuarto_Animacion6.gif'"<<endl;
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
    }

void TermineCuadro(void){
    cout<<endl;
    }

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Particula[N];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=3.0,V0=1.0, Q0=1.0;
  int i,j;
  double T=2*M_PI*sqrt(m0/Kcentral);
  double t,tdibujo,tmax=5*T, tcuadro=tmax/1000,dt=0.01;
  double Theta;
  double Ranillo=40.0;
 
  //Inicializar las partículas: 
  for (i=0;i<N;i++){
    Theta=2*M_PI*ran64.r();
    //--------------------(                     x0,                      y0,Vx0,Vy0, m0, R0, Q0)
    Particula[i].Inicie(Ranillo*cos(2*M_PI*i/N), Ranillo*sin(2*M_PI*i/N),  V0*cos(Theta),  V0*sin(Theta), m0, R0, Q0);
  }
  
  InicieAnimacion(); 
  
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
