// Simular el movimiento de un planeta por PEFRL
#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

const double GM=1.0;

const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Lambda);

class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt,double Coeficiente);
  void Mueva_V(double dt,double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0; 
} 
void Cuerpo::CalculeFuerza(void){
  double aux=GM*m*pow(r.norm2(),-1.5);
  F=-aux*r;
} 
void Cuerpo::Mueva_r(double dt,double Coeficiente){
  r+=V*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double Coeficiente){
  V+=F*(Coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Balon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-12:12]"<<endl;
  cout<<"set yrange[-12:12]"<<endl;
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
  Cuerpo Planeta;
  double r0=10, m0=1;
  double omega,V0,T;
  omega=sqrt(GM/(r0*r0*r0)); T=2*M_PI/omega; V0=r0*omega;
  double t,tdibujo,tmax=3.3*T,tcuadro=T/100,dt=0.1;

  InicieAnimacion(); //Dibujar

  //------------(x0,y0,Vx0,Vy0, m0, R0)
  Planeta.Inicie(r0, 0,  0,0.5*V0, m0,0.5);
  
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      /*
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      */
      cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    }         
    //Muevase por PEFRL (Omelyan 2002)
    Planeta.Mueva_r(dt,Zeta);
    Planeta.CalculeFuerza(); Planeta.Mueva_V(dt,Coeficiente1);
    Planeta.Mueva_r(dt,Chi);
    Planeta.CalculeFuerza(); Planeta.Mueva_V(dt,Lambda);
    Planeta.Mueva_r(dt,Coeficiente2);
    Planeta.CalculeFuerza(); Planeta.Mueva_V(dt,Lambda);
    Planeta.Mueva_r(dt,Chi);
    Planeta.CalculeFuerza(); Planeta.Mueva_V(dt,Coeficiente1);
    Planeta.Mueva_r(dt,Zeta);
  }   

  
  return 0;
}

  
