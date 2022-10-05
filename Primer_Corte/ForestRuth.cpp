#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes globales
const double GM=1.0;

//Constantes de Forest Roth
const double Theta=1/(2-pow(2.0,1.0/3));
const double ThetaU2=Theta/2;
const double UmTheta_U2=(1-Theta)/2;
const double Um2Theta=1-2*Theta;

//Declaración de las clases
class Cuerpo;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  vector3D r,V,F;
  double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0, double Vy0,double Vz0,double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt,double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();};
};
void Cuerpo::Inicie(double x0,double y0,double z0,
		    double Vx0, double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0), m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  double aux=GM*m*pow(r.norm2(),-1.5);
  F=(-aux)*r;
}

void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(dt*Coeficiente);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(dt*Coeficiente/m);
}


//----------- Funciones Globales -----------
int main(){
  Cuerpo Planeta;
  double r0=10,m0=1;
  double omega,V0,T;
  double t,dt=0.0001;

  omega=sqrt(GM/(r0*r0*r0)); V0=omega*r0; T=2*M_PI/omega;
  
  //------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  Planeta.Inicie(r0, 0, 0, 0, V0/2,0,m0,0.5);
  
  for(t=0;t<1.1*T;t+=dt){
    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    //Mover por Forest-Roth
    Planeta.Mueva_r(dt,ThetaU2);
    
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, Theta);
    
    Planeta.Mueva_r(dt,UmTheta_U2);
    
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, Um2Theta);
    
    Planeta.Mueva_r(dt,UmTheta_U2);
    
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, Theta);
    
    Planeta.Mueva_r(dt,ThetaU2);
     
   
    
  }

  return 0;
}
