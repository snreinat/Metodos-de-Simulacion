#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes globales
const double GM=1.0;

//Declaración de las clases
class Cuerpo;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  vector3D r,rold,V,F;
  double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0, double Vy0,double Vz0,double m0,double R0);
  void CalculeFuerza(void);
  void Arranque(double dt);
  void Muevase(double dt);
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
void Cuerpo::Arranque(double dt){
  rold=r-V*dt+F*(dt*dt/(2*m));
}
void Cuerpo::Muevase(double dt){
  vector3D rnew;
  rnew=2*r-rold+(dt*dt*F)/m;
    V=(r-rold)/(2*dt);
    rold=r; r=rnew;
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
  Planeta.CalculeFuerza();
  Planeta.Arranque(dt);
  
  for(t=0;t<1.1*T;t+=dt){
    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }

  return 0;
}
