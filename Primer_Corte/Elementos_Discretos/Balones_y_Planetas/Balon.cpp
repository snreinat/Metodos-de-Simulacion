#include <iostream>
#include <cmath>
using namespace std;

//constantes Globales
const double g=9.8;

//Declaración de las clases
class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  double Getx(void){return x;}; //Función Inline
  double Gety(void){return y;};
  
};

void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
  Fx=0; Fy=-m*g;
}

void Cuerpo::Muevase(double dt){
  x+=Vx*dt; y+=Vy*dt;
  Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}

//------------Funciones Globales---------
int main(){
  Cuerpo Balon;
  double t, dt=0.1;
  
  //..........(x0,y0,Vx0,Vy0,m0,R0)
  Balon.Inicie(0,0,4,3,0.453,0.15);
  
  for(t=0; t<3; t+=dt){
    cout<<Balon.Getx()<<"\t"<<Balon.Gety()<<endl;
    Balon.CalculeFuerza();
    Balon.Muevase(dt);
  }
  return 0;
}
