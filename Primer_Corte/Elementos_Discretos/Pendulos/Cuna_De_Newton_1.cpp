// Simular el movimiento de 2 planetas por PEFRL
#include <iostream>
#include <cmath>
#include "vector.h" 
using namespace std;


//------------------------------Declarar constantes------------------

const int N=1;  //numero de cuerpos
const double g=980;
const double K=1e013;

const double E=0.1786178958448091e00;
const double L=-0.2123418310626054e0;
const double X=-0.6626458266981849e-1;
const double coeficiente1=(1-2*L)/2;
const double coeficiente2=(1-2*(X+E))*2/2;

//------------------------------Declarar clases----------------------
class Cuerpo;
class Colisionador;

//-----------------------------Implementar clases--------------------

//------------------------------Clase Cuerpo-------------------------
class Cuerpo{
private:
  double  Theta, Omega, Tau;   double m, R, I, l, x0;
public:
  void Inicie(double Theta0, double Omega0, double m0,double R0, double L0, double x00);
  void BorreTorque(void){Tau=0;};
  void AdicioneTorque(double Tau0){Tau+=Tau0;};
  void Mueva_Theta(double dt, double coeficiente);
  void Mueva_Omega(double dt, double coeficiente);
  void Dibujese(void);
  double GetTau(void){return Tau;};  //inline
  double GetX(void){return x0+l*std::sin(Theta);};
  double GetY(void){return -l*std::cos(Theta);};
  double GetTheta(void){return Theta;}; 

  

  friend class Colisionador;
};
void Cuerpo::Inicie(double Theta0, double Omega0, double m0,double R0, double L0,double x00){
  Theta=Theta0; Omega=Omega0; m=m0;  R=R0; l=L0; I=m*l*l; x0=x00;
}

void Cuerpo::Mueva_Theta(double dt, double coeficiente){
  Theta+=Omega*dt*coeficiente;
}

void Cuerpo::Mueva_Omega(double dt, double coeficiente){
  Omega+=(Tau*dt*coeficiente)/I;
}

void Cuerpo::Dibujese(void){
  cout<<" , "<< GetX() <<"+"<<R<<"*cos(t),"<< GetY()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<Theta<<"),-"<<l/7<<"*t*cos("<<Theta<<")";
}

//---------------------------------Clase Colisionador-----------------------

class Colisionador{

private:

public:
  void CalculeTorques(Cuerpo * Pendulo);
  void CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2);
};

void Colisionador::CalculeTorques(Cuerpo * Pendulo){

  int i;
  //---------------------------borrar todas las fuerzas----------------------------- 
  for(i=0; i<N; i++){
  Pendulo[i].BorreTorque();
  }

  //---------------------------fuerzas individuales----------------------------------
  
  for(i=0; i<N; i++) {
    Pendulo[i].AdicioneTorque(-Pendulo[i].l*Pendulo[i].m*g*std::sin(Pendulo[i].Theta));
  }
  //---------------------------Calcular todas las fuerzas de colision----------------
		       

  for (i=0; i<N-1; i++){
   CalculeTorqueEntre(Pendulo[i], Pendulo[i+1]);
   }

}
void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2){

  double s=(Pendulo1.GetX()+Pendulo1.R)-(Pendulo2.GetX()-Pendulo2.R);

  if (s>0){
    double F=K*std::pow(s,1.5);
    double T2=F*Pendulo2.l;
    
  Pendulo1.AdicioneTorque(-T2); Pendulo2.AdicioneTorque(T2);
  }

  
}


//-------------------------- Funciones de Animacion -------------------

void InicieAnimacion(void){
  // cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'CunaDeNewton.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-14:14]"<<endl;
  cout<<"set yrange[-18:0]"<<endl;
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

//-------------------------------Programa Principal-----------------------------------------------------

int main(void){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  int i;
  double m0=100, R0=2, L0=12;
  double T=2*M_PI*std::sqrt(L0/g);
  
  
  double t, tdibujo, tmax=3*T, tcuadro=T/75, dt=0.0000001;

 InicieAnimacion(); //Dibujar

//-------------------------------Iniciar los pendulos--------------------------------------------------------
 
//-------------------------------(theta0, omega0, m0,R0, L0, x00)
                    Pendulo[0].Inicie( -0.6, 0, m0,R0, L0, 0);
 for(i=1; i<N; i++) Pendulo[i].Inicie(    0, 0, m0,R0, L0,R0*(2*i));
 
 for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
   
 //-------------------------------Dibujar animacion----------------------------------------------------
    
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      
      for(int i=0; i<N; i++)
	Pendulo[i].Dibujese();
      
      TermineCuadro();
      
      // hacer un plot
      //std::cout<< Pendulo[0].Getx()<<"\t"<<Pendulo[0].Gety()<<"\t"<< Pendulo[1].Getx()<<"\t"<<Pendulo[1].Gety()<<endl;
      tdibujo=0;
    }
 //-----------------------------------muevase por OMELYAN PEFRL----------------------------------------------
    
    for(i=0; i<N; i++)Pendulo[i].Mueva_Theta(dt,E);
    Newton.CalculeTorques(Pendulo);  for(i=0; i<N; i++)Pendulo[i].Mueva_Omega(dt,coeficiente1);
    for(i=0; i<N; i++)Pendulo[i].Mueva_Theta(dt,X);
    Newton.CalculeTorques(Pendulo);  for(i=0; i<N; i++)Pendulo[i].Mueva_Omega(dt,L);
    for(i=0; i<N; i++)Pendulo[i].Mueva_Theta(dt,coeficiente2);
    Newton.CalculeTorques(Pendulo);  for(i=0; i<N; i++)Pendulo[i].Mueva_Omega(dt,L);
    for(i=0; i<N; i++)Pendulo[i].Mueva_Theta(dt,X);
    Newton.CalculeTorques(Pendulo);  for(i=0; i<N; i++)Pendulo[i].Mueva_Omega(dt,coeficiente1);
    for(i=0; i<N; i++)Pendulo[i].Mueva_Theta(dt,E);
    
  }   
  return 0;
}
