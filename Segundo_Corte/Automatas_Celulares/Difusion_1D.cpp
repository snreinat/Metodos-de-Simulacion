#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=1024;
const double p=0.5;//Probabilidad de quedarse quieto

const int Q=2;//Número de flechas


//----------Clase LatticeGas------------
class LatticeGas{
private:
  int V[Q]; //Matriz de Q direcciones.V[i] i=0 (derecha) i=1 (izquierda)
  int n[Lx][Q], nnew[Lx][Q];//n[ix][i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double mu, double sigma, Crandom &ran64);//Importante pasar por referencia para solo usar un generador aleatorio
  void Show();
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double rho(int ix){return n[ix][0]+n[ix][1];};
  double Varianza(void);
  void GrafiqueRho(void);
};

LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;  
}

double LatticeGas::Varianza(void){
  int ix; double N, Xprom, Sigma2;
  
  //Calcular N
  for (N=0, ix=0; ix<Lx; ix++)
    N+=rho(ix);
  
  //Calcular Xprom
  for(Xprom=0, ix=0; ix<Lx; ix++)
    Xprom+=ix*rho(ix);
  Xprom/=N;
  
  //Calcular Sigma2
   for(Sigma2=0, ix=0; ix<Lx; ix++)
     Sigma2+=pow(ix*-Xprom,2.0)*rho(ix);
   Sigma2/=(N-1);

   return Sigma2;
  
}
void LatticeGas::Borrese(void){
  for (int ix=0; ix<Lx;ix++)
    for(int i=0;i<Q;i++)
      n[ix][i]=nnew[ix][i]=0;
}

void LatticeGas::Inicie(int N, double mu, double sigma, Crandom &ran64){
  int ix,i;
  while (N>0){ 
    ix=(int)ran64.gauss(mu,sigma);//Escoger una celda al azar, el (int) lo convierte a entero
    if(ix<0) ix=0; if(ix>Lx-1) ix=Lx-1;//COrrgir en los bordes si es necesario
    i=(int)Q*ran64.r(); //Escoger una dirección al azar, elige un número entre 0 y Q
    
    if (n[ix][i]==0) //si (ese sitio está)
      { n[ix][i]=1; N--;}
   //pongo una bolita ahí y decrezco N
    // Busco otro sitio 
	 }
}

void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){ // El for va al revés porque primero imprimo todos los que van a la derecha y después todos los que van a la izquierda
    for (int ix=0; ix<Lx;ix++)
      cout<<n[ix][i];
    cout<< endl;
  }
}

void LatticeGas::GrafiqueRho(void){
  //for(int ix=0;ix<Lx;ix++)
  // cout<<ix<<" "<<rho(ix)<<endl;
}

void LatticeGas::Colisione(Crandom & ran64){

  for (int ix=0; ix<Lx;ix++){ //Para cada celda
    if(ran64.r()>p){//Genero un número al azar y si es mayor que p volteo
      nnew[ix][0]=n[ix][1]; nnew[ix][1]==n[ix][0];//intercambia los contenidos
    }
      else
	{ nnew[ix][0]=n[ix][0]; nnew[ix][1]==n[ix][1];}// No los intercambia
  }
}
void LatticeGas::Adveccione(void){
   for (int ix=0; ix<Lx;ix++)
	for(int i=0;i<Q;i++)
	  n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i];// el %Lx genera condiciones de frontera periodicas 
  }


//-----------Programa Principal------------

int main(void){
  
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=400;
  double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;

  Difusion.Borrese();
  Difusion.Inicie(N, mu, sigma, ran64);

  for(t=0; t<tmax;t++){
    cout<<t<<" "<<Difusion.Varianza()<<endl;
      Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  
  return 0;
}
