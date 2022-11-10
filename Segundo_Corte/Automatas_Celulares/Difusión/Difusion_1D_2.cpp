#include <iostream>
#include <cmath>
using namespace std;

const int Lx=1024;
const double p=0.5;//Probabilidad de quedarse quieto

const int Q=2;//Número de flechas


//----------Clase LatticeGas------------
class LatticeGas{
private:
  int V[Q]; //Matriz de Q direcciones.V[i] i=0 (derecha) i=1 (izquierda)
  int f[Lx][Q], fnew[Lx][Q];//f[ix][i]
public:
  LatticeGas(void);
  void Inicie(int N, double mu, double sigma);//Importante pasar por referencia para solo usar un generador aleatorio
  void Colisione(void);
  void Adveccione(void);
  double rho(int ix);
  double Varianza(void);
  void GrafiqueRho(void);
};

LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;  
}

void LatticeGas::Inicie(int N, double mu, double sigma){
  for(int ix=0;ix<Lx;ix++){
    double rho=(N/(sigma*sqrt(2.0*M_PI)))*exp(-0.5*pow((ix-mu)/sigma,2.0));
  for(int i=0;i<Q;i++)
    f[ix][i]=rho/Q;
  }
}


double LatticeGas::rho(int ix){
  double suma; int i;
  
  for(suma=0,i=0;i<Q;i++)
    suma+=f[ix][i];
  
  return suma;
}


void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
   cout<<ix<<" "<<rho(ix)<<endl;
}

void LatticeGas::Colisione(void){
  int ix,i,j;
 
  for (ix=0; ix<Lx;ix++)//Para cada celda
    for(i=0;i<Q;i++){
      j=(i+1)%Q;
      fnew[ix][i]=f[ix][i]+(1-p)*(f[ix][j]-f[ix][i]);
    }
}


void LatticeGas::Adveccione(void){
   for (int ix=0; ix<Lx;ix++)
	for(int i=0;i<Q;i++)
	  f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i];// el %Lx genera condiciones de frontera periodicas 
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
     Sigma2+=pow((ix-Xprom),2.0)*rho(ix);
   Sigma2/=(N-1);

   return Sigma2;
  
}

//-----------Programa Principal------------

int main(void){
  
  LatticeGas Difusion;
  int N=6;
  double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;


  Difusion.Inicie(N, mu, sigma);
  
  for(t=0; t<tmax;t++){
    cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }
  // Difusion.GrafiqueRho();
  
  return 0;
}
