#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;


const int Lx=256, Ly=256; //Arreglo bidimensional
const double p0=0.25, p=0.25;//p_0 es la probabilidad de quedarse quieto y p la probabilidad de girar a la derecha 90 grados

const int Q=4;//Número de flechas, en este caso tenemos 4


//----------Clase LatticeGas------------
class LatticeGas{
private:
  int Vx[Q], Vy[Q]; //Vectores de cooredenadas de las direcciones 
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q];//n[ix][i]
public:
  LatticeGas(void);
  void Inicie(int N, double mu, double sigma);
  void Show(void);
  void Shownew(void);
  void Colisione(void);
  void Adveccione(void);
  double rho(int ix, int iy); 
  double Varianza(void);
  void GrafiqueRho(void);
};


LatticeGas::LatticeGas(void){
  
  Vx[0]=1; Vx[1]=0; Vx[2]=-1; Vx[3]=0; 
  Vy[0]=0; Vy[1]=1; Vy[2]=0;  Vy[3]=-1;

}


void LatticeGas::Inicie(int N, double mu, double sigma){
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      double rho=(N/(sigma*sigma*2.0*M_PI))*(exp(-0.5*pow((ix-mu)/sigma,2.0)-0.5*pow((iy-mu)/sigma,2.0)));// rho es la densidad de probabilidad Gaussiana
    for(int i=0;i<Q;i++)
      f[ix][iy][i]=rho/Q;
    }
  }
  }

/*
void LatticeGas::Inicie(int N, double mu, double sigma){
  for(int ix=0;ix<Lx;ix++)
     for(int iy=0;iy<Ly;iy++)
       for(int i=0;i<Q;i++)
	 f[ix][iy][i]=N/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2)-0.5*pow((iy-mu)/sigma,2));
	 }*/


void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){ // El for va al revés porque primero imprimo todos los que van a la derecha y después todos los que van a la izquierda
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	      cout<<f[ix][iy][i];
        //cout<< endl;
  }
    }
    cout<<endl;
  }
}

void LatticeGas::Shownew(void){
  for(int i=0;i<Q;i++){ // El for va al revés porque primero imprimo todos los que van a la derecha y después todos los que van a la izquierda
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	      cout<<fnew[ix][iy][i];
        //cout<< endl;
  }
    }
    cout<<endl;
  }
}


void LatticeGas::Colisione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
	fnew[ix][iy][i]=p0*f[ix][iy][i]+p*f[ix][iy][(i+1)%4]+(1-p0-2*p)*f[ix][iy][(i+2)%4]+p*f[ix][iy][(i+3)%4];
}

//Moverse a las siguientes celdas
void LatticeGas::Adveccione(void){
   for (int ix=0; ix<Lx;ix++){
     for (int iy=0; iy<Ly; iy++){
       for(int i=0;i<Q;i++){
	 f[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Ly][i]=fnew[ix][iy][i];//Condiciones de frontera periodicas 
      }
     }
   }
}

double LatticeGas::rho(int ix, int iy){
  double suma; int i;
  for(suma=0, i=0; i<Q; i++)
    suma+=f[ix][iy][i];
  return suma;
}


double LatticeGas::Varianza(void){
  int ix, iy; double N, R, Rprom, Sigma2;
  
  //Calcular N
  for (N=0,ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      N+=rho(ix,iy);
    }
  }

  //Calcular Xprom
  for(Rprom=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      R=0;
      // R=sqrt(ix*ix+iy*iy);
      R=std::hypot(ix,iy);
      Rprom+=R*rho(ix,iy);
    }
  }
  Rprom=Rprom/N;
  
  //Calcular Sigma2 para X
  for(Sigma2=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      R=0;
      // R=sqrt(ix*ix+iy*iy);
       R=std::hypot(ix,iy);
      Sigma2+=pow((R-Rprom),2.0)*rho(ix,iy);
    }
  }
  Sigma2=Sigma2/(N-1);
  
  //Calcular Sigma2 para y
   return Sigma2;
  
}

void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      cout<<ix<<iy<<" "<<rho(ix,iy)<<endl;
      }


//-----------Programa Principal------------

int main(void){
  
  LatticeGas Difusion;
  int N=2400;
  double mu=Lx/2, sigma=16;
  int t, tmax=350;

 
  Difusion.Inicie(N, mu, sigma);

 for(t=0; t<tmax;t++){
   cout<<t<<" "<<Difusion.Varianza()<<endl;
   Difusion.Colisione();
   Difusion.Adveccione();
   //Difusion.Shownew();
   
   //cout<<endl;

  }
 // Difusion.GrafiqueRho();
 // Difusion.Show();

   
 return 0;
}
