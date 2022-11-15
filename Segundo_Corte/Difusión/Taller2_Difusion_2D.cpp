#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=256, Ly=256; //Arreglo bidimensional
const double p_0=0.25, p=0.25;//p_0 es la probabilidad de quedarse quieto y p la probabilidad de girar a la derecha 90 grados

const int Q=4;//Número de flechas, en este caso tenemos 4


//----------Clase LatticeGas------------
class LatticeGas{
private:
  int Vx[Q], Vy[Q]; //Vectores de cooredenadas de las direcciones 
  int n[Lx][Ly][Q], nnew[Lx][Ly][Q];//n[ix][i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double mu, double sigma, Crandom &ran64, Crandom &ran264);
  void Show(void);
  void Shownew(void);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double rho(int ix, int iy){ return (n[ix][iy][0]+n[ix][iy][1]+n[ix][iy][2]+n[ix][iy][3]);};
  //double rhox(int ix, int iy);
  //double rhoy(int ix, int iy);
  double Varianza(void);
  void GrafiqueRho(void);
};


LatticeGas::LatticeGas(void){
  
  Vx[0]=1; Vx[1]=0; Vx[2]=-1; Vx[3]=0; 
  Vy[0]=0; Vy[1]=1; Vy[2]=0;  Vy[3]=-1;

}

void LatticeGas::Borrese(void){
  for (int ix=0; ix<Lx;ix++){
    for(int iy=0; iy<Ly; iy++){
      for(int i=0;i<Q;i++){
	n[ix][iy][i]=nnew[ix][iy][i]=0;
      }
    }
  }
}

void LatticeGas::Inicie(int N, double mu, double sigma, Crandom &ran64, Crandom &ran264){
  int ix,iy,i;
  
  while (N>0){ 
    ix=(int) ran64.gauss(mu,sigma);//Escoger una celda al azar
    iy=(int) ran264.gauss(mu,sigma);
    if(ix<0) ix=0; if(ix>Lx-1) ix=Lx-1;//Corregir en los bordes si es necesario
    if(iy<0) iy=0; if(iy>Ly-1) iy=Ly-1;
    i=(int) Q*ran64.r(); //Escoger una dirección al azar, elige un número entre 0 y Q
    
    if (n[ix][iy][i]==0) //si ese sitio está vacío
      { n[ix][iy][i]=1; N--;}//pongo una bolita ahí y decrezco N
    	 }
}


void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){ // El for va al revés porque primero imprimo todos los que van a la derecha y después todos los que van a la izquierda
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	      cout<<n[ix][iy][i];
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
	      cout<<nnew[ix][iy][i];
        //cout<< endl;
  }
    }
    cout<<endl;
  }
}




void LatticeGas::Colisione(Crandom & ran64){
  
  for (int ix=0; ix<Lx; ix++){ 
    for(int iy=0; iy<Ly; iy++){
double PB=ran64.r();//Generar un número al azar entre 0 y 1
      if(PB>=0 && PB<=p_0){
	nnew[ix][iy][0]=n[ix][iy][0];//Girar 0°
	nnew[ix][iy][1]=n[ix][iy][1];
	nnew[ix][iy][2]=n[ix][iy][2];
	nnew[ix][iy][3]=n[ix][iy][3];
	}
      
      if(PB>p_0 && PB<=p_0+p){
	
	nnew[ix][iy][0]=n[ix][iy][3];//Gira 90° a la derecha
	nnew[ix][iy][1]=n[ix][iy][0];
	nnew[ix][iy][2]=n[ix][iy][1];
	nnew[ix][iy][3]=n[ix][iy][2];
      }
      
      if(PB>p_0+p && PB<=p_0+2*p){
	
	nnew[ix][iy][0]=n[ix][iy][1];//Gira 90° a la izquierda	
  nnew[ix][iy][1]=n[ix][iy][2];
	nnew[ix][iy][2]=n[ix][iy][3];
	nnew[ix][iy][3]=n[ix][iy][0];
      }
      
    
      if(PB>(p_0+2*p) && PB<=1){
	nnew[ix][iy][0]=n[ix][iy][2];// Gira 180°
	nnew[ix][iy][1]=n[ix][iy][3];
	nnew[ix][iy][2]=n[ix][iy][0];
	nnew[ix][iy][3]=n[ix][iy][1];
      }         
    }
  }
}

//Moverse a las siguientes celdas
void LatticeGas::Adveccione(void){
   for (int ix=0; ix<Lx;ix++){
     for (int iy=0; iy<Ly; iy++){
       for(int i=0;i<Q;i++){
	 n[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Ly][i]=nnew[ix][iy][i];//Condiciones de frontera periodicas 
      }
     }
   }
}


/*double LatticeGas::rhox(int ix, int iy){
  return n[ix][iy][0]+n[ix][iy][2];
}

double LatticeGas::rhoy(int ix, int iy){
  return n[ix][iy][1]+n[ix][iy][3];
}

*/

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
      R=sqrt(ix*ix+iy*iy);
      Rprom+=R*rho(ix,iy);
    }
  }
  Rprom=Rprom/N;
  
  //Calcular Sigma2 para X
  for(Sigma2=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      R=0;
      R=sqrt(ix*ix+iy*iy);
      Sigma2+=pow((R-Rprom),2.0)*rho(ix,iy);
    }
  }
  Sigma2=Sigma2/(N-1);
  
  //Calcular Sigma2 para y
   return Sigma2;
  
}
/*
double LatticeGas::Varianza(void){
  int ix, iy; double N, Nx, Ny, Xprom, Yprom, Sigma2, Sigma2X, Sigma2Y, Cov;
  
  //Calcular N
  for (Nx=0,Ny=0,ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++){
      Nx+=rhox(ix,iy);
      Ny+=rhoy(ix,iy);
    }
  N=Nx+Ny;
  
  //Calcular Xprom
  for(Xprom=0, Yprom=0, ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++){
      Xprom+=ix*rhox(ix,iy);
      Yprom+=iy*rhoy(ix,iy);
    }
  Xprom/=Nx;
  Yprom/=Ny;
    
  
  //Calcular Sigma2 para X
  for(Sigma2X=0,Sigma2Y=0,Cov, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      Sigma2X+=pow((ix-Xprom),2.0)*rhox(ix,iy);//rhox(ix,iy);
      Sigma2Y+=pow((iy-Yprom),2.0)*rhoy(ix,iy);//rhoy(ix,iy);
      Cov+=(ix-Xprom)*(iy-Yprom);

    }
  }
  
  //Sigma2X/=(Nx-1);
  //Sigma2Y/=(Ny-1);
  Sigma2=((Sigma2X+Sigma2Y)-2*Cov)/(N-1);
 
  
  //Calcular Sigma2 para y
  return Sigma2;
  
}
  
*/

  /*
void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      cout<<ix<<" "<<rhox(ix,iy)<<endl;
      }*/





//-----------Programa Principal------------

int main(void){
  
  LatticeGas Difusion;
  Crandom ran64(2);//Se puede iniciar con cualquier semilla que no sea 0
  Crandom ran264(3);
  int N=2400;
  double mu=Lx/2, sigma=16;
  int t, tmax=350;

  Difusion.Borrese();
  Difusion.Inicie(N, mu, sigma, ran64, ran264);

 for(t=0; t<tmax;t++){
   cout<<t<<" "<<Difusion.Varianza()<<endl;
   Difusion.Colisione(ran64);
   Difusion.Adveccione();
   //Difusion.Show();
   //Difusion.Shownew();
   
   //cout<<endl;

  }
 //Difusion.GrafiqueRho();
   
 return 0;
}
