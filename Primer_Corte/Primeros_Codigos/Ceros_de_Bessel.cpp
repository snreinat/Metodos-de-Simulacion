#include <iostream>
#include <cmath>
using namespace std;

const double ErrMax=1e-7;

double f(double alpha, double x, double t){
  return cos(alpha*t-x*sin(t));
}

//reemplazamos la variable de integración por t
double IntegralPorSimpson(double alpha, double x, double a, double b, int n){//se añaden las variables que falten de la función
  double t,h,suma;
  int i;
  n*=2;
  h=(b-a)/n;
  suma=0;
  
  for(i=0; i<=n; i++){
    t=a+i*h;
    
    if(i==0 || i==n)//Si es el primero o el último
      suma+=f(alpha,x,t);
    
    else if(i%2==0)//Si es par
      suma+=2*f(alpha,x,t);
    
    else
      suma+=4*f(alpha,x,t);//Si es impar
    
	}
  return suma*h/3;
  
  
}

//Devuelve el valor de la integral para un alpha y un x
double Bessel(double alpha, double x){
  double a=0, b=M_PI; //limites de integración
  int n=50;
  return IntegralPorSimpson(alpha,x,a,b,n)/M_PI;
}

//Calcular los ceros de la función Bessel
double CerosPorBiseccion(double alpha, double a, double b){
  double m, fa, fm;
  
  fa=Bessel(alpha,a);
  while(b-a>ErrMax){
    //Reduzca el intervalo a la mitad
    m=(a+b)/2; fm=Bessel(alpha,m);
    if(fa*fm>0) //Si son del mismo signo
      {a=m;  fa=fm;} //Corro a hasta m
    else b=m; //Corro b hasta m
  }
  return (a+b)/2;
}

int main(){
  
  double alpha=0;
  double a=2, b=4;
  
  cout<<"La función de Bessel con alpha="<<alpha<<" tiene un 0 en x="<<CerosPorBiseccion(alpha,a,b)<<endl;
  
  return 0;
}
