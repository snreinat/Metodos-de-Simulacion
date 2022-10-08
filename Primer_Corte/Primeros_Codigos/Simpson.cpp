#include <iostream>
#include <cmath>
using namespace std;


double f(double x){
  return cos(x);
}

double IntegralPorSimpson(double a, double b, int n){
  double x,h,suma;
  int i;
  n*=2;
  h=(b-a)/n;
  suma=0;
  
  for(i=0; i<=n; i++){
    x=a+i*h;
    
    if(i==0 || i==n)//Si es el primero o el último
      suma+=f(x);
    
    else if(i%2==0)//Si es par
      suma+=2*f(x);
    
    else
      suma+=4*f(x);//Si es impar
    
	}
  return suma*h/3;
  
  
}

int main(){
  
  double a=0, b=M_PI/2;
  int n=50;
  
  cout<<"La integral es "<<IntegralPorSimpson(a,b,n)<<endl;
  
  return 0;
}
