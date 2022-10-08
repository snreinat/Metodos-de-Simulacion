#include <iostream>
#include <cmath>
using namespace std;


const double ErrMax=1e-7;

double f(double x){
  return sin(x)/x;
}

double CerosPorBiseccion(double a, double b){
  double m, fa, fm;
  
  fa=f(a);
  while(b-a>ErrMax){
    //Reduzca el intervalo a la mitad
    m=(a+b)/2; fm=f(m);
    if(fa*fm>0) //Si son del mismo signo
      {a=m;  fa=fm;} //Corro a hasta m
    else b=m; //Corro b hasta m
  }
  return (a+b)/2;
}



int main(){
  
  double a=2, b=4;
  cout<<"El cero es "<<CerosPorBiseccion(a,b)<<endl;
  
  return 0;
}
