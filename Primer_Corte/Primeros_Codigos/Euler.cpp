#include <iostream>
#include <cmath>

using namespace std;

double f (double t,double x){
  return x;
}

double UnPasoDeEuler(double &t, double &x, double dt){
  double dx;
  dx=f(t,x)*dt;
  t+=dt;
  x+=dx;
}

int main(){
  double t, x, dt=0.001;

  for(t=0,x=1; t<2+dt/2; ){
    cout<<t<<"\t"<<x<<"\t"<<exp(t)<<endl;
    UnPasoDeEuler(t,x,dt);
  }

  return 0;
}
