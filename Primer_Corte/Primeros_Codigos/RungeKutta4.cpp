#include <iostream>
#include <cmath>

using namespace std;

double f (double t,double x){
  return x;
}

double UnPasoDeRungeKutta4(double &t0, double &x0, double dt){
  double dx1,dx2,dx3,dx4;
  dx1=f(t0,x0)*dt;
  dx2=f(t0+dt/2,x0+dx1/2)*dt;
  dx3=f(t0+dt/2,x0+dx2/2)*dt;
  dx4=f(t0+dt,x0+dx3)*dt;
  t0+=dt;
  x0+=(dx1+2*(dx2+dx3)+dx4)/6;
}

int main(){
  double t, x, dt=0.001;

  for(t=0,x=1; t<2+dt/2; ){
    cout<<t<<"\t"<<x<<"\t"<<exp(t)<<endl;
    UnPasoDeRungeKutta4(t,x,dt);
  }

  return 0;
}
