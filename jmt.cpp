#include <iostream>
#include <functional>
#include <cmath>

double secant(const std::function<double(double)>& f,double min,double max,double tol) {
  double x0(min),x1(max);
  double f0(f(x0)),f1(f(x1));
  while(true) {
    //std::cout<<" x0 : "<<x0<<" f0 : "<<f0<<" x1 : "<<x1<<" f1 : "<<f1<<std::endl;
    double xs = (f0*x1-f1*x0)/(f0-f1);
    double fs = f(xs);
    if(fabs(fs)<tol)
      return xs;
    else {
      if(f0*fs<0) {
        f1 = fs;
        x1 = xs;
      } else if (f1*fs<0) {
        f0 = fs;
        x0 = xs;
      } else {
        if(fabs(f0)<fabs(f1))
          return x0;
        else
          return x1;
      }
    }
  }
}

double newton_raphson(const std::function<double(double)>& f,double x0,double tol) {
  double f0 = f(x0);
  double h = 1e-9;

  while(fabs(f0)>tol) {
    double fdash0 = (f(x0+h)-f0)/h;
    x0-=f0/fdash0;
    f0 = f(x0);
    //std::cout<<" x0 : "<<x0 <<" f0 : "<<f0<<" fdash0 : "<<fdash0<<std::endl;
  }
  return x0;
}

// assuming maximum jerk is constant
// assuming a negative jerk to begin with
// assuming there is no maximum acceleration
double deltaVelocity(double a0,double j,double t1,double t2) {
  double ret= a0*t1+j*t1*t1*0.5+(a0+j*t1)*(t2-t1)-j*(t2-t1)*(t2-t1)*0.5;
  return ret;
}


double maxAcc(double a0,double dv) {
  double j = 10;
  double amin=-10,amax=10;
  auto f = [](double a0,double dv,double extreme_jerk){
             return [a0,extreme_jerk,dv]  (double t1) {
                      double t2 = 2*t1+ a0/extreme_jerk ;// t1 -(a0-j*t1)/j;
                      double delta_v = deltaVelocity(a0,extreme_jerk,t1,t2);
                      return (delta_v-dv);
                    };
           };
  double t;
  auto fpj = f(a0,dv,j);
  auto fnj = f(a0,dv,-j);
  double t_fpj = secant(fpj,0,4,1e-9);
  double t_fnj = secant(fnj,0,4,1e-9);
  //double t_fpj_nr = newton_raphson(fpj,2.0,1e-9);
  //double t_fnj_nr = newton_raphson(fnj,2.0,1e-9);
  //std::cout<<" t_fpj_nr : "<<t_fpj_nr<<" t_fnj_nr : "<<t_fnj_nr<<std::endl;
  std::cout<<" t_fpj : "<<t_fpj<<" t_fnj : "<<t_fnj<<std::endl;
  std::cout<<" f(t_fpj) : "<<fpj(t_fpj)<<" t_fnj : "<<fnj(t_fnj)<<std::endl;
  double a_extreme;
  if(fabs(fpj(t_fpj))<fabs(fnj(t_fnj))) {
    a_extreme = a0+j*t_fpj;
  } else {
    a_extreme = a0-j*t_fnj;
  }
  if(a_extreme<=amin)
    a_extreme = amin;
  else if(a_extreme>=amax)
    a_extreme = amax;
  return a_extreme;
}


std::function<double(double)> polyfn(std::initializer_list<double> roots) {
  return [roots](double x) {
           double ret = 1.0;
           for(auto r : roots) {
             ret *= (x-r);
           }
           return ret;
         };
}

int main() {
  double a0(0),dv(5),j(10),amin(-10),amax(10),tol(1e-7);
  double apeak;
  for(double a0:{-10,-7,-5,-1,0,1,5,7,10})
    for(double dv:{20,15,10,5,1,0,-1,-5,-10,-20}) {
      apeak = maxAcc(a0,dv);
      std::cout<<" apeak : "<<apeak<<" a0 : "<<a0<<" dv : "<<dv<<std::endl;
    }
  /*
  for(double a0:{-10,-5,0,5,10})
    for(double dv:{20.,10.,5.,1.,0.1,0.01,0.001,-0.001,-0.01,-0.1,-1.,-5.,-10.,-20.})
      maxAcc(a0,dv);
  */
}
