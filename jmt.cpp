#include <iostream>
#include <functional>
#include <cmath>
#include <fstream>
#include <gsl/gsl_ieee_utils.h>

const double jmax = 10;
const double amin = -10;
const double amax = 10;

double secant(const std::function<double(double)>& f,double min,double max,double tol,const char* call_name) {
  double x0(min),x1(max);
  double f0(f(x0)),f1(f(x1));
  while(true) {
    //std::cout<<" x0 : "<<x0<<" f0 : "<<f0<<" x1 : "<<x1<<" f1 : "<<f1<<std::endl;
    double xs = (f0*x1-f1*x0)/(f0-f1);
    double fs = f(xs);
    if(fabs(fs)<tol) {
      std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return xs "<<std::endl;
      return xs;
    } else {
      if(f0*fs<0) {
        f1 = fs;
        x1 = xs;
      } else if (f1*fs<0) {
        f0 = fs;
        x0 = xs;
      } else {
        if(fabs(f0)<fabs(f1)) {
	  std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return x0 "<<std::endl;
	  return x0;
        } else {
	  std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return x1 "<<std::endl;
          return x1;
	}
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

// acceleration and velocity are initial values
double dist(double v,double a,double j,double t) {
  return v*t+a*t*t*0.5+j*t*t*t/6.0;
}

double delta_v(double a,double j,double t) {
  return a*t+j*t*t*0.5;
}



double maxAcc(double a0,double dv) {
  auto f = [](double a0,double dv,double extreme_jerk){
             return [a0,extreme_jerk,dv]  (double t1) {
                      double t2 = 2*t1+ a0/extreme_jerk ;// t1 -(a0-j*t1)/j;
                      double delta_v = deltaVelocity(a0,extreme_jerk,t1,t2);
                      return (delta_v-dv);
                    };
           };
  double t;
  auto fpj = f(a0,dv,jmax);
  auto fnj = f(a0,dv,-jmax);
  static char buffer[1000];
  sprintf(buffer,"maxAcc-fpj a0 %f dv %f",a0,dv);
  double t_fpj = secant(fpj,0,4,1e-9,buffer);
  sprintf(buffer,"maxAcc-fnj a0 %f dv %f",a0,dv);
  double t_fnj = secant(fnj,0,4,1e-9,buffer);
  double a_extreme;
  if(fabs(fpj(t_fpj))<fabs(fnj(t_fnj))) {
    a_extreme = a0+jmax*t_fpj;
  } else {
    a_extreme = a0-jmax*t_fnj;
  }
  if(a_extreme<=amin)
    a_extreme = amin;
  else if(a_extreme>=amax)
    a_extreme = amax;
  return a_extreme;
}

std::tuple<double,std::function<double(double)>,double> distDuringVelocityChange(double a0,double v0,double v1) {
  double a_extreme = maxAcc(a0,v1-v0);
  double t1 = fabs(a_extreme-a0)/jmax;
  double t3 = fabs(a_extreme)/jmax;
  double j1,j3;
  if(a_extreme>a0) {
    j1 = jmax;
    j3 = -jmax;
  } else {
    j1 = -jmax;
    j3 = jmax;
  }
  double u1 = v0+delta_v(a0,j1,t1);
  double u2 = v1-delta_v(a_extreme,j3,t3);
  double t2 = (u2-u1)/a_extreme;
  if(t2<0)
    std::cout<<"ERROR"<<std::endl;
  auto jfn=[t1,t2,t3,j1,j3](double t) {
    if(t<0)
      std::cout<<" argument to jerk function out-of-bound "<<t<<std::endl;
    else if(t<t1) {
      return j1;
    } else if (t<t1+t2) {
      return 0.0;
    } else if (t<t1+t2+t3) {
      return j3;
    }
  };
  return std::make_tuple(dist(v0,a0,j1,t1)+dist(u1,a_extreme,0,t2)+dist(u2,a_extreme,j3,t3),jfn,t1+t2+t3);
}

// final velocity is assumed to be zero
double peakVel(double a0,double v0,double vmin,double vmax,double delta_d) {
  // assuming there is no max or min velocity
  auto f = [a0,v0,delta_d](double vpeak) {
    return std::get<0>(distDuringVelocityChange(a0,v0,vpeak))+std::get<0>(distDuringVelocityChange(0,vpeak,0))-delta_d;
  };
  static char buffer[1000];
  sprintf(buffer,"peakVel a0 %f v0 %f delta_d %f",a0,v0,delta_d);
  double vpeak = secant(f,-200,200,1e-5,buffer);
  if(vpeak>vmax)
    vpeak = vmax;
  else if(vpeak<vmin)
    vpeak = vmin;
  return vpeak;  
}

std::tuple<std::function<double(double)>,double> achieveZeroAcelAndVel(double a0,double v0,double vmin,double vmax,double delta_d) {
  double v_extreme = peakVel(a0,v0,vmin,vmax,delta_d);
  auto vChangeData1 = distDuringVelocityChange(a0,v0,v_extreme);
  auto vChangeData3 = distDuringVelocityChange(0,v_extreme,0);
  double t3 = std::get<2>(vChangeData3);
  double t1 = std::get<2>(vChangeData1);
  double d3 = std::get<0>(vChangeData3);
  double d1 = std::get<0>(vChangeData1);
  double d2 = delta_d-d1-d3;
  double t2 = d2/v_extreme;
  auto jfn3 = std::get<1>(vChangeData3);
  auto jfn1 = std::get<1>(vChangeData1);
  std::cout<<" a0 : "<<a0<<" v0 : "<<v0<<" jmax : "<<jmax <<" vmin : "<<vmin<<" vmax : "<<vmax<<" delta_d : "<<delta_d<<std::endl;
  std::cout<<"achieveZeroAcelAndVel :  v_extreme : "<<v_extreme<<" t1 : "<<t1<<" t2 : "<<t2<<" t3 : "<<t3
	   <<" d1 : "<<d1<<" d2 : "<<d2<<" d3 : "<<d3<<std::endl;
  auto jfn = [t1,t2,t3,jfn1,jfn3](double t) {
    if(t>=0 && t < t1) {
      return jfn1(t);
    } else if (t<=t1+t2) {
      return 0.0;
    } else if (t<=t1+t2+t3) {
      return jfn3(t-t1-t2);
    }
  };
  double total_t = t1+t2+t3;
  return std::make_tuple(jfn,total_t);
}



void path(double a0,double v0,double vmin,double vmax,double delta_d,const char* fname) {
  auto jfn_t = achieveZeroAcelAndVel(a0,v0,vmin,vmax,delta_d);
  auto jfn = std::get<0>(jfn_t);
  auto total_t = std::get<1>(jfn_t);
  std::cout<<"total_t : "<<total_t<<std::endl;
  double ap(a0),vp(v0),sp(0);
  double dt = 0.02;
  std::ofstream fout(fname);
  fout<<"t,j,a,v,s"<<std::endl;
  fout<<"0,0,"<<a0<<","<<v0<<",0"<<std::endl;
  for(double t = dt; t<=total_t;t+=dt) {
    double ac,vc,sc,jc;
    jc = jfn(t);
    ac = ap+jc*dt;
    vc = vp+(ap+ac)*dt*0.5;
    sc = sp+(vp+vc)*dt*0.5;
    fout<<t<<","<<jc<<","<<ac<<","<<vc<<","<<sc<<std::endl;
    ap = ac;
    vp = vc;
    sp = sc;
  }
}

int main() {

  gsl_ieee_env_setup(); /* read GSL_IEEE_MODE */
  double a0(0),dv(5),tol(1e-7);
  //void path(double a0,double v0,double vmin,double vmax,double delta_d,const char* fname) 
  path(0.0,0.0, 0,23, 100,"trial1.csv");
  path(5.0,0.0,0,23,100,"trial2.csv");
  path(10.0,2.0,0,23,100,"trial3.csv");
  path(10.0,10.0,0,23,100,"trial4.csv");
  path(0.0,23.0,0,23,50,"trial5.csv");
}
