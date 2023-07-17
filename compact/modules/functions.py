#!/usr/bin/env python
import sys
from math import log,sqrt,acos,exp

from ROOT import TMath


# asymetric gaussian
class GaussAss:
   def __call__( self, x, par ):

       out=0.0
       gaus1=0.0
       gaus2left=0.0
       gaus2right=0.0
       gaus1=par[0]*exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]))
       gaus2left=par[5]*exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[3]*par[3]))
       gaus2right=par[5]*exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[4]*par[4]))

       cc=par[2]*par[3]*par[4]

       if x[0]<=par[1] and cc>0:
                       out=gaus1+gaus2left;
       if x[0]>par[1] and cc>0:
                       out=gaus1+gaus2right;

       return out

# P1
class P1:
   def __call__( self, x, par ):
       out=par[0] + par[1]*x[0]; 
       return out;



# GaussAss+P1
class GaussAssBkg:
   def __call__( self, x, par ):

       out=0.0
       gaus1=0.0
       gaus2left=0.0
       gaus2right=0.0

       gaus1=par[0]*exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]))
       gaus2left=par[5]*exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[3]*par[3]))
       gaus2right=par[5]*exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[4]*par[4]))

       cc=par[2]*par[3]*par[4]

       if x[0]<=par[1] and cc>0:
                       out=gaus1+gaus2left;
       if x[0]>par[1] and cc>0:
                       out=gaus1+gaus2right;


       back = par[6] + par[7]*x[0];
       # add backg 
       backout=out+ back;

       return backout;

# threshold
class Bkg:
   def __call__( self, x, par ):
      xx=x[0]-1.435944
      f3g2=0.
      if par[1] < 0: return 1.e30;
      if (xx > 0.0 ):
           f3g2=par[0]* (xx**par[1])*(1+par[2]*xx)
           return f3g2;

class MGauss:
   def __call__( self, x, par ):
      f3g2=0.
      cost=1./math.sqrt(2*acos(-1.));
      f3g2=(par[0]*cost*exp(-((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2]) ));
      return f3g2;

# 3-parameter function for dijet analysis of 2015
class  ThreeParam2015:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x=x[0] / fCenterOfMassEnergy;
    ff=par[0]*TMath.Power((1.0-x),par[1])*TMath.Power(x,par[2]);
    return ff;

# 3-parameter function for dijet analysis of 2015
class  ThreeParam2015TEV:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x=x[0] / fCenterOfMassEnergy;
    ff=par[0]*TMath.Power((1.0-x),par[1])*TMath.Power(x,par[2]);
    return ff;

# 3-parameter function for dijet analysis of 2015 plus Gaussian 
class  ThreeParam2015Gauss:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x=x[0] / fCenterOfMassEnergy;
    ff=par[0]*TMath.Power((1.0-x),par[1])*TMath.Power(x,par[2]);
    ff1=0
    if (par[5]>0 and par[3]>0):
        ff1=par[3] * TMath.Gaus(x,par[4],par[5])
    return ff+ff1;



# 4-parameter function for dijet analysis of 2015
# FourParamFitFunction in src_util/SearchPhase.cxx
class  FourParam2015:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff=par[0]*TMath.Power((1.0-x),par[1])*TMath.Power(x,(par[2]+par[3]*log(x)));
    return ff;

class  FiveParam2015:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)))
    ff=ff1*ff2;
    return ff;

class  FiveParam2015TEV:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)))
    ff=ff1*ff2;
    return ff;

class  FiveParam2015TEVGauss:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)))
    ff=ff1*ff2;
    ff3=0;
    if (par[7]>0 and par[6]>0):
        ff3=par[5] * TMath.Gaus(x,par[6],par[7])
    return ff+ff3;

class  FiveParam2015Gauss:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)))
    ff=ff1*ff2;
    gaus1=0;
    if (par[7]>0): gaus1=par[5] * TMath.Gaus(x,par[6],par[7])
    return ff+gaus1;

class  FiveParam2015TEValt1:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]/sqrt(x) ))
    ff=ff1*ff2;
    return ff;


class  FourParam2015TEV:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)))
    ff=ff1*ff2;
    return ff;

class  SixParam2015TEV:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x) + par[5]*log(x)*log(x)*log(x)))
    ff=ff1*ff2;
    return ff;


class  FourParam2015TEV:
  def __call__( self, x, par ):
    Ecm=13.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)))
    ff=ff1*ff2;
    return ff;


class  FiveParam2015alt1:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]/sqrt(x) ))
    ff=ff1*ff2;
    return ff;


class  SixParam2015:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)+par[5]*log(x)*log(x)*log(x)))
    ff=ff1*ff2;
    return ff;

class  SixParam2015Gauss:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    ff1=par[0]*TMath.Power((1.0-x),par[1])
    ff2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)+par[5]*log(x)*log(x)*log(x)))
    ff=ff1*ff2;
    A=2.5066272 # sqrt(2*3.14159)
    gaus1=(1.0/(par[7]*A)) * par[6]*exp((-0.5)*(x-par[8])*(x-par[8])/(par[7]*par[7]))
    return ff+gaus1;


class FourParamUA2:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    f1=par[0]*TMath.Power(x,par[1])
    f2=TMath.Exp(-par[2]*x - par[3]*x*x)
    return f1*f2

class FourGammaGamma:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    RR=0.33333333333333333333
    f1=par[0]*TMath.Power(1-TMath.Power(x,RR),par[1])
    f2=TMath.Power(x,(par[2]+par[3]*log(x)+par[4]*log(x)*log(x)))
    return f1*f2
 


class GaussX:
  def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    x = x[0] / fCenterOfMassEnergy;
    gaus=par[0]*exp((-0.5)*(x-par[2])*(x-par[2])/(par[1]*par[1]))
    return gaus;

# par[1] - mean
# par[2] - sigma
class Gauss:
   def __call__( self, x, par ):
       out=par[0] * TMath.Gaus(x[0],par[1],par[2])
       return out;


# par[1] - mean
# par[2] - sigma
class GaussECM:
   def __call__( self, x, par ):
    Ecm=13000.;
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = False;
    x = x[0] / fCenterOfMassEnergy;
    out=par[0] * TMath.Gaus(x,par[1],par[2])
    return out;


class FGauss: 
  def __call__( self, x, par ):
    A=2.5066272 # sqrt(2*3.14159)
    gaus=(1.0/(par[2]*A)) * par[0]*exp((-0.5)*(x-par[1])*(x-par[1])/(par[2]*par[2]))
    return gaus;
 

class Massimo2:
   def __call__( self, x, par ):
      xx=x[0]-1.435944
      f3g2=0.
      xcost=1.0/math.sqrt(2.0*math.acos(-1.))
      if ( par[7] < 0): return 1.e30;
      if (xx > 0.0 ):
          f3g2=par[6]*(xx**par[7])*(1+par[8]*xx);
          f3g2=f3g2+(par[0]*xcost*exp(-((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2]) ));
          f3g2=f3g2+(par[3]*xcost*exp(-((x[0]-par[4])*(x[0]-par[4]))/(2*par[5]*par[5]) ));
          return f3g2;
