# S.Chekanov. ANL

import sys
from array import array
from ROOT import gROOT,TMath  
from ROOT import TH2F,TCut,TPad,TH1F,TF1,TObject
from ROOT import gSystem,gDirectory
from math import *


# square root of s value
s=7000

# fit functions
class Bkg8: #from "Search for New Physics in Dijet Distributions with the ATLAS Detector" by The ATLAS Collaboration
   def __call__( self, x, par ):
       out= (par[0]*(1-(x[0]/s))**par[1])*(  (x[0]/s)**( par[2]+par[3]*log((x[0]/s)) )  )
       return out;


class Bkg8a: #from "Search for New Physics in Dijet Distributions with the ATLAS Detector" by The ATLAS Collaboration
   def __call__( self, x, par ):
       out= (par[0]*(1-(x[0]/s))**par[1]) * exp(-par[2]*x[0]) 
       return out;


class Bkg7: #from "Top Jets at the LHC" by Almeida, Lee, Perez, Sung, and Virzi
   def __call__( self, x, par ):
       out= (par[0]/x[0])*log(par[1]/x[0],10)
       return out;

class Bkg6:
   def __call__( self, x, par ):
       out= ( (par[0]*(x[0]**-par[1]) )*((1-(x[0]/s))**par[2]) )*( 1+(par[3]*(x[0]/s))+(par[4]*((x[0]/s)**2))+(par[5]*((x[0]/s)**3))+(par[6]*((x[0]/s)**4)) )
       return out;

class Bkg5:
   def __call__( self, x, par ):
       out= ( (par[0]*(x[0]**-par[1]) )*((1-(x[0]/s))**par[2]) )*( 1+(par[3]*(x[0]/s))+(par[4]*((x[0]/s)**2))+(par[5]*((x[0]/s)**3)) )
       return out;

class Bkg4:
   def __call__( self, x, par ):
       out= ( (par[0]*(x[0]**-par[1]) )*((1-(x[0]/s))**par[2]) )*(1+(par[3]*(x[0]/s))+(par[4]*((x[0]/s)**2)))
       return out;

class Bkg3:
   def __call__( self, x, par ):
       out= ( (par[0]*(x[0]**-par[1]) )*((1-(x[0]/s))**par[2]) )*(1+(par[3]*(x[0]/s)))
       return out;


class NLO:
      def __call__( self, x, par ):
         if (x[0]>1):
           return (par[0]/x[0]) * log(par[1]/x[0])
         return 0;

class NLO1:
      def __call__( self, x, par ):
         if (x[0]>1):
           return (par[0]/x[0]) * log(par[1]/x[0]) + par[2]*exp(-par[3]*x[0]) 
         return 0;


class Bkg2:
   def __call__( self, x, par ):
       out= (par[0]*(x[0]**-par[1]) )*((1-(x[0]/s))**par[2])
       return out;


class Bkg1:
   def __call__( self, x, par ):
       out= par[0]*(x[0]**(-par[1]))
       return out;

class BPoly:
   def __call__( self, x, par ):
       out= par[0]*x[0]*x[0]*x[0]-par[1]*x[0]*x[0]+par[2]*x[0]+par[3] 
       return out;



class Bkg1Gauss:
   def __call__( self, x, par ):
       out= par[0]*(x[0]**(-par[1])) + par[2]*TMath.Gaus(x[0],par[3],par[4]) 
       return out;


class Bkg11:
   def __call__( self, x, par ):
       out= par[0]*(1-x[0]+par[1]*x[0]*x[0])*(x[0]**(-par[2]))
       return out;

class Bkg11Gauss:
   def __call__( self, x, par ):
       out= par[0]*(1-x[0]+par[1]*x[0]*x[0])*(x[0]**(-par[2])) +  par[3]*TMath.Gaus(x[0],par[4],par[5]) 
       return out;

class Bkg13:
   def __call__( self, x, par ):
       # out= par[0]*(1-par[1]*x[0])*(x[0]**(-par[2]))
       x1=par[0]*(1-par[1]*x[0])
       x2=TMath.Power(x[0],-par[2])
       return x1*x2;

class Bkg14:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0],-par[1]) 
       x2=TMath.Exp(-par[2]*x[0])
       return x1*x2;

class Bkg15:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0],-par[1])
       x2=TMath.Exp(-par[2]*x[0])
       x3=TMath.Exp(-par[3]*x[0]*x[0])
       # x3=TMath.Power(x[0]*x[0],-par[1])
       # x3=TMath.Exp(-par[3]*TMath.Sqrt(x[0]))
       # x3=par[3]*TMath.Gaus(x[0],par[4],par[5])
       return x1*x2*x3;


class Bkg15Gauss:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0],-par[1])
       x2=TMath.Exp(-par[2]*x[0])
       x3=TMath.Exp(-par[3]*x[0]*x[0])
       return x1*x2*x3+par[4]*TMath.Gaus(x[0],par[5],par[6]);



class Bkg14UA1:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0],-par[1])
       x2=TMath.Exp(-par[2]*x[0])
       x3=TMath.Exp(-par[3]*(TMath.Power(x[0],2)))
       return x1*x2*x3;


# modified UA1
class Bkg16:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0]-60,-par[1])
#       x2=TMath.Exp(-1.0*(TMath.Power(par[2]*(x[0]-60),2)))
#       x2=TMath.Exp(-par[2]*(x[0]-60))*TMath.Exp(-1.0*(TMath.Power(par[3]*(x[0]-60),2)))
       x2=TMath.Exp(-1.0*(TMath.Power(par[2]*(x[0]-60),2)))
#       x2=TMath.Exp(-par[2]*(x[0]-60))


       return x1*x2;

class Bkg13Gauss:
   def __call__( self, x, par ):
       out= par[0]*(1-par[1]*x[0])*(x[0]**(-par[2])) + par[3]*TMath.Gaus(x[0],par[4],par[5]) 
       return out;

class Bkg16Gauss:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0]-60,-par[1]) 
       x2=TMath.Exp(-1.0*(TMath.Power(par[2]*(x[0]-60),2)))
       out= x1*x2 + par[3]*TMath.Gaus(x[0],par[4],par[5])
       return out;

class Bkg14Gauss:
   def __call__( self, x, par ):
       x1=par[0]*TMath.Power(x[0],-par[1])
       x2=TMath.Exp(-par[2]*x[0])
       out= x1*x2 + par[3]*TMath.Gaus(x[0],par[4],par[5])
       return out;


class Bkg12:
   def __call__( self, x, par ):
       out=0.
       xx=x[0]/7000.
       t=1.0-par[1]*xx*xx
       if (t>0):
            out= par[0]* (t**par[2]) *(x[0]**(-par[3]))

       return out;


class Bkg:
   def __call__( self, x, par ):
       out= par[0]+ (par[1]/x[0]) * log(1.0/x[0])
       return out;

class P1:
   def __call__( self, x, par ):
       out= par[0]+ par[1]*x[0]
       return out;

class P2:
   def __call__( self, x, par ):
       out= par[0]+ par[1]*x[0]+par[2]*x[0]*x[0] 
       return out;

class Exp:
   def __call__( self, x, par ):
       out=0
       try :
          out= par[0]* exp(-par[1]*x[0]) 
       except :
          pass
       return out;

class Exp2:
   def __call__( self, x, par ):
       out= par[0]* exp(-par[1]*x[0]*x[0])
       return out;

class ExpLog:
   def __call__( self, x, par ):
       out=0
       if x[0]>20:
            out= par[0]* exp(-par[1]*x[0])* (par[2]/x[0]);  # log (par[2]/(x[0]*x[0])) 
       return out;


class Expp:
   def __call__( self, x, par ):
       out= par[0]* exp(-par[1]*x[0])+par[2]*x[0]
       return out;

class GaussDouble:
   def __call__( self, x, par ):

       out=0.0
       gaus1=0.0
       gaus2=0.0
       gaus1=par[0]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]))
       gaus2=par[4]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[3]*par[3]))
       cc=par[2]*par[3]
       if cc>0:
                       out=gaus1+gaus2;
       return out

"""
class CrystalBall:
   def __call__( self, x, par ):
      out=0
      t = (x[0]-par[2])/par[3]
      if (par[0] < 0): t = -t;
      absAlpha = TMath.Abs(par[0]);
      if (t>=-absAlpha): 
        out=par[4]*TMath.Exp(-0.5*t*t);
      else: 
        a =TMath.Power(par[1]/absAlpha,par[1])*TMath.Exp(-0.5*absAlpha*absAlpha);
        b= par[1]/absAlpha - absAlpha;
        out=par[4]*(a/TMath.Power(b - t, par[1]));
      return out
"""

# http://en.wikipedia.org/wiki/Crystal_Ball_function
# par[0] - normalisation
# par[1] - mean
# par[2] - sigma
# par[3] - alpha
# par[4]  - N 

class CrystalBall:
   def __call__( self, x, par ):
      out=0
      t = (x[0]-par[1])/par[2]
      if (par[0] < 0): t = -t;
      absAlpha = TMath.Abs(par[3]);
      if (t>=-absAlpha):
        out=par[0]*TMath.Exp(-0.5*t*t);
      else:
        a =TMath.Power(par[4]/absAlpha,par[4])*TMath.Exp(-0.5*absAlpha*absAlpha);
        b= par[4]/absAlpha - absAlpha;
        out=par[0]*(a/TMath.Power(b - t, par[4]));
      return out


class Bkg1CrystalBall:
   def __call__( self, x, par ):

      out=0
      t = (x[0]-par[3])/par[4]
      if (par[0] < 0): t = -t;
      absAlpha = TMath.Abs(par[5]);
      if (t>=-absAlpha):
        out=par[2]*TMath.Exp(-0.5*t*t);
      else:
        a =TMath.Power(par[6]/absAlpha,par[6])*TMath.Exp(-0.5*absAlpha*absAlpha);
        b= par[6]/absAlpha - absAlpha;
        out=par[2]*(a/TMath.Power(b - t, par[6]));

      return par[0]*(x[0]**(-par[1])) + out 


class Bkg1CrystalBall:
   def __call__( self, x, par ):

      out=0
      t = (x[0]-par[3])/par[4]
      if (par[0] < 0): t = -t;
      absAlpha = TMath.Abs(par[5]);
      if (t>=-absAlpha):
        out=par[2]*TMath.Exp(-0.5*t*t);
      else:
        a =TMath.Power(par[6]/absAlpha,par[6])*TMath.Exp(-0.5*absAlpha*absAlpha);
        b= par[6]/absAlpha - absAlpha;
        out=par[2]*(a/TMath.Power(b - t, par[6]));

      return par[0]*(x[0]**(-par[1])) + out







# http://en.wikipedia.org/wiki/Crystal_Ball_function
# par[0] - normalisation
# par[1] - mean
# par[2] - sigma
# par[3] - alpha
# par[4]  - N 
# par[5]  - norm Gauss2 
# par[6]  - sigma Gauss2 

class CrystalBallGauss:
   def __call__( self, x, par ):
      out=0
      t = (x[0]-par[1])/par[2]
      if (par[0] < 0): t = -t;
      absAlpha = TMath.Abs(par[3]);
      if (t>=-absAlpha):
        out=par[0]*TMath.Exp(-0.5*t*t);
      else:
        a =TMath.Power(par[4]/absAlpha,par[4])*TMath.Exp(-0.5*absAlpha*absAlpha);
        b= par[4]/absAlpha - absAlpha;
        out=par[0]*(a/TMath.Power(b - t, par[4]));
      out=out+par[5]*TMath.Gaus(x[0],par[1],par[6]) 
      return out



class GaussAss:
   def __call__( self, x, par ):

       out=0.0
       gaus1=0.0
       gaus2left=0.0
       gaus2right=0.0
       gaus1=par[0]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]))
       gaus2left=par[5]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[3]*par[3]))
       gaus2right=par[5]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[4]*par[4]))

       cc=par[2]*par[3]*par[4]

       if x[0]<=par[1] and cc>0:
                       out=gaus1+gaus2left;
       if x[0]>par[1] and cc>0:
                       out=gaus1+gaus2right;

       return out


class BkgP1:
   def __call__( self, x, par ):
       out=par[0] + par[1]*x[0];
       return out;

class BkgP2:
   def __call__( self, x, par ):
       out=par[0] + par[1]*x[0]+par[2]*par[2]*x[0];
       return out;


class GaussAssBkg:
   def __call__( self, x, par ):

       out=0.0
       gaus1=0.0
       gaus2left=0.0
       gaus2right=0.0

       gaus1=par[0]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]))
       gaus2left=par[5]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[3]*par[3]))
       gaus2right=par[5]*TMath.Exp((-0.5)*(x[0]-par[1])*(x[0]-par[1])/(par[4]*par[4]))

       cc=par[2]*par[3]*par[4]

       if x[0]<=par[1] and cc>0:
                       out=gaus1+gaus2left;
       if x[0]>par[1] and cc>0:
                       out=gaus1+gaus2right;


       back = par[6] + par[7]*x[0];
       # add backg 
       backout=out+ back;

       return backout;

# const, mean, gamma
class BW:
   def __call__( self, x, par ):
       out=par[0] * TMath.BreitWigner(x[0],par[1],par[2])  
       return out;

class BWbkg:
   def __call__( self, x, par ):
       out= ( par[0]*TMath.BreitWigner(x[0],par[1],par[2]) + par[3]*(x[0]**-par[4]) ) 
       return out;

class GaussBkg:
   def __call__( self, x, par ):
       out= ( par[0]*TMath.Gaus(x[0],par[1],par[2])) +  par[3]*(x[0]**-par[4]);  
       return out;

# par[1] - mean
# par[2] - sigma
class Gauss:
   def __call__( self, x, par ):
       out=par[0] * TMath.Gaus(x[0],par[1],par[2])
       return out;

