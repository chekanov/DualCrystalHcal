#!/usr/bin/env python
import sys
sys.path.append("./modules/")
from initialize  import *

# is it batch?
# if batch, set input to batch
myinput="no any input"
if (len(sys.argv) > 1):
   myinput = sys.argv[1]
print "Mode=",myinput

ymin=0.0
ymax=5000;


# gSystem.Load('macros/Compare_C')
# from ROOT import Compare 
# m=Compare() # first histogram is data!



epsfig="figs/see_kshort.eps"
# file
f1=TFile( "out/data2004.root" )



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

class Bkg:
   def __call__( self, x, par ):
       out=par[0] + par[1]*x[0]; 
       return out;


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


MyMinX=0.46
MyMaxX=0.56
# background
back=TF1("back",Bkg(),MyMinX,MyMaxX,2);
back.SetNpx(100); back.SetLineColor(5); back.SetLineStyle(3)
back.SetLineWidth(3)


peakback=TF1("peakback",GaussAssBkg(),MyMinX,MyMaxX,8)
peakback.SetNpx(100); peakback.SetLineColor(2); peakback.SetLineStyle(2)
peakback.SetLineWidth(3); peakback.SetParNames("Const","Peak","Sigma","Sigma-left","Sigma-right", "Const1", "P0", "P1" )


peak=TF1("peak",GaussAss(),MyMinX,MyMaxX,6)
peak.SetNpx(100); peak.SetLineColor(2); peak.SetLineStyle(2)
peak.SetLineWidth(3); peak.SetParNames("Const","Peak","Sigma","Sigma-left","Sigma-right","Const1");


gROOT.Reset(); 
gROOT.SetStyle("Plain");
gStyle.SetEndErrorSize(0.0); gStyle.SetErrorX(0);
gStyle.SetTitleX(0.1)
gStyle.SetTitleW(0.8)
gStyle.SetTitleBorderSize(0)
gStyle.SetTitleBorderSize(0)

gStyle.SetTitleSize(0.05,"x")
# gStyle.SetTitleSize(0.06,"y")
# gStyle.SetTitleX(0.2) 
# gStyle.SetTitleH(0.07)
# gStyle.SetTitleW(0.3);

gStyle.SetTitleColor(4);

gStyle.SetLabelSize(0.045,"y");
gStyle.SetLabelOffset(0.005,"y");
gStyle.SetLabelSize(0.045,"x");
gStyle.SetLabelOffset(0.005,"x");


c1=TCanvas("c","BPRE",10,10,600,400);
ps1 = TPostScript( epsfig,113)
# c1.Divide(2,2,0.005,0.005);



gStyle.SetOptStat(1111111);
gStyle.SetStatY(0.9);                
gStyle.SetStatX(0.9);                
# gStyle.SetStatW(0.4);                
# gStyle.SetStatH(0.2);      



gPad.SetLeftMargin(0.2);
# c1.SetGridx(2);
c1.SetTickx()
c1.SetTicky()
c1.SetTitle("")
c1.SetLineWidth(3)
c1.SetBottomMargin(0.1)
c1.SetTopMargin(0.1)
c1.SetFillColor(0)


name="kshort"
c1.cd(1);
h=f1.Get(name)
h.SetTitle("2004");
# h.Draw()

back.SetParameter(0,300)
back.SetParameter(1,0.4)
# h.Fit(back,'R0')
# par = back.GetParameters()

peakback.SetParameter(0,100)
peakback.SetParameter(1,0.48)
peakback.SetParLimits(1,0.46,0.51)
peakback.SetParameter(2,0.004)
peakback.SetParameter(3,0.004)
peakback.SetParameter(4,0.004)
peakback.SetParameter(5,10.0)
# peakback.SetParameter(6,par[0])
# peakback.SetParameter(7,par[1])

h.Fit(peakback,'RV+','pe')


ax=h.GetXaxis(); ax.SetTitle( "M (GeV)" );
# h.SetAxisRange(xmin,xmax,"x")
h.SetAxisRange(ymin,ymax,"y")
ax.SetTitleOffset(1.2)
ay=h.GetYaxis(); ay.SetTitle( "Events" );
ay.SetTitleOffset(1.2)
ax.Draw("same"); ay.Draw("same"); 


c1.Update()
ps1.Close()
if (myinput != "batch"):
              if (raw_input("Press any key to exit") != "-9999"): sys.exit(1);



