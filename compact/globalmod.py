# Global module

from math import sqrt,acos,exp
from array import array
from ROOT import gROOT,gPad,gStyle,TCanvas,TSpline3,TFile,TLine,TLatex,TAxis,TLegend,TPostScript
from ROOT import TH2D,TArrow,TCut,TPad,TH1D,TF1,TObject,TPaveText,TGraph,TGraphErrors,TGraphAsymmErrors
from ROOT import TMath, TGraph2D,TTree,TMultiGraph,TBranch,gSystem,gDirectory
from ROOT import TPaveStats
import random

# par[1] - mean
# par[2] - sigma
class Gauss:
   def __call__( self, x, par ):
       out=par[0] * TMath.Gaus(x[0],par[1],par[2])
       return out;

# convert histogram to TGraph
# also apply scale if needed
def TH1toTGraphError(h1, shiftX=1.0):

    g1 = TGraphErrors()
    for i in range(h1.GetNbinsX()):
        y = h1.GetBinContent(i+1)
        ey = h1.GetBinError(i+1)
        x = h1.GetBinCenter(i+1)
        ex = h1.GetBinWidth(i+1)/2.0
        g1.SetPoint(i, x*shiftX, y)
        g1.SetPointError(i, ex, ey)

    g1.SetMarkerColor(1)
    g1.SetMarkerStyle(20)
    g1.SetMarkerSize(0.5)

    # g1->Print();
    return g1


def TH1correctError(h1, scale):
    for i in range(h1.GetNbinsX()):
        y = h1.GetBinContent(i+1)
        ey=h1.GetBinError(i+1)
        x=h1.GetBinCenter(i+1)
        # ex=h1.GetBinWidth(i+1)
        # h1.SetBinContent(i+1,y)
        # h1.SetBinError(i+1,ey)
        if (x>0.1): continue
        if (y > 0):
            h1.SetBinError(i+1, ey*scale )
        else:
            h1.SetBinError(i+1, 0)

def getFit(hh, color=1):
    mean=hh.GetMean()
    sigma=hh.GetRMS()
    f1 = TF1("g0","gaus",mean-3*sigma, mean+3*sigma);
    f1.SetParameter(0,200)
    f1.SetParameter(1,mean)
    f1.SetParameter(2,sigma)
    f1.SetLineColor(color)
    hh.Fit(f1,"RM");
    par1 = f1.GetParameters()
    err1 = f1.GetParErrors()
    b1='%.3f #pm %.3f'%( par1[2]/par1[1], err1[2]/par1[1]  )
    return b1


# 90% RMS
# http://agenda.linearcollider.org/event/2703/contributions/8926/attachments/6867/11488/RMS90.pdf
# https://pdfs.semanticscholar.org/d98d/1dbf35236e45edbc4dc68aede5db6a82d294.pdf
# RMS90 is defined as the RMS falculated for a
# distribution in which 10% of the events in the tails
# are excluded from the RMS calculation
def rms90(h):
   axis = h.GetXaxis();
   nbins = axis.GetNbins();
   imean = axis.FindBin(h.GetMean());
   entries =0.9*h.GetEntries();
   w = h.GetBinContent(imean);
   x = h.GetBinCenter(imean);
   sumw = w;
   sumwx = w*x;
   sumwx2 = w*x*x;
   for i in range(1,nbins):
      if (i> 0):
         w = h.GetBinContent(imean-i);
         x = h.GetBinCenter(imean-i);
         sumw += w;
         sumwx += w*x;
         sumwx2 += w*x*x;
      if (i<= nbins):
         w = h.GetBinContent(imean+i);
         x = h.GetBinCenter(imean+i);
         sumw += w;
         sumwx += w*x;
         sumwx2 += w*x*x;
      if (sumw > entries): break;
   x = sumwx/sumw;
   rms2 = TMath.Abs(sumwx2/sumw -x*x)
   result = TMath.Sqrt(rms2)
   # print ("RMS of central 90% = ",result, " RMS total =",h.GetRMS());
   return 1.25*result


def myText(x,y,color=1,size=0.08,text=""):
  l=TLatex()
  l.SetTextSize(size);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);

def myTextUser(x,y,color=1,size=0.08,text=""):
  l=TLatex()
  l.SetTextSize(size);
  #l.SetUser();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);

