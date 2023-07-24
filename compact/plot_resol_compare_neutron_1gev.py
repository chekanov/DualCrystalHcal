import sys
from globalmod import *

myinput="interactive"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]


print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)')

############# Configs ##############
nameX="E_{meas}/E"
nameY="Events"

Xmin=0.
Xmax=2.0 
Ymin=0.
Ymax=549

def getRes(hh):
      RMS=rms90(hh);
      #RMS=hh.GetRMS()
      ERR=hh.GetRMSError()
      MEAN=hh.GetMean()
      MEAN_ERR=hh.GetMeanError()
      res=RMS/MEAN
      err=MEAN_ERR/MEAN
      a1='%.2f'%( res )
      b1='%.2f'%( err )
      return a1+"#pm"+b1
 

import sys
epsfig="figs/"+sys.argv[0].replace(".py",".eps")


gROOT.SetStyle("Plain");
xwin=600
ywin=600
c1=TCanvas("c","Mass",10,10,xwin,ywin);
c1.SetFrameBorderMode(0);
ps1 = TPostScript( epsfig,113)
c1.SetTickx()
c1.SetTicky()
c1.SetTitle("")
c1.SetLineWidth(3)
c1.SetBottomMargin(0.12)
c1.SetTopMargin(0.05)
c1.SetRightMargin(0.1)
c1.SetFillColor(0)
c1.SetLeftMargin(0.15)

# c1.SetLogy(1)
c1.SetLogx(0)

h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
h.Draw()

ax=h.GetXaxis();
ax.SetTitle( nameX );
ay=h.GetYaxis();
ay.SetTitle( nameY );

ax.SetLabelSize(0.04)
ax.SetTitleSize(0.05)
ax.SetTitleOffset(1.0)
ay.SetTitleOffset(1.38)

ay.SetLabelSize(0.04)
ay.SetTitleSize(0.05)
gPad.SetTickx()
gPad.SetTicky()
ay.Draw("same")
ax.Draw("same")

hfile=TFile("histos/hist_neutron_1gev.root")
histo=hfile.Get("heest_meass")
histo1=hfile.Get("heest_scint")
histo2=hfile.Get("heest_cherenk")

histo.SetMarkerColor(1)
histo.SetLineWidth(3)
histo.SetMarkerSize(1.2)
histo.SetMarkerStyle(21)
histo.Draw("same pe")

histo1.SetLineStyle(1)
histo1.SetLineWidth(3)
histo1.SetMarkerColor(2)
histo1.SetMarkerSize(1.2)
histo1.SetMarkerStyle(22)
histo1.SetLineColor(2)
histo1.Draw("same pe")

histo2.SetLineWidth(3)
histo2.SetMarkerColor(3)
histo2.SetMarkerSize(1.2)
histo2.SetMarkerStyle(20)
histo2.SetLineColor(3)
histo2.Draw("same pe")


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


b1=getFit(histo1, color=2)
print("Scint: Width/Mean=",b1)

b2=getFit(histo2, color=3)
print("Cheren: Width/Mean=",b2)

b3=getFit(histo, color=1)
print("Scint + Cheren: Width/Mean=",b3)


leg2=TLegend(0.18, 0.68, 0.89, 0.94);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.033);
leg2.SetHeader("1 GeV neutrons")
#leg2.AddEntry(histo,"Cherenk+Scint: RMS90/E="+getRes(histo),"l")
#leg2.AddEntry(histo1,"Scint: RMS90/E=" + getRes(histo1),"l")
#leg2.AddEntry(histo2,"Cherenk: RMS90/E="+getRes(histo2) ,"l")
leg2.AddEntry(histo1,"Scint: #sigma/E=" + b2,"lp")
leg2.AddEntry(histo2,"Cherenk: #sigma/E="+b3 ,"lp")
leg2.AddEntry(histo,"Cherenk+Scint: #sigma/E="+b1,"lp")

leg2.Draw("same");


gPad.RedrawAxis()
c1.Update()


if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);

