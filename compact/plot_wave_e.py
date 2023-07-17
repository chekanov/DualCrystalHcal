import sys
from globalmod import *

myinput="interactive"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]


print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)')

############# Configs ##############
nameX="W [nm]"
nameY="Events / Total"


Xmin=270
Xmax=640.0
Ymin=0.0
Ymax=0.279


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
c1.SetTopMargin(0.06)
c1.SetRightMargin(0.1)
c1.SetFillColor(0)
c1.SetLeftMargin(0.15)

c1.SetLogy(0)
#c1.SetLogx(1)

h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
h.Draw()

ax=h.GetXaxis();
ax.SetTitle( nameX );
ay=h.GetYaxis();
ay.SetTitle( nameY );

ax.SetLabelSize(0.04)
ax.SetTitleSize(0.05)
ax.SetTitleOffset(1.0)
ay.SetTitleOffset(1.3)

ay.SetLabelSize(0.04)
ay.SetTitleSize(0.05)
gPad.SetTickx()
gPad.SetTicky()

ay.Draw("same")
ax.Draw("same")

e="1"
fname="histos/hist_e-_"+str(e)+"gev.root"
xfile=TFile( fname )
hh1=xfile.Get("wave_cherenk")
hh2=xfile.Get("wave_scintil") 

hh1.SetMarkerColor( 1 )
hh1.Scale(1.0/hh1.Integral())
hh1.SetLineColor( 1 )
hh1.Draw("same h")

hh2.SetLineStyle( 2 )
hh2.Scale(  1.0/hh2.Integral() )
hh2.SetMarkerColor( 2 )
hh2.SetLineColor( 2 )
hh2.Draw("same h")



leg2=TLegend(0.5, 0.75, 0.89, 0.9);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.03);
leg2.AddEntry(hh1,"Cherenkov light","pl")
leg2.AddEntry(hh2,"Scintillation light","pl")
leg2.Draw("same");


gPad.RedrawAxis()
c1.Update()


if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);

