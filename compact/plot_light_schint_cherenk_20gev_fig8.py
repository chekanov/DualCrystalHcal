import sys
from globalmod import *

myinput="interactive"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]


print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)')

############# Configs ##############
nameX="S/E"
nameY="C/E"


Xmin=0
Xmax=1.1
Ymin=0.0
Ymax=1.1 


import sys
epsfig="figs/"+sys.argv[0].replace(".py",".eps")


gROOT.SetStyle("Plain");
xwin=650
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
c1.SetRightMargin(0.07)
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

name="schint_cherenk_calib_rel"
pp="#pi^{-}"
e="20"
fname="histos/hist_pi-_"+str(e)+"gev.root"
xfile=TFile( fname )
hist1=xfile.Get( name )
hist1.SetDirectory(0)
xfile.Close()

hist1.SetMarkerColor(2);
hist1.SetMarkerSize(0.2);
hist1.Draw("SCAT same")

################### electrons ##########################
e="20"
fname="histos/hist_e-_"+str(e)+"gev.root"
xfile=TFile( fname )
hist1a=xfile.Get( name )
hist1a.SetDirectory(0)
xfile.Close()

hist1a.SetMarkerColor(1);
hist1a.SetMarkerSize(0.2)
hist1a.Draw("same")

leg2=TLegend(0.6, 0.8, 0.89, 0.94);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.03);
#leg2.AddEntry(hh1,"Cherenkov light","pl")
#leg2.AddEntry(hh2,"Scintillation light","pl")
leg2.Draw("same");


histAll=hist1.Clone()
histAll.Add(hist1a)

f1 = TF1("g0","pol1",Xmin,Xmax);
f1.SetLineColor(4)
histAll.Fit(f1,"RM0");
f1.Draw("same")


f2=TF1("g1","pol1",Xmin,Xmax);
f2.SetLineStyle(2)
f2.SetLineWidth(2)
f2.SetLineColor(1)
f2.SetParameter(0,0.0)
f2.SetParameter(1,1.0)
f2.Draw("same")


# f(x) = p0 + p1*x
par1 = f1.GetParameters()
p0='%.3f'%( par1[0] )
p1='%.3f'%( par1[1] )
s1=p1+" x  (S/E) "+p0

leg2=TLegend(0.16, 0.6, 0.5, 0.94);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.05);
leg2.SetHeader("E=20 GeV")
leg2.AddEntry(f1,s1,"l")
leg2.AddEntry(hist1,"#pi^{-}","fp")
leg2.AddEntry(hist1a,"e^{-}","fp")
leg2.Draw("same");

# myText(0.2,0.7,2,0.05,pp+" 1,5,10,20 GeV")



gPad.RedrawAxis()
c1.Update()


if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);

