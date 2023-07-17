import sys
from globalmod import *
import random

#####################

myinput="interactive"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]


print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)')



############# Configs ##############
nameX="Energy [GeV]"
nameY="p0 value"

Xmin=0.0001 
Xmax=50 
Ymin=0. 
Ymax=80*1000-10 

ptmin=0.016
ptmax=40

import sys
epsfig="figs/"+sys.argv[0].replace(".py",".eps")

gROOT.SetStyle("Plain");
xwin=700
ywin=600
c1=TCanvas("c","Mass",10,10,xwin,ywin);
c1.SetFrameBorderMode(0);
c1.Divide(2,2);
#c1.Divide(2,4,0.02,0.02);
gPad.SetLeftMargin(0.25);
c1.SetLeftMargin(0.25);
gPad.SetRightMargin(0.1);
c1.SetRightMargin(0.1);

ps1 = TPostScript( epsfig,113)
c1.SetTickx()
c1.SetTicky()
c1.SetTitle("")
c1.SetLineWidth(3)
c1.SetFillColor(0)


# particles=["e-","pi0","pi-","kaon-", "neutron","proton"] 
# particles=["pi-","kaon-", "neutron","proton"] 

Particle="pi-"
IND=0
# get index files
if (Particle == "pi"): IND=0
if (Particle == "K"): IND=1
if (Particle == "n"): IND=2
if (Particle == "p"): IND=3

legends=[]
functions1=[]
functions2=[]

#eq= "[0]+ [1]*(x) + [2]*(x*x)+[3]*(x*x*x)"
#eq= "[0]*sqrt([1]*(x)) + [2]*(x)"
#eq= "[0]*pow(x,[1])+[2]*x"
#eq= "[0]*pow(x,[1]+[2]*x)"
eq= "[0]+x*[1]+[2]*x*x"
#eq= "[0]+x*[1]"


############ scint 

scint={}
Energy=1
xfileS=open("figs/gev"+str(Energy)+"_scintl_mc_em_fit.txt","r")
with xfileS as file:
    lines = [line.rstrip() for line in file]
scint[Energy]=lines
Energy=5
xfileS=open("figs/gev"+str(Energy)+"_scintl_mc_em_fit.txt","r")
with xfileS as file:
    lines = [line.rstrip() for line in file]
scint[Energy]=lines
Energy=10
xfileS=open("figs/gev"+str(Energy)+"_scintl_mc_em_fit.txt","r")
with xfileS as file:
    lines = [line.rstrip() for line in file]
scint[Energy]=lines
Energy=20
xfileS=open("figs/gev"+str(Energy)+"_scintl_mc_em_fit.txt","r")
with xfileS as file:
    lines = [line.rstrip() for line in file]
scint[Energy]=lines


cheren={}
############ cherenkov
Energy=1
xfileC=open("figs/gev"+str(Energy)+"_cheren_mc_em_fit.txt","r")
with xfileC as file:
    lines = [line.rstrip() for line in file]
cheren[Energy]=lines
Energy=5
xfileC=open("figs/gev"+str(Energy)+"_cheren_mc_em_fit.txt","r")
with xfileC as file:
    lines = [line.rstrip() for line in file]
cheren[Energy]=lines
xfileC=open("figs/gev"+str(Energy)+"_cheren_mc_em_fit.txt","r")
Energy=10
with xfileC as file:
    lines = [line.rstrip() for line in file]
cheren[Energy]=lines
xfileC=open("figs/gev"+str(Energy)+"_cheren_mc_em_fit.txt","r")
Energy=20
with xfileC as file:
    lines = [line.rstrip() for line in file]
cheren[Energy]=lines


def fillCher(pp, Energy,kk, index):
   cy=cheren[Energy]
   cy=cy[IND] # pion
   cy=cy.split()
   d1=1 # value
   d2=4 # errpor
   if index==0:
          d1=1 # value
          d2=4 # errpor
   if index==1:
          d1=2 # value
          d2=5 # errpor
   if index==2:
          d1=3 # value
          d2=6 # errpor
   pp.SetPoint(kk,Energy,float(cy[d1]))
   pp.SetPointError(kk,0,float(cy[d2]));

def fillScin(pp, Energy,kk, index):
   cy=scint[Energy]
   cy=cy[IND] # pion
   cy=cy.split()
   d1=1 # value
   d2=4 # errpor
   if index==0:
          d1=1 # value
          d2=4 # errpor
   if index==1:
          d1=2 # value
          d2=5 # errpor
   if index==2:
          d1=3 # value
          d2=6 # errpor
   pp.SetPoint(kk,Energy,float(cy[d1]))
   pp.SetPointError(kk,0,float(cy[d2]));


nameX="Energy [GeV]"

ID=0
c1.cd(ID+1);
nameY="p0 value"
################ get p0 
p0_c=TGraphErrors()
p0_c.SetMarkerColor(2)
p0_c.SetMarkerSize(1)
p0_c.SetMarkerStyle(20)
p0_c.SetLineColor(1)
fillCher(p0_c,1,0, ID)
fillCher(p0_c,5,1, ID)
fillCher(p0_c,10,2, ID)
fillCher(p0_c,20,3, ID)


p0_s=TGraphErrors()
p0_s.SetMarkerColor(3)
p0_s.SetMarkerSize(1)
p0_s.SetMarkerStyle(20)
p0_s.SetLineColor(1)
fillScin(p0_s,1,0, ID)
fillScin(p0_s,5,1, ID)
fillScin(p0_s,10,2, ID)
fillScin(p0_s,20,3, ID)

gPad.SetLeftMargin(0.22)
gPad.SetBottomMargin(0.18);
gPad.SetTopMargin(0.05);
gPad.SetRightMargin(0.02);

h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
h.Draw()
ax=h.GetXaxis();
ax.SetTitle( nameX );
ay=h.GetYaxis();
ay.SetTitle( nameY );
p0_s.Draw("same pe")
p0_c.Draw("same pe")


ID=1
Ymax=10000-10
c1.cd(ID+1);
nameY="p1 value"
################ get p0 
p1_c=TGraphErrors()
p1_c.SetMarkerColor(2)
p1_c.SetMarkerSize(1)
p1_c.SetMarkerStyle(20)
p1_c.SetLineColor(1)
fillCher(p1_c,1,0, ID)
fillCher(p1_c,5,1, ID)
fillCher(p1_c,10,2,ID )
fillCher(p1_c,20,3, ID)

p1_s=TGraphErrors()
p1_s.SetMarkerColor(3)
p1_s.SetMarkerSize(1)
p1_s.SetMarkerStyle(20)
p1_s.SetLineColor(1)
fillScin(p1_s,1,0, ID)
fillScin(p1_s,5,1, ID)
fillScin(p1_s,10,2, ID)
fillScin(p1_s,20,3, ID)

gPad.SetLeftMargin(0.22)
gPad.SetBottomMargin(0.18);
gPad.SetTopMargin(0.05);
gPad.SetRightMargin(0.02);

h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
h.Draw()
ax=h.GetXaxis();
ax.SetTitle( nameX );
ay=h.GetYaxis();
ay.SetTitle( nameY );
p1_s.Draw("same pe")
p1_c.Draw("same pe")


ID=2
Ymin=-500
Ymax=500-1
c1.cd(ID+1);
nameY="p2 value"
################ get p0 
p2_c=TGraphErrors()
p2_c.SetMarkerColor(2)
p2_c.SetMarkerSize(1)
p2_c.SetMarkerStyle(20)
p2_c.SetLineColor(1)
fillCher(p2_c,1,0, ID)
fillCher(p2_c,5,1, ID)
fillCher(p2_c,10,2, ID)
fillCher(p2_c,20,3, ID)

p2_s=TGraphErrors()
p2_s.SetMarkerColor(3)
p2_s.SetMarkerSize(1)
p2_s.SetMarkerStyle(20)
p2_s.SetLineColor(1)
fillScin(p2_s,1,0, ID)
fillScin(p2_s,5,1, ID)
fillScin(p2_s,10,2, ID)
fillScin(p2_s,20,3, ID)

gPad.SetLeftMargin(0.22)
gPad.SetBottomMargin(0.18);
gPad.SetTopMargin(0.05);
gPad.SetRightMargin(0.02);

h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
h.Draw()
ax=h.GetXaxis();
ax.SetTitle( nameX );
ay=h.GetYaxis();
ay.SetTitle( nameY );
p2_s.Draw("same pe")
p2_c.Draw("same pe")



  #myTextUser(0.58,700,1,0.07,"<S>="+s0)
  #myTextUser(0.58,500,1,0.07,"<C>="+s1)

gPad.RedrawAxis()
c1.Update()


if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);

