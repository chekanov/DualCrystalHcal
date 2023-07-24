import sys
from globalmod import *
import random

### define energy here:
Energy=10 
#####################

myinput="interactive"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]


print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)')



############# Configs ##############
nameX="E_{EM} [GeV]"
nameY="<Y>"

Xmin=1 
Xmax=50 
Ymin=10000. 
Ymax=200000*1000-10

ptmin=1
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
particles=["pi-","kaon-", "neutron","proton"] 

e=Energy;

legends=[]
functions1=[]
functions2=[]

#eq= "[0]+ [1]*(x) + [2]*(x*x)+[3]*(x*x*x)"
#eq= "[0]*sqrt([1]*(x)) + [2]*(x)"
#eq= "[0]*pow(x,[1])+[2]*x"
#eq= "[0]*pow(x,[1]+[2]*x)"
eq= "[0]+x*[1]+[2]*x*x"
#eq= "[0]+x*[1]"


for i in range(len(particles)):
  leg2=TLegend(0.24, 0.7, 0.5, 0.92);
  leg2.SetBorderSize(0);
  leg2.SetTextFont(62);
  leg2.SetFillColor(10);
  leg2.SetTextSize(0.07);
  legends.append(leg2)

  f1=TF1("f1",eq,ptmin,ptmax);
  f1.SetLineStyle(1)
  f1.SetLineWidth(2)
  f1.SetLineColor(1)
  f1.SetParameter(0,0.4)
  f1.SetParameter(1,0.01)
  f1.SetParameter(2,0.1)
  functions1.append(f1)

  f2=TF1("f2",eq,ptmin,ptmax);
  f2.SetLineStyle(1)
  f2.SetLineWidth(2)
  f2.SetLineColor(4)
  f2.SetParameter(0,0.004)
  f2.SetParameter(1,0.01)
  f2.SetParameter(2,0.1)
  functions2.append(f2)

xfileC=open("figs/gev"+str(Energy)+"_cheren_mc_em_fit.txt","w")
xfileS=open("figs/gev"+str(Energy)+"_scintl_mc_em_fit.txt","w")


errors1=[]
# get graphs with SD errors
for i in range( len(particles) ):
  part=particles[i]
  fname="histos/hist_"+part+"_"+str(e)+"gev.root"
  xfile0=TFile( fname )
  hh0sd=xfile0.Get("scint_vs_EM_E_before_SD")
  hh0sd.SetDirectory(0)
  g1=TH1toTGraphError(hh0sd.Clone())
  g1.SetFillColor(5)
  errors1.append(g1)
  xfile0.Close()


for i in range( len(particles) ):
  c1.cd(i+1);


  gPad.SetLeftMargin(0.22)
  print("Draw pad=",i)
  gPad.SetBottomMargin(0.18);
  gPad.SetTopMargin(0.05);
  gPad.SetRightMargin(0.02);
  gPad.SetLogy(1)

  h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
  h.Draw()
  ax=h.GetXaxis();
  ax.SetTitle( nameX );
  ay=h.GetYaxis();
  ay.SetTitle( nameY );

  if i>0:
          if (i%2 !=0):
                     gPad.SetLeftMargin(0.1);
                     gPad.SetRightMargin(0.05);

  gPad.SetBottomMargin(0.06);
  gPad.SetTopMargin(0.005);
  if (i>1): gPad.SetBottomMargin(0.20);
  if (i<2): gPad.SetTopMargin(0.05);


  ax.SetLabelSize(0.05)
  ax.SetTitleSize(0.09)
  ax.SetTitleOffset(1.0)
  ay.SetTitleOffset(1.3)

  ay.SetLabelSize(0.05)
  ay.SetTitleSize(0.08)
  gPad.SetTickx()
  gPad.SetTicky()

  ay.Draw("same")
    
  part=particles[i]

  fname="histos/hist_"+part+"_"+str(e)+"gev.root"
  xfile0=TFile( fname )
  hh0=xfile0.Get("scint_vs_EM_E_before_SD")
  hh0.SetDirectory(0)

  hh0.SetMarkerColor( 1 )
  hh0.SetMarkerStyle(20)
  hh0.SetMarkerSize(0.9)
  hh0.SetLineColor( 1 )
  hh0.SetLineWidth(2)

  hh0.Draw("same pe")
  mean0='%.2f'%( hh0.GetMean()  )
  rms0='%.2f'%( hh0.GetRMS()  )

  ######## cherenkov
  hh1=xfile0.Get("cher_vs_EM_E_before_SD")
  hh1.SetDirectory(0)
  xfile0.Close()

  #lab=getFit(hh0, color=1)

  hh1.SetMarkerColor( 4 )
  hh1.SetMarkerStyle(23)
  hh1.SetMarkerSize(0.9)
  hh1.SetLineColor( 4 )
  hh1.SetLineWidth(2)
  hh1.SetLineStyle(1)
  hh1.Draw("same pe")
  mean1='%.2f'%( hh1.GetMean()  )
  rms1='%.2f'%( hh1.GetRMS()  )

  leg2=legends[i]

  pp=part
  if (pp=="mu-"): pp="#mu^{-}"
  if (pp=="gamma"): pp="#gamma"
  if (pp=="kaon-"): pp="K^{-}"
  if (pp=="proton"): pp="p"
  if (pp=="neutron"): pp="n"
  if (pp=="pi0"): pp="#pi^{0}"
  if (pp=="pi-"): pp="#pi^{-}"

  if (i==0):
    leg2.AddEntry(hh0,"Scint","lp")
    leg2.AddEntry(hh1,"Cherenkov","lp")
    leg2.Draw("same");
    #myTextUser(1.0,25000,2,0.07,str(Energy)+" GeV")
    myText(0.71,0.86,2,0.08,str(Energy)+" GeV")

  #myText(0.8,0.9,2,0.12,pp)
 
  if (i==0): myText(0.75,0.76,2,0.09,pp)
  if (i>0): myText(0.75,0.82,2,0.09,pp)
 
  #myTextUser(0.5,25000,2,0.12,pp)
  #myTextUser(0.58,0.12,1,0.07,"<S>="+str(mean0)+"> RMS_{S}="+str(rms0))
  #myTextUser(0.58,0.10,1,0.07,"<C>="+str(mean1)+"> RMS_{C}="+str(rms1))
  #myTextUser(0.1,0.72,2,0.12,pp+" "+str(e)+" GeV")

  f1=functions1[i]
  ## apply fit with 2nd order poly

  chi2=0
  nn=0
  chi2min=10000
  for ii in range(100):
     fitr=hh0.Fit(f1,"SR0E+","",ptmin,ptmax);
     print ("Status=",int(fitr), " is valid=",fitr.IsValid())
     if (fitr.IsValid()==True): 
             chi2=f1.GetChisquare()/f1.GetNDF()
             if chi2<chi2min:
                    chi2min=chi2;
                    if nn>4:
                           break;
                    f1.SetParameter(0,random.randint(0,10000))
                    f1.SetParameter(1,random.randint(0,50))
                    f1.SetParameter(2,random.randint(-2,2))
                    nn=nn+1

  fitr.Print()
  print ("Is valid=",fitr.IsValid())

  par1 = f1.GetParameters()
  a0='%.3f'%( par1[0] )
  a1='%.3f'%( par1[1] )
  a2='%.3f'%( par1[2] )
  s0=a0+"+("+ a1+") x +  ("+a2+") x*x"
  f1.Draw("same")

  par1e = f1.GetParErrors ()
  a0e='%.3f'%( par1e[0] )
  a1e='%.3f'%( par1e[1] )
  a2e='%.3f'%( par1e[2] )
  xfileS.write(pp+" "+a0+" "+a1+" "+a2+ " "+a0e+" "+a1e+" "+a2e+"\n")

  f2=functions2[i]
  nn=0
  chi2min=10000
  for ii in range(100):
     fitr=hh1.Fit(f2,"SR0E+","",ptmin,ptmax);
     print ("Status=",int(fitr), " is valid=",fitr.IsValid())
     if (fitr.IsValid()==True):
             chi2=f2.GetChisquare()/f2.GetNDF()
             if chi2<chi2min:
                    chi2min=chi2;
                    if nn>4:
                           break;
                    f2.SetParameter(0,random.randint(0,10000))
                    f2.SetParameter(1,random.randint(0,50))
                    f2.SetParameter(2,random.randint(-2,2))
                    nn=nn+1

  fitr.Print()
  print ("Is valid=",fitr.IsValid())

  par2 = f2.GetParameters()
  a0='%.3f'%( par2[0] )
  a1='%.3f'%( par2[1] )
  a2='%.3f'%( par2[2] )
  s2=a0+"+("+ a1+") x +  ("+a2+") x*x"
  f2.Draw("same")

  par2e = f2.GetParErrors ()
  a0e='%.3f'%( par2e[0] )
  a1e='%.3f'%( par2e[1] )
  a2e='%.3f'%( par2e[2] )

  xfileC.write(pp+" "+a0+" "+a1+" "+a2+ " "+a0e+" "+a1e+" "+a2e+"\n")


  #myTextUser(0.58,700,1,0.07,"<S>="+s0)
  #myTextUser(0.58,500,1,0.07,"<C>="+s1)

gPad.RedrawAxis()
c1.Update()


if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);

