from ROOT import *

ROOT.gSystem.Load('roostats_cl95_C.so')


ilum=200 #  - Nominal integrated luminosity (pb^-1)         \n"
slum=0.01   #  - Absolute error on the integrated luminosity   \n"
eff=0.8  #  - Nominal value of the efficiency times         \n"
         #   acceptance (in range 0 to 1)                  \n"
seff=0.1 #  - Absolute error on the efficiency times        \n"
         #    acceptance                                    \n"
bck=100 #  - Nominal value of the background estimate      \n"
sbck=20  #  - Absolute error on the background              \n"
n=100   #  - Number of observed events (not used for the   \n"
         #   expected limit)                               \n"
ntoys=500 #    - Number of pseudoexperiments to perform for    \n"
           #      expected limit calculation)                   \n"
gauss=False  #   - if true, use Gaussian statistics for signal   \n"
            #     instead of Poisson; automatically false       \n"
            #                       for n = 0.                                    \n"
#                     Always false for expected limit calculations  \n"
nuisanceModel=0 #  - distribution function used in integration over\n"
               #        nuisance parameters:                          \n"
               #        0 - Gaussian (default), 1 - lognormal,        \n"
               #       2 - gamma;                                    \n"
#                       (automatically 0 when gauss == true)          \n"
method="mcmc"  #      - method of statistical inference:              \n"
               #         \"bayesian\"  - Bayesian with numeric         \n"
"""
                                       integration (default),        \n"
                       \"mcmc\"      - another implementation of     \n"
                                       Bayesian, not optimized,      \n"
                                       to be used for cross checks   \n"
                                       only!                         \n"
                       \"cls\"       - CLs observed limit. We suggest\n"
                                       using the dedicated interface \n"
                                       roostats_cls() instead        \n"
                       \"fc\"        - Feldman Cousins with numeric  \n"
                                     integration,                    \n"
                       \"workspace\" - only create workspace and save\n"
                                     to file, no interval calculation\n"
"""
plotFileName="plot.pdf" #   - file name for the control plot to be created  \n"

seed=0          # seed for random number generation,            \n"
                # specify 0 for unique irreproducible seed 


limit=roostats_cl95(ilum,slum,eff,seff,bck,sbck, n,gauss, nuisanceModel, method, plotFileName, seed);
#limit = roostats_limit(ilum, slum, eff, seff, bck, sbck, n, gauss, nuisanceModel, method, plotFileName, seed); 

limit = roostats_clm(ilum, slum, eff, seff, bck, sbck, ntoys, nuisanceModel, method, seed);
#print "expected=",expected_limit


#average_limit = roostats_cla(ilum, slum, eff, seff, bck, sbck, nuisanceModel, method, seed);
#print "print=", average_limit

obs_limit = limit.GetObservedLimit()
exp_limit = limit.GetExpectedLimit()
exp_up    = limit.GetOneSigmaHighRange()
exp_down  = limit.GetOneSigmaLowRange()
exp_2up   = limit.GetTwoSigmaHighRange()
exp_2down = limit.GetTwoSigmaLowRange()

print "Observed 95% = ",obs_limit
print "Expected 95% (median)  = ",exp_limit
print "Expected +1 sigma = ",exp_up
print "Expected -1 sigma = ",exp_down
print "Expected +2 sigma = ",exp_2up
print "Expected -2 sigma = ",exp_2down 

 
