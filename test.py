from DELCgen import *
import scipy.stats as st

# File Route
route = "./"
datfile = "NGC4051.dat"

# Bending power law params
A,v_bend,a_low,a_high,c = 0.03, 2.3e-4, 1.1, 2.2, 0 
# Probability density function params
kappa,theta,lnmu,lnsig,weight = 5.67, 5.96, 2.14, 0.31,0.82
# Simulation params
RedNoiseL,RandomSeed,aliasTbin, tbin = 100,12,1,100 

#--------- Commands ---------------

# load data lightcurve
datalc = Load_Lightcurve(route+datfile,tbin)

# plot the data lightcurve and its PDF and PSD
#datalc.Plot_Lightcurve()

# estimate underlying variance od data light curve
datalc.STD_Estimate()


# simulate artificial light curve with Emmanoulopoulos method, using the PSD and PDF of the data
delc = datalc.Simulate_DE_Lightcurve() # defaults to bending PL and mix of gamma and lognormal dist.

# simulate artificial light curve with Timmer & Koenig method
tklc = Simulate_TK_Lightcurve(datalc,BendingPL, (A,v_bend,a_low,a_high,c),
                                RedNoiseL,aliasTbin,RandomSeed)

# simulate artificial light curve with Emmanoulopoulos method, scipy distribution
delc2 = Simulate_DE_Lightcurve(datalc,BendingPL, (A,v_bend,a_low,a_high,c),
                                ([st.gamma,st.lognorm],[[kappa,0, theta],\
                                    [lnsig,0, np.exp(lnmu)]],[weight,1-weight]))

# simulate artificial light curve with Emmanoulopoulos method, custom distribution
delc3 = Simulate_DE_Lightcurve(datalc,BendingPL, (A,v_bend,a_low,a_high,c),
                                ([[Gamma,LogNormal],[[kappa, theta],\
                                  [lnmu, lnsig]],[weight,1-weight]]),MixtureDist)                                

# plot lightcurves and their PSDs ands PDFs for comparison
Comparison_Plots([datalc,tklc,delc])

# Save lightcurve and Periodogram as text files
delc.Save_Lightcurve('lightcurve.dat')
delc.Save_Periodogram('periodogram.dat')