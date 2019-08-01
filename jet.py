from DELCgen import *
import scipy.stats as st
import matplotlib.pyplot as plt 
from constants import *

# File Route
route = "./"
datfile = "NGC4051.dat"

# Bending power law params
pl_index = 1
A,v_bend,a_low,a_high,c = 1, 300, pl_index, pl_index, 1
# Probability density function params
kappa,theta,lnmu,lnsig,weight = 5.67, 5.96, 2.14, 0.31,0.82
# Simulation params
# let's do everything in units of kyr
# run for 100 Myr (1e5 kyr) in bins of 100 kyr
lognorm_params = (1.5,0,np.exp(1.5))
RedNoiseL,RandomSeed,aliasTbin, tbin, Age = 100,12,100,100,100
N = Age / tbin

#delc = Simulate_DE_Lightcurve(BendingPL, (A,v_bend,a_low,a_high,c),
#				st.lognorm,(0.3, 0.0, 1e43),
#                                  tbin = 1, LClength = 1000)

#lc = Simulate_TK_Lightcurve(BendingPL, (A,v_bend,a_low,a_high,c),
#                                RedNoiseL=RedNoiseL,aliasTbin=aliasTbin,randomSeed=RandomSeed,length=Age, tbin=tbin,
#                                mean=1e43, std=1e44)

lc = Simulate_DE_Lightcurve(BendingPL, (A,v_bend,a_low,a_high,c),st.lognorm,lognorm_params,
                                RedNoiseL=RedNoiseL,aliasTbin=aliasTbin,randomSeed=RandomSeed,LClength=Age, tbin=tbin)





#Comparison_Plots([lc,lc])
#lc.Plot_Lightcurve()
lc.Save_Lightcurve('lightcurve.dat')
flux_scale = 1e43
flux = lc.flux * flux_scale 
#flux = lc.flux
# time is in kyr, so convert to Myr

#plt.clf()


def power_threshold(rigidity, v_over_c=0.1, eta=0.1):
	power = (rigidity / 1e19)**2 * (0.1/v_over_c) * 1e44
	return power 

def max_energy(power, v_over_c=0.1, eta=0.1):
	return 1e19*np.sqrt( (power/1e44) * (v_over_c/0.1))

energies = np.logspace(12,21,num=1000)
powers = power_threshold(energies, v_over_c=0.5)
fraction = np.zeros_like(energies)
for_cdf = powers / flux_scale

for i, energy in enumerate(energies):
	fraction[i] = np.sum(flux > powers[i]) / (1.0*N)


BETA = 2
ncr = np.zeros_like(energies)
escape_time = (1e19 / energies) * 1e7 * YR 
frac_elem = [1.0,0.1,1e-4,3.16e-05]
z_elem = [1,2,7,26]
#flux[lc.time <= 10000] = 1e45
#flux[lc.time > 10000] = 1e30

for i in range(len(flux)):
	#for j in range(len(f_elem)):
	Emax = max_energy(flux[i], v_over_c=0.5)

	print (Emax)
	# normalisation of cosmic ray distribution
	Q0 = 0.1 * flux[i] / (np.log(Emax/1e9))

	delta_t = tbin * 1000.0 * YR

	#print (delta_t, escape_time)

	add_on = Q0 * ((energies / 1e9)**-BETA)
	add_on[(energies>=Emax)] = 0.0
	escaping = ncr / escape_time
	change = (add_on - (escaping))*delta_t

	# if i > 5:
	# 	select = (np.fabs(change) < ncr)
	# 	ncr[select] += change[select]
	# 	ncr[~select] = 0.0
	# else:
	ncr += change
	ncr[(ncr<0)] = 0

	#if plot_all:
	plt.clf()
	plt.figure()
	plt.subplot(211)
	plt.plot(lc.time[:i+1] / 1000.0, flux[:i+1])
	plt.semilogy()
	plt.xlabel("t (Myr)")
	plt.ylabel("Q_j (erg/s)")

	plt.subplot(223)
	#plt.xlabel("E (eV)")
	#plt.ylabel("n(E)")

	#lt.plot(energies, ncr)
	for j,frac in enumerate(frac_elem):
		energy = energies * z_elem[j]
		plt.plot(energy, energy*energy*(frac*ncr))
	plt.loglog()
	plt.xlim(1e17,1e21)
	#plt.ylim(1e13,1e19)

	plt.subplot(224)
	for j,frac in enumerate(frac_elem):
		energy = energies * z_elem[j]
		plt.plot(energy, energy*energy*frac*escaping*delta_t)

	# plt.plot(energies, 1e15*((energies / 1e9)**-2), c="k", ls="--")
	plt.xlim(1e17,1e21)
	#plt.ylim(0.1,1e5)
	#plt.plot(energies, 1e43*((energies / 1e9)**-3), c="k", ls=":")
	plt.loglog()
	plt.savefig("spectra{:03d}.png".format(i))	
	plt.close("all")

#import sys
#sys.exit()


plt.subplot(211)
plt.plot(lc.time[:i+1] / 1000.0, flux[:i+1])
plt.semilogy()
plt.xlabel("t (Myr)")
plt.ylabel("Q_j (erg/s)")

plt.subplot(234)
#plt.plot(energies, fraction*100.0)
powers2 = np.arange(40,47,0.1)
for_pdf = powers2 / flux_scale
plt.hist(np.log10(flux), log=True, density=True, bins=powers2)
plt.plot(np.log10(powers2), st.lognorm.pdf(for_pdf, lognorm_params[0], lognorm_params[1], lognorm_params[2]), c="r")
#plt.semilogx()
plt.xlabel("$E_{max}$ (eV)")
plt.ylabel("PDF")
plt.xlim(40,47)

plt.subplot(235)
plt.plot(energies, fraction*100.0)
plt.plot(energies, 1e4*st.lognorm.sf(for_cdf, lognorm_params[0], lognorm_params[1], lognorm_params[2]))
plt.loglog()
plt.xlabel("$E_{max}$ (eV)")
plt.ylabel("Survival Function")

plt.subplot(236)
plt.xlabel("E (eV)")
plt.ylabel("n(E)")

#lt.plot(energies, ncr)
plt.plot(energies, ncr/escape_time)
plt.plot(energies, 1e15*((energies / 1e9)**-2), c="k", ls="--")
plt.xlim(1e15,1e21)
#plt.plot(energies, 1e43*((energies / 1e9)**-3), c="k", ls=":")
plt.loglog()
plt.show()	






#tklc = Simulate_TK_Lightcurve(BendingPL, (1.0,300,2.1,2.3,0.1),
#				st.lognorm,(0.3, 0.0, 7.4),
#                                  tbin = 1, LClength = 1000, mean = 1, std = 0.1)
#--------- Commands ---------------                            

#lc.Plot_Lightcurve()

#lc.Plot_Periodogram()

# # Save lightcurve and Periodogram as text files
# delc.Save_Lightcurve('lightcurve.dat')
# delc.Save_Periodogram('periodogram.dat')