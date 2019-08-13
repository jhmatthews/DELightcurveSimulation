import matplotlib
matplotlib.use("TkAgg")
import numpy as np 
import matplotlib.pyplot as plt
from constants import *
from scipy.interpolate import interp1d
from rectilinear_losses import *

convert = C / PARSEC / 1e6 
fname_cmb = "../escape/CRPropa3-data/data/Photodisintegration/rate_CMB.txt"
fname_irb = "../escape/CRPropa3-data/data/Photodisintegration/rate_IRB_Gilmore12.txt"
plt.rcParams["text.usetex"] = "True"
def combine(x1,y1,x2,y2, log=True):

	#print (x1.shape, y1.shape, x2.shape, y2.shape)
	if x1.min() > x2.min():
		grid = x1
		y = y2
		x = x2
		ykeep = y1
	else:
		grid = x2
		x = x1
		y = y1
		ykeep = y2

	if log: 
		y = np.log10(y)

	#print (x.shape, y.shape)
	interp_func = interp1d(x, y)
	y_new = 10.0**interp_func(grid) + ykeep
	return grid, y_new

data = np.loadtxt(fname_irb, unpack=True).T
gamma_global = np.linspace(6,14,num=len(data[0][2:]))

data_cmb = np.loadtxt(fname_cmb, unpack=True).T

Z = [1,2,7,26]
N = [0,2,7,30]
label=["H","He", "N", "Fe"]

convert = 3.086e18 * 1e6 / 3e10 / YR / 1e6
#
Zuse = [1,1,1,1]

for ielem in range(len(Z)):
	#print (ielem)
	plt.figure(figsize=(7,5))
	energies, ppLoss, pdLoss, epLoss = calc_losses(units="time", A=Z[ielem]+N[ielem], Z=Z[ielem])
	energies*=1e18
	if N[ielem] > 0: 
		plt.plot(energies, pdLoss, label="Photodisintegration", lw=3, c="C2")
	plt.plot(energies, ppLoss, label="Photopion", lw=3, c="C1")
	plt.plot(energies, epLoss, label="Pair Production", lw=3, c="C0")
	total = 1.0 / ((1.0 / pdLoss) + (1.0 / ppLoss) + (1.0 / epLoss))
	plt.plot(energies, total, label="Pair Production", lw=4, c="k", alpha=0.5)


	array_to_save = np.column_stack( (energies, total, ppLoss, pdLoss, epLoss) )
	np.savetxt("tau_"+label[ielem]+".dat", array_to_save)

	L_scale = 100.0 * 1000.0 * PARSEC
	L_lobe = 100.0 * 1000.0 * PARSEC
	energies = np.logspace(18,21,num=500)
	E =	4.8035e-10
	# B field in Tesla
	B = 1.0 * 1e-6 
	#print (energies)
	#energ = 1e18 * EV2ERGS
	r_g = energies * EV2ERGS / B / 4.8035e-10 / Z[ielem]


	print (r_g / PARSEC, r_g)
	v_d = r_g / L_scale * C / 3
	tau = L_lobe / v_d
	tau = L_lobe * L_lobe / r_g / C * 3

	MYR = 1e6 * YR
	plt.plot(energies/Zuse[ielem], tau/MYR, lw=3, c="k", label=r"$\tau_{\mathrm{esc}}, L=100$kpc")
	lower = np.zeros_like(energies) + 10.0
	#plt.plot(energies, lower, c="C1", alpha=0.8, lw=3)

	from matplotlib.patches import Ellipse, Polygon
	#plt.gca().add_patch(Polygon([[1e18, 0], [1e20, 0], [1e20, 10], [1e18, 10]], closed=True,
	#                      fill=True, alpha=0.5, color="C1"))
	plt.plot(energies, 22.0 * energies / 1e19 * 1e3 * 1e-2/Z[ielem], ls="--", c="k", label=r"$\tau_{F2}$")
	plt.plot(energies, np.ones_like(energies)*500.0, ls="-.", label=r"$\tau_{\mathrm{buoy}}$", c="k")
	plt.plot(energies, np.ones_like(energies)*24.0, ls=":", label=r"$\tau_{\mathrm{sync}}$", c="k")

	plt.loglog()
	plt.legend()
	#plt.text(5e19,2e3,"GZK", fontsize=16)
	#plt.text(1.5e18,5.5,"Reservoir", fontsize=12)
	plt.xlabel("$E~(eV)$", fontsize=18)
	plt.ylabel(r"$\tau~(\mathrm{Myr})$", fontsize=18)

	plt.xlim(1e18,1e21)
	plt.ylim(1,1e4)
	plt.legend(loc=3)
	plt.title(label[ielem], fontsize=18)
	plt.savefig("new_tau_"+label[ielem]+".png", dpi=300)
	plt.clf()
	# plt.show()
