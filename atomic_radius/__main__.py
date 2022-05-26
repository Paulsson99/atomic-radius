import numpy as np
import pandas as pd
from scipy.optimize import minimize, fsolve, leastsq
from scipy.integrate import quad
import matplotlib.pyplot as plt

from atomic_radius.scattering import experiment
from atomic_radius.nucleus import charge_density, rms_radius, protons, nuclear_units_to_SI
import atomic_radius.config as conf
import atomic_radius.constants as const


def quadratic_sum(X, theta, cross_section, error):
	total_sum = 0

	X = nuclear_units_to_SI(X)
	excpected_cross_section = experiment(conf.Z, conf.E, theta, X) * const.m_sq2barn * 1e3
	return np.sum(((excpected_cross_section - cross_section) / error)**2)
	# for angle, cs, err in zip(theta, cross_section, error):
	# 	# Calculate excpected cross-section in mb/sr
	# 	excpected_cross_section = experiment(conf.Z, conf.E, angle, X) * const.m_sq2barn * 1e3
	# 	total_sum += ((excpected_cross_section - cs) / err)**2
	# return total_sum

def charge_condition(X):
	"""
	The total charge must sum to Z
	"""
	X = nuclear_units_to_SI(X)
	return conf.Z - protons(X)

def main():
	# Open scattering data
	scattering_data = pd.read_csv(conf.scattering_data_file)
	theta = np.deg2rad(scattering_data['theta (degree)'].to_numpy())
	cross_section = scattering_data['cross-section (mb/sr)'].to_numpy()
	error = scattering_data['error (mb/sr)'].to_numpy()

	# Fit rho_0, a and b to data
	constraint = {
		'type': 'eq',
		'fun': charge_condition
	}
	bounds = ((0.06, 0.08), (3, 5), (.3, .7))

	res = minimize(
		quadratic_sum, 
		x0=[0.07, 4, 5], 
		args=(theta, cross_section, error), 
		constraints=constraint, 
		bounds=bounds
	)

	print(quadratic_sum(res.x, theta, cross_section, error))
	print(charge_condition(res.x))
	print(protons(nuclear_units_to_SI(res.x)))

	print(res.x)

	rms_r = rms_radius(conf.Z, nuclear_units_to_SI(res.x))

	print(rms_r)

	# Visualize
	fig1 = plt.figure(1)
	fig2 = plt.figure(2)

	ax1 = fig1.add_subplot(111)
	ax1.scatter(np.rad2deg(theta), cross_section)
	# ax1.errorbar(theta, cross_section, yerr=error, fmt='none', ecolor='red')
	lin_theta = np.linspace(np.min(theta), np.max(theta), 200)
	ax1.plot(
		np.rad2deg(lin_theta), 
		experiment(conf.Z, conf.E, lin_theta, nuclear_units_to_SI(res.x)) * const.m_sq2barn * 1e3,
		color='black'
	)
	ax1.set_yscale('log')
	ax1.grid(True)
	ax1.set_ylabel(r"$d\sigma/d\Omega\quad[mb/sr]$")
	ax1.set_xlabel(r"$\theta\quad[grader]$")
	ax1.legend(["Experimentell data", r"Teori baserad på $\mathcal{X}^2(\mathbf{X})$"])

	ax2 = fig2.add_subplot(111)
	r = np.linspace(0, 8, 200)
	rho = charge_density(r*const.fm, nuclear_units_to_SI(res.x)) * const.fm**3 / const.e
	ax2.plot(r, rho, color='black')
	ax2.axvline(x=rms_r / const.fm, ymin=0, ymax=1, linestyle='--', color='black')
	ax2.fill_between(r, rho, color='red', alpha=0.5)
	ax2.grid(False)
	ax2.set_ylabel(r"$\rho_\mathrm{ch}(r, \mathbf{X})\quad[efm^{-3}]$")
	ax2.set_xlabel(r"$r\quad[fm]$")
	ax2.set_xlim(0, 8)
	ax2.set_ylim(0, 0.08)

	ax2.legend([r"Protondensitet $\rho_\mathrm{ch}(r, \mathbf{X})$", r"rms-laddningsradie $\sqrt{\left\langle r^2\right\rangle}$"])
	
	plt.show()

if __name__ == '__main__':
	main()
