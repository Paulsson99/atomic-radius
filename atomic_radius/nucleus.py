import numpy as np
from scipy.integrate import quad

import atomic_radius.constants as const
import atomic_radius.config as conf


def charge_density(r, X):
	"""
	Charge density is assumed to be shericaly symetric

	rho(r, rho_0, a, b) = rho_0 / {1 + exp([r-a]/b)}

	X:
		rho_0: [C/m**3]
		a: [m]
		b: [m]

	returns:
		rho: [C/m**3]
	"""
	return X[0] / ( 1 + np.exp( (r - X[1]) / X[2]) )

def rms_radius(Z, X):
	"""
	Calculate the rms radius of the nucleus

	args:
		Z: Number of protons in the nucleus
		X: Constants for the charge_density
	
	Z: [-]
	X:
		rho_0: [C/m**3]
		a: [m]
		b: [m]

	returns:
		r: [m]
	"""
	charge_integral = integrate_charge(X, power=4)
	factor = 4 * np.pi / (Z * const.e)
	return np.sqrt(factor * charge_integral)

def protons(X):
	"""
	Calculate the total charge of the nucleus

	X:
		rho_0: [C/m**3]
		a: [m]
		b: [m]

	returns:
		charge: [C]

	"""
	return 4 * np.pi / const.e * integrate_charge(X, power=2)

def integrate_charge(X, power=0, func=lambda r: 1):
	"""
	Integrate the function r^power * charge_density(X)
	from 0 to inf
	"""
	# Charge density must be small outside this value
	inf = 5*X[2] + X[1]
	r = np.reshape(np.linspace(0, inf, conf.N), (conf.N, 1))

	integral_func = r**power * charge_density(r, X) * func(r)
	return np.trapz(integral_func, x=r, axis=0)

def nuclear_units_to_SI(X):
	return np.array([X[0] * const.e / const.fm**3, X[1] * const.fm, X[2] * const.fm])
