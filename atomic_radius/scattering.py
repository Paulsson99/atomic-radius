import numpy as np
from scipy.integrate import quad

from atomic_radius.nucleus import charge_density, integrate_charge
from atomic_radius.relativity import relativistic_velocity, relativistic_momentum
import atomic_radius.constants as const


def rutherford(Z, E, theta):
	"""
	Calculate the Rutherford cross-section

	dsigma/dOmega = Z**2 * alpha**2 * (h_bar*c)**2 / (4 * E**2 * sin(theta/2)**4)
	
	args:
		Z: Charge of the projectile (in units of e)
		E: Energy of the projectile
		theta: Scattering angle of the projectile
	"""
	beta = relativistic_velocity(E, const.m_e) / const.c
	return (Z * const.alpha * const.h_bar * const.c / (2 * beta * E * np.sin(theta/2)**2) )**2

def mott(Z, E, theta):
	"""
	Calculate the Mott cross-section
	
	dsigma/dOmega = rutherford(Z, E, theta) * [1 - beta**2 * sin(theta/2)**2)]

	args:
		Z: Charge of the projectile (in units of e)
		E: Energy of the projectile
		theta: Scattering angle of the projectile
	"""
	beta = relativistic_velocity(E, const.m_e) / const.c
	return rutherford(Z, E, theta) * (1 - beta**2 * np.sin(theta/2)**2)

def experiment(Z, E, theta, X):
	"""
	Calculate the experimental cross-section

	dsigma/dOmega = mott(Z, E, theta) * |F(q**2)|**2

	args:
		Z: Charge of the projectile (in units of e)
		E: Energy of the projectile
		theta: Scattering angle of the projectile
		X: Constants for the charge_density
	"""
	F = form_factor(Z, E, theta, X)
	return mott(Z, E, theta) * np.abs(F)**2

def form_factor(Z, E, theta, X):
	"""
	Calculate the form factor

	args:
		Z: Charge of the projectile (in units of e)
		E: Energy of the projectile
		theta: Scattering angle of the projectile
		X: Constants for the charge_density
	"""
	p = relativistic_momentum(E, const.m_e)
	q = np.sqrt(4 * p**2 * np.sin(theta/2)**2)

	func = lambda r: np.sin(q * r / const.h_bar)
	charge_integral = integrate_charge(X, power=1, func=func)

	return charge_integral * 4 * np.pi * const.h_bar / (Z * const.e * q)



