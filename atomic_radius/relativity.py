import numpy as np

from .constants import *

def relativistic_velocity(E, m):
	"""
	Return the relativistic velocity of a particule with energy E and mass m
	"""
	return c * np.sqrt(1 - m**2 * c**4 / E**2)

def relativistic_momentum(E, m):
	"""
	Return the relativistic momentum of a particule with energy E and mass m
	"""
	return np.sqrt(E**2 / c**2 - m**2 * c**2)