from atomic_radius.relativity import relativistic_velocity, relativistic_momentum
from atomic_radius.constants import c

def test_relativistic_velocity_for_photon():
	"""
	Test the relativistic velocity for a massless particle
	"""
	assert relativistic_velocity(E=1, m=0) == c

def test_relativistic_momentum_for_photon():
	"""
	Test the relativistic velocity for a massless particle
	"""
	assert relativistic_momentum(E=1, m=0) == 1/c

def test_relativistic_velocity_for_stationary_particle():
	"""
	Test the relativistic velocity for a massless particle
	"""
	assert relativistic_velocity(E=1*c**2, m=1) == 0

def test_relativistic_momentum_for_stationary_particle():
	"""
	Test the relativistic velocity for a massless particle
	"""
	assert relativistic_momentum(E=1*c**2, m=1) == 0
