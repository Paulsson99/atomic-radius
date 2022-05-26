from .constants import *

Z = 20	# Number of protons in the nucleus
E = 250	# [Mev] electron energy

scattering_data_file = '/Users/Paulsson/Documents/Private/Simon/Skola/Högskola/Kurser/Åk3/LP4/Subatomic/Datoruppgift/atomic-radius/atomic_radius/scattering_data.csv'

# Other
N = 1000	# Number of divisions in the integral

# Recalculate to SI-units
E = E * eV2J * 1e6