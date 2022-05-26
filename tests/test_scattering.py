import numpy as np

from atomic_radius.scattering import form_factor


def test_form_factor_with_float_angle():
	assert isinstance(form_factor(20, 1, 1, [1, 1e-15, 1e-15]), np.ndarray)