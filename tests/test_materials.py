import unittest

import numpy as np
from numpy.testing import assert_allclose

from hamopy import ham_library as ham
from hamopy.classes import Material, Mesh
from hamopy.materials import standard


class MaterialConductivityTests(unittest.TestCase):
    def test_constant_conductivity_does_not_require_sorption_data(self):
        material = Material()
        material.set_conduc(lambda_0=0.038)

        self.assertAlmostEqual(material.conduc(), 0.038)
        pressure = ham.p_c(np.array([[0.25, 0.5], [0.75, 0.95]]))
        conductivity = material.conduc(pressure)
        self.assertEqual(conductivity.shape, pressure.shape)
        assert_allclose(conductivity, 0.038)

    def test_moisture_dependent_conductivity_still_uses_sorption(self):
        material = standard.concrete
        pressure = ham.p_c(np.array([0.25, 0.75, 0.95]))
        expected = 1.75 + material.w(pressure) / 1000 * 4.5
        assert_allclose(material.conduc(pressure), expected)
        self.assertGreater(expected[-1], expected[0])


class SteicoInternalTests(unittest.TestCase):
    def test_heat_only_mesh_has_published_capacity_and_heat_flux(self):
        # A 60 mm board stores rho * cp * thickness = 20160 J/(m2.K).
        # A 20 K difference drives lambda * 20 / 0.06 = 12.6667 W/m2.
        mesh = Mesh([standard.steico_internal], [0.06], [2])
        temperature = 273.15 + 20 * mesh.x / 0.06
        capacity, conductivity = mesh.system_matrices_thermo(temperature)

        self.assertAlmostEqual(float(capacity.sum()), 20160.)
        flux = conductivity.tocsr().dot(temperature)
        assert_allclose(flux[[0, -1]], [-12.6666666667, 12.6666666667])
        assert_allclose(flux[1:-1], 0., atol=1e-10)

    def test_published_mu_is_constant_across_humidity(self):
        pressure = np.r_[-np.inf, ham.p_c(np.array([0.25, 0.5, 0.95, 1.0]))]
        for temperature in (273.15, 293.15, 313.15):
            permeability = standard.steico_internal.delta_p(pressure, temperature)
            air_permeability = 26.1e-6 / (ham.Rv * temperature)
            # The datasheet gives sd = 0.30 m for the 60 mm board.
            assert_allclose(air_permeability / permeability * 0.06, 0.30)


if __name__ == '__main__':
    unittest.main()
