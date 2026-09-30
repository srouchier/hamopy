from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
from xml.etree import ElementTree as ET

import numpy as np
from numpy.testing import assert_allclose

from hamopy import ham_library as ham
from hamopy.classes import Material, Mesh
from hamopy.materials.wufi import load_wufi_material


FIXTURE = Path(__file__).parent / 'fixtures' / 'wufi_synthetic.xml'


class TabulatedSorptionTests(unittest.TestCase):
    def setUp(self):
        self.material = Material()
        self.material.set_isotherm('tabulated', HR=[0., .5, .8, 1.],
                                  W=[0., 5., 10., 100.])

    def test_preserves_points_without_overshoot(self):
        for temperature in (273.15, 293.15, 313.15):
            pressure = ham.p_c(np.array([.5, .8, 1.]), temperature)
            assert_allclose(self.material.w(pressure, temperature), [5., 10., 100.])
        pressure = ham.p_c(np.linspace(.0001, 1, 1000))
        moisture = self.material.w(pressure)
        self.assertTrue(np.all(np.diff(moisture) >= 0))
        self.assertTrue(np.all((moisture >= 0) & (moisture <= 100)))
        self.assertEqual(self.material.w(-np.inf), 0.)
        self.assertEqual(self.material.c_w(-np.inf), 0.)
        assert_allclose(self.material.w([0., 1., 1e12]), 100.)
        assert_allclose(self.material.c_w([1., 1e12]), 0.)

    def test_capacity_matches_pressure_derivative(self):
        pressure = ham.p_c(np.array([[.2, .6], [.85, .999]]))
        difference = (self.material.w(pressure + 1) - self.material.w(pressure - 1)) / 2
        assert_allclose(self.material.c_w(pressure), difference, rtol=1e-6)
        self.assertTrue(np.all(difference > 0))

    def test_invalid_tables_are_rejected(self):
        for humidity, moisture in [
            ([0, .5, .5, 1], [0, 5, 10, 100]),
            ([0, .5, .8, 1], [0, 10, 5, 100]),
            ([0, .5, .8, 1], [0, 5, np.nan, 100]),
            ([.1, .5, .8, 1], [0, 5, 10, 100]),
        ]:
            with self.subTest(humidity=humidity, moisture=moisture):
                with self.assertRaises(ValueError):
                    self.material.set_isotherm('tabulated', HR=humidity, W=moisture)


class WufiMaterialTests(unittest.TestCase):
    def setUp(self):
        self.material = load_wufi_material(FIXTURE)

    def test_converts_decimal_commas_and_ten_degree_conductivity_reference(self):
        self.assertEqual(self.material.rho * self.material.cp(), 300000.)
        self.assertAlmostEqual(self.material.conduc(T=283.15), .05)
        self.assertAlmostEqual(self.material.conduc(T=303.15), .07)
        self.assertAlmostEqual(self.material.conduc(ham.p_c(.8), 283.15), .055)
        assert_allclose(self.material.delta_p(ham.p_c(np.array([.1, .9]))),
                        26.1e-6 / (4 * ham.Rv * 293.15))
        self.assertIn('Artificial material', self.material.source['info'])

    def test_generated_liquid_curves_use_reference_and_saturation_points(self):
        suction = self.material.liquid_diffusivity(np.array([0, 9.99, 10, 55, 100]),
                                                   mode='suction')
        redistribution = self.material.liquid_diffusivity(
            np.array([0, 9.99, 10, 55, 100]), mode='redistribution')
        # Published generation formula at w80=10, wf=100, A=0.1.
        start = 7.581996796881743e-9
        assert_allclose(suction, [0, 0, start, np.sqrt(start * 3.8e-6), 3.8e-6])
        assert_allclose(redistribution,
                        [0, 0, start, np.sqrt(start * 3.8e-7), 3.8e-7])

    def test_pressure_permeability_preserves_liquid_flux(self):
        pressure = ham.p_c(np.array([[.9, .95], [.98, .999]]))
        moisture = self.material.w(pressure)
        pressure_gradient = 25000.
        difference = (self.material.w(pressure + 1) - self.material.w(pressure - 1)) / 2
        flux_from_moisture = self.material.liquid_diffusivity(moisture) * difference * pressure_gradient
        assert_allclose(self.material.k_l(pressure) * pressure_gradient,
                        flux_from_moisture, rtol=1e-6)
        self.assertTrue(np.all(flux_from_moisture > 0))

    def test_transport_mode_is_explicit_and_copyable(self):
        suction = load_wufi_material(FIXTURE, liquid_transport='suction')
        self.assertAlmostEqual(suction.k_l(0.) / self.material.k_l(0.), 10.)
        copied = deepcopy(suction)
        self.assertEqual(copied.k_l(ham.p_c(.95)), suction.k_l(ham.p_c(.95)))
        with self.assertRaises(ValueError):
            load_wufi_material(FIXTURE, liquid_transport='automatic')

    def test_user_can_replace_liquid_permeability(self):
        self.material.set_perm_liquid('exp', a=[-30.])
        self.assertAlmostEqual(self.material.k_l(ham.p_c(.9)), np.exp(-30.), places=20)

    def test_rejects_incompatible_units_and_incomplete_models(self):
        changes = [
            ('.', 'UnitSystem', '1'),
            ('Material/Scalars/Scalar[@id="10"]', 'unit', 'kg/m²√h'),
            ('Material/Scalars/Scalar[@id="11"]', 'value', '0'),
            ('Material/Scalars/Scalar[@id="4"]', 'value', 'nan'),
            ('Material/Scalars/Scalar[@id="6"]', 'value', '11'),
            ('Material/Functions/Function[@id="2"]', 'version', '2'),
            ('Material/Functions/Function[@id="1"]/Items', 'Count', '3'),
            ('Material/Functions/Function[@id="5"]/Items/Item', 'y', '6'),
        ]
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'material.xml'
            for xpath, key, value in changes:
                root = ET.parse(FIXTURE).getroot()
                root.find(xpath).set(key, value)
                ET.ElementTree(root).write(target, encoding='utf-8')
                with self.subTest(xpath=xpath, key=key):
                    with self.assertRaises(ValueError):
                        load_wufi_material(target)
            root = ET.parse(FIXTURE).getroot()
            functions = root.find('Material/Functions')
            functions.remove(functions.find('Function[@id="1"]'))
            ET.ElementTree(root).write(target, encoding='utf-8')
            with self.assertRaises(ValueError):
                load_wufi_material(target)
            for identifier, field, value in [('2', 'Generated', 'FALSE'),
                                             ('3', 'GenerationMode', '1')]:
                root = ET.parse(FIXTURE).getroot()
                root.find('Material/Functions/Function[@id="{}"]/{}'.format(
                    identifier, field)).text = value
                ET.ElementTree(root).write(target, encoding='utf-8')
                with self.subTest(function=identifier, field=field):
                    with self.assertRaises(ValueError):
                        load_wufi_material(target)


class ThermalAssemblyTests(unittest.TestCase):
    def test_temperature_is_not_passed_as_capillary_pressure(self):
        material = Material()
        material.set_conduc(lambda_0=.02, lambda_t=.001)
        mesh = Mesh([material], [.06], [2])
        _, conductivity = mesh.system_matrices_thermo(np.full(mesh.nbr_nodes, 283.15))
        flux = conductivity.tocsr().dot(mesh.x)
        assert_allclose(flux[[0, -1]], [-.03, .03])


if __name__ == '__main__':
    unittest.main()
