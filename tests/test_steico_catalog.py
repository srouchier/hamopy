from contextlib import redirect_stdout
from io import StringIO
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from xml.etree import ElementTree as ET

import numpy as np
from numpy.testing import assert_allclose

from hamopy import ham_library as ham
from hamopy.algorithm import calcul
from hamopy.classes import Boundary, Material, Mesh, Time
from hamopy.materials import steico


FIXTURE = Path(__file__).parent / 'fixtures' / 'wufi_synthetic.xml'
DATA = Path(__file__).resolve().parents[1] / 'data' / 'wufi' / 'steico'


class SteicoCatalogTests(unittest.TestCase):
    def test_catalog_covers_the_manifest(self):
        manifest = json.loads((DATA / 'manifest.json').read_text())
        self.assertEqual(len(steico.MATERIAL_FILES), 19)
        self.assertEqual(set(steico.MATERIAL_FILES.values()),
                         {Path(item['xml']).name for item in manifest['materials']})
        self.assertNotIn('steico_internal', steico.MATERIAL_FILES)

    def test_catalog_and_named_imports_return_independent_materials(self):
        with tempfile.TemporaryDirectory() as directory:
            for filename in steico.MATERIAL_FILES.values():
                (Path(directory) / filename).write_bytes(FIXTURE.read_bytes())
            materials = steico.load_materials(directory)
            self.assertEqual(len(materials), 19)
            self.assertTrue(all(isinstance(m, Material) for m in materials.values()))
            with patch.dict(os.environ, HAMOPY_STEICO_XML_DIR=directory):
                from hamopy.materials.steico import steico_therm
                self.assertIsNot(steico_therm, steico.steico_therm)
                steico_therm.set_density(123)
                self.assertEqual(steico.steico_therm.rho, 200.)
                suction = steico.load_material('steico_therm', liquid_transport='suction')
                self.assertEqual(suction.liquid_transport, 'suction')

    def test_missing_files_explain_how_to_obtain_them(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, 'data/wufi/steico/README.md'):
                steico.load_materials(directory)
            with self.assertRaises(KeyError):
                steico.load_material('steico_internal', directory)


def run_moisture_step(material, initial_humidity, ambient_humidity):
    mesh = Mesh([material], [.06], [4])
    boundary = [Boundary('Fourier', T=293.15, HR=ambient_humidity)] * 2
    time = Time('variable', delta_t=30., t_max=3600., delta_max=300.)
    output = StringIO()
    with redirect_stdout(output):
        result = calcul(mesh, boundary, {'T': 293.15, 'HR': initial_humidity}, time)
    return result, output.getvalue()


class CoupledSimulationTests(unittest.TestCase):
    def test_synthetic_material_wets_and_dries_in_coupled_solver(self):
        from hamopy.materials.wufi import load_wufi_material
        material = load_wufi_material(FIXTURE)
        for initial, ambient in ((.85, .95), (.95, .85)):
            with self.subTest(initial=initial):
                result, output = run_moisture_step(material, initial, ambient)
                self.assertEqual(output, '')
                self.assertEqual(result['t'][-1], 3600.)
                self.assertTrue(np.all(np.isfinite(result['U'])))
                self.assertTrue(np.all((result['HR'] > 0) & (result['HR'] < 1)))
                change = material.w(result['PC'][-1]).mean() - material.w(result['PC'][0]).mean()
                self.assertGreater(change * (ambient - initial), 0)


@unittest.skipUnless(all((DATA / 'xml' / name).is_file()
                        for name in steico.MATERIAL_FILES.values()),
                     'Manufacturer XMLs are optional local downloads.')
class ManufacturerIntegrationTests(unittest.TestCase):
    def test_all_materials_preserve_source_tables_and_have_finite_properties(self):
        for key, material in steico.load_materials(DATA / 'xml').items():
            with self.subTest(material=key):
                root = ET.parse(material.source['path'])
                items = root.findall('Material/Functions/Function[@id="1"]/Items/Item')
                points = np.array([[float(p.get(v).replace(',', '.')) for v in ('x', 'y')]
                                   for p in items])
                assert_allclose(material.w(ham.p_c(points[1:, 0])), points[1:, 1])
                pressure = ham.p_c(np.linspace(.001, 1, 1000))
                self.assertTrue(np.all(np.diff(material.w(pressure)) >= 0))
                for function in (material.w, material.c_w, material.k_l,
                                 material.conduc, material.delta_p):
                    values = function(pressure)
                    self.assertTrue(np.all(np.isfinite(values)))
                    self.assertTrue(np.all(values >= 0))
                assert_allclose(material.k_l(0.), material.liquid_diffusivity(material.w_sat)
                                * material.c_w(0.))

    def test_all_materials_complete_wetting_and_drying_simulations(self):
        for key, material in steico.load_materials(DATA / 'xml').items():
            for initial, ambient in ((.85, .95), (.95, .85)):
                with self.subTest(material=key, initial=initial):
                    result, output = run_moisture_step(material, initial, ambient)
                    self.assertEqual(output, '')
                    self.assertEqual(result['t'][-1], 3600.)
                    self.assertTrue(np.all(np.isfinite(result['U'])))
                    self.assertTrue(np.all((result['HR'] > 0) & (result['HR'] < 1)))
                    change = material.w(result['PC'][-1]).mean() - material.w(result['PC'][0]).mean()
                    self.assertGreater(change * (ambient - initial), 0)


if __name__ == '__main__':
    unittest.main()
