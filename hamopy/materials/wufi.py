"""Import the SI WUFI XML profile used by STEICO's published material files.

This is a deliberately limited importer, not a reader for arbitrary WUFI
databases. Unsupported units and function definitions raise ValueError.
Manufacturer XML files are read locally and are not distributed with Hamopy.
"""

import hashlib
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

from hamopy.classes import Material
from hamopy import ham_library as ham


_SCALARS = {
    '1': ('rho', 'kg/m³'),
    '2': ('porosity', 'm³/m³'),
    '3': ('cp', 'J/kgK'),
    '4': ('lambda_10', 'W/mK'),
    '5': ('mu', '-'),
    '6': ('w80', 'kg/m³'),
    '9': ('wf', 'kg/m³'),
    '10': ('absorption_coefficient', 'kg/m²√s'),
    '11': ('drying_factor', '-'),
    '12': ('moisture_supplement', '%/M.-%'),
    '13': ('construction_moisture', 'kg/m³'),
    '14': ('temperature_supplement', 'W/mK²'),
}
_FUNCTIONS = {
    '1': ('FALSE', '-', 'kg/m³'),
    '2': ('TRUE', 'kg/m³', 'm²/s'),
    '3': ('TRUE', 'kg/m³', 'm²/s'),
    '4': ('TRUE', 'kg/m³', 'W/(m K)'),
    '5': ('FALSE', '-', '-'),
    '6': ('TRUE', '°C', 'W/(m K)'),
    '7': ('FALSE', '°C', 'J/kg'),
}


def _number(value):
    try:
        result = float(value.replace(',', '.'))
    except (AttributeError, ValueError) as error:
        raise ValueError('Missing or invalid WUFI numeric value: {!r}'.format(value)) from error
    if not np.isfinite(result):
        raise ValueError('WUFI numeric values must be finite.')
    return result


def _indexed(elements):
    result = {}
    for element in elements:
        key = element.get('id')
        if key is None or key in result:
            raise ValueError('Missing or duplicate WUFI scalar/function ID.')
        result[key] = element
    return result


def _mode(mode):
    if mode not in ('suction', 'redistribution'):
        raise ValueError("liquid_transport must be 'suction' or 'redistribution'.")
    return mode


class WufiMaterial(Material):
    """Material with tabulated sorption and generated WUFI liquid transport.

    Create instances with :func:`load_wufi_material`. ``liquid_transport``
    selects one curve for the whole simulation; Hamopy does not switch curves
    automatically with rain. ``source`` retains provenance and modelling notes.
    """

    def __init__(self, name, scalars, sorption, source, liquid_transport):
        super().__init__(name, rho=scalars['rho'], cp=scalars['cp'])
        self.source = source
        self.porosity = scalars['porosity']
        self.w80 = scalars['w80']
        self.free_saturation = scalars['wf']
        self.absorption_coefficient = scalars['absorption_coefficient']
        self.drying_factor = scalars['drying_factor']
        self.liquid_transport = _mode(liquid_transport)
        self.set_isotherm('tabulated', HR=sorption[:, 0], W=sorption[:, 1])
        self.set_perm_vapor('interp_mu', HR=[0., 1.], MU=[scalars['mu']] * 2)
        # WUFI lambda is referenced to 10 °C; Hamopy's intercept is at 0 °C.
        self.set_conduc(
            lambda_0=scalars['lambda_10'] - 10 * scalars['temperature_supplement'],
            lambda_m=1000 * scalars['lambda_10'] * scalars['moisture_supplement'] / self.rho,
            lambda_t=scalars['temperature_supplement'])

    def liquid_diffusivity(self, moisture, mode=None):
        """Generated D_w [m²/s] at water content [kg/m³].

        WUFI 4.1 help, pp. 31–34: D_ws = 3.8 (A/wf)^2 1000^(w/wf-1).
        Generated tables use w80 and wf as their positive anchors, with log-D
        interpolation. Redistribution divides only the wf value by the drying
        factor. Transport below w80 is zero; above wf the last value is held.
        https://wufi.de/download/WUFI41_help_e.pdf
        """
        mode = _mode(self.liquid_transport if mode is None else mode)
        moisture = np.asarray(moisture, dtype=float)
        if self.absorption_coefficient == 0:
            return np.zeros_like(moisture)
        saturation = 3.8 * (self.absorption_coefficient / self.free_saturation)**2
        reference = saturation * 1000**(self.w80 / self.free_saturation - 1)
        if mode == 'redistribution':
            saturation /= self.drying_factor
        fraction = np.clip((moisture - self.w80) / (self.free_saturation - self.w80), 0, 1)
        diffusivity = np.exp(np.log(reference) + fraction * np.log(saturation / reference))
        return np.where(moisture < self.w80, 0., diffusivity)

    def k_l(self, p_c):
        """Pressure-based liquid permeability [s], k_l = D_w * dw/dp_c."""
        if hasattr(self, 'kl_method'):
            return super().k_l(p_c)
        return self.liquid_diffusivity(self.w(p_c)) * self.c_w(p_c)


def load_wufi_material(path, *, liquid_transport='redistribution'):
    """Load a local SI WUFI XML file as a coupled heat/moisture material.

    Supported profile: tabulated moisture storage, constant mu and dry heat
    capacity, generated linear conductivity, and liquid coefficients generated
    from A, w80, wf and the drying factor (generation mode 0).
    """
    _mode(liquid_transport)
    path = Path(path)
    content = path.read_bytes()
    root = ET.fromstring(content)
    records = root.findall('Material')
    if root.tag != 'WUFIMaterial' or root.get('UnitSystem') != '0' or len(records) != 1:
        raise ValueError('Expected one WUFIMaterial record in SI units (UnitSystem=0).')
    record = records[0]
    if record.get('version') != '2' or not (record.findtext('Name') or '').strip():
        raise ValueError('Expected a named WUFI material, version 2.')
    scalars = _indexed(record.findall('Scalars/Scalar'))
    values = {}
    for key, (name, unit) in _SCALARS.items():
        scalar = scalars.get(key)
        if scalar is None or scalar.get('unit', '').strip() != unit:
            raise ValueError('Missing scalar {} ({}) or incompatible units; expected {}.'
                             .format(key, name, unit))
        values[name] = _number(scalar.get('value'))
    for name in ('rho', 'cp', 'lambda_10', 'mu', 'w80', 'wf', 'drying_factor'):
        if values[name] <= 0:
            raise ValueError('{} must be positive.'.format(name))
    for name in ('absorption_coefficient', 'moisture_supplement',
                 'temperature_supplement', 'construction_moisture'):
        if values[name] < 0:
            raise ValueError('{} must be nonnegative.'.format(name))
    if not 0 < values['porosity'] <= 1 or values['w80'] >= values['wf']:
        raise ValueError('Expected porosity in (0, 1] and 0 < w80 < wf.')

    functions = _indexed(record.findall('Functions/Function'))
    if set(functions) != set(_FUNCTIONS):
        raise ValueError('Expected WUFI functions 1 through 7.')
    tables = {}
    for key, (generated, xunit, yunit) in _FUNCTIONS.items():
        function = functions[key]
        items = function.find('Items')
        if (function.get('version') != '1' or
                (function.findtext('Generated') or '').strip() != generated or
                (function.findtext('GenerationMode') or '').strip() != '0' or
                items is None or items.get('version') != '2' or
                items.get('xunit', '').strip() != xunit or
                items.get('yunit', '').strip() != yunit):
            raise ValueError('Unsupported definition or units for WUFI function {}.'.format(key))
        points = items.findall('Item')
        if items.get('Count') != str(len(points)):
            raise ValueError('Incorrect item count for WUFI function {}.'.format(key))
        tables[key] = np.array([[_number(p.get('x')), _number(p.get('y'))]
                                for p in points], dtype=float).reshape(-1, 2)

    expected = {'2': (0, 0), '3': (0, 0), '4': (0, values['lambda_10']),
                '5': (0, values['mu']), '6': (10, values['lambda_10'])}
    for key, point in expected.items():
        if tables[key].shape != (1, 2) or not np.allclose(tables[key][0], point, rtol=1e-10, atol=0):
            raise ValueError('Unsupported table for WUFI function {}.'.format(key))
    if tables['7'].size and (tables['7'].shape != (1, 2) or np.any(tables['7'] != 0)):
        raise ValueError('Temperature-dependent enthalpy tables are not supported.')

    source = {
        'path': str(path.resolve()), 'sha256': hashlib.sha256(content).hexdigest(),
        'name': record.findtext('Name'), 'info': record.findtext('Info', ''),
        'date': record.findtext('DateDB'), 'db_id': record.get('dbID'),
        'version': record.get('version'), 'scalars': values.copy(),
    }
    material = WufiMaterial(record.findtext('Name').strip(), values, tables['1'],
                            source, liquid_transport)
    if (not np.isclose(material.w(ham.p_c(.8)), values['w80'], rtol=1e-6, atol=0) or
            not np.isclose(material.w_sat, values['wf'], rtol=1e-6, atol=0)):
        raise ValueError('Sorption table disagrees with the w80 or free-saturation scalar.')
    return material
