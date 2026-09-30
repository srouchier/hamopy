"""STEICO material catalog built from the manufacturer's local WUFI XML files.

Use ``load_materials(xml_directory)`` for the full catalog or import named
materials, e.g. ``from hamopy.materials.steico import steico_therm``. Each load
returns independent objects. No files are downloaded during import.
"""

import os
from pathlib import Path

from .wufi import load_wufi_material


MATERIAL_FILES = {
    'steico_' + suffix.replace('-', '_'): 'STEICO_WUFI-data-set_' + suffix + '_EUR_de_i.xml'
    for suffix in (
        'install', 'isorel', 'phaltex', 'protect-H', 'protect-H-dry',
        'protect-L-dry', 'protect-M', 'roof-dry', 'safe-120-160mm',
        'safe-40mm', 'safe-60-100mm', 'special-dry', 'therm', 'therm-dry',
        'therm-SD', 'top', 'universal', 'universal-dry-35mm',
        'universal-dry-40-100mm',
    )
}
# Public Python names are lowercase; source filenames retain their exact case.
MATERIAL_FILES = {name.lower(): filename for name, filename in MATERIAL_FILES.items()}


def load_material(name, xml_directory=None, *, liquid_transport='redistribution'):
    """Build one named material from a local XML, optionally selecting suction.

    Files are read from ``xml_directory``, then ``HAMOPY_STEICO_XML_DIR``, or
    ``data/wufi/steico/xml`` in a source checkout. See that folder's README for
    manufacturer download and checksum verification instructions.
    """
    filename = MATERIAL_FILES[name]
    if xml_directory is None:
        xml_directory = os.environ.get(
            'HAMOPY_STEICO_XML_DIR',
            Path(__file__).resolve().parents[2] / 'data' / 'wufi' / 'steico' / 'xml')
    path = Path(xml_directory).expanduser() / filename
    if not path.is_file():
        raise FileNotFoundError(
            'Missing STEICO XML: {}. Download and extract the manufacturer files '
            'as described in data/wufi/steico/README.md, then pass xml_directory '
            'or set HAMOPY_STEICO_XML_DIR.'.format(path))
    return load_wufi_material(path, liquid_transport=liquid_transport)


def load_materials(xml_directory=None, *, liquid_transport='redistribution'):
    """Build all 19 materials, keyed by Python name; fail if any file is absent."""
    return {name: load_material(name, xml_directory, liquid_transport=liquid_transport)
            for name in MATERIAL_FILES}


def __getattr__(name):
    if name in MATERIAL_FILES:
        return load_material(name)
    raise AttributeError('module {!r} has no attribute {!r}'.format(__name__, name))


def __dir__():
    return sorted(set(globals()) | set(MATERIAL_FILES))


__all__ = ['MATERIAL_FILES', 'load_material', 'load_materials'] + list(MATERIAL_FILES)
