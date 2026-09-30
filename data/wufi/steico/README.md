# STEICO WUFI material files

Retrieved on 2026-09-29 from STEICO's official
[technical downloads page](https://www.steico.com/de/downloads/dokumente/technik-verarbeitung),
under **WUFI Datensätze als XML (ab WUFI Pro 5.0)**.

The manufacturer's collection contains 19 material XML files and 19 accompanying
PDF descriptions. This repository records their provenance and checksums. Local
downloads are excluded from version control because no redistribution licence
was identified in the archive.

- [Material index, source URL, dates and SHA-256 checksums](manifest.json)
- [Manufacturer download](https://www.steico.com/fileadmin/user_upload/importer/downloads/wufi_datenstze_als_xml_ab_wufi_pro_0/STEICO_WUFI-data-set_complete_EUR_de_i.zip)

## Download locally

Run from this directory:

```sh
curl --fail --location \
  'https://www.steico.com/fileadmin/user_upload/importer/downloads/wufi_datenstze_als_xml_ab_wufi_pro_0/STEICO_WUFI-data-set_complete_EUR_de_i.zip' \
  --output steico-wufi-datasets.zip
```

Then verify and extract the original bytes using the manifest:

```python
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json

manifest = json.loads(Path('manifest.json').read_text())
archive_path = Path(manifest['archive'])
assert hashlib.sha256(archive_path.read_bytes()).hexdigest() == manifest['archive_sha256']
with ZipFile(archive_path) as archive:
    for material in manifest['materials']:
        for kind in ('xml', 'pdf'):
            member = str(Path(material['archive_member']).with_suffix('.' + kind))
            content = archive.read(member)
            assert hashlib.sha256(content).hexdigest() == material[kind + '_sha256']
            destination = Path(material[kind])
            destination.parent.mkdir(exist_ok=True)
            destination.write_bytes(content)
```

The XML files appear in `xml/` and matching descriptions in `pdf/`. If the archive
checksum differs, STEICO may have revised the download; inspect the new version
before updating the manifest.

## Materials

| XML material name | Python library identifier |
| --- | --- |
| STEICOinstall | `steico_install` |
| STEICOisorel | `steico_isorel` |
| STEICOphaltex | `steico_phaltex` |
| STEICOprotect H | `steico_protect_h` |
| STEICOprotect H dry | `steico_protect_h_dry` |
| STEICOprotect L dry | `steico_protect_l_dry` |
| STEICOprotect M | `steico_protect_m` |
| STEICOroof dry | `steico_roof_dry` |
| STEICOsafe (120-160mm) | `steico_safe_120_160mm` |
| STEICOsafe (40mm) | `steico_safe_40mm` |
| STEICOsafe (60-100mm) | `steico_safe_60_100mm` |
| STEICOspecial dry | `steico_special_dry` |
| STEICOtherm | `steico_therm` |
| STEICOtherm dry | `steico_therm_dry` |
| STEICOtherm SD | `steico_therm_sd` |
| STEICOtop | `steico_top` |
| STEICOuniversal | `steico_universal` |
| STEICOuniversal dry (35mm) | `steico_universal_dry_35mm` |
| STEICOuniversal dry (60-100mm) | `steico_universal_dry_40_100mm` |

Names above are taken from the XML records. The last record's filename and
website listing say 40-100mm, while its internal name says 60-100mm; both are
preserved in the source files and index.

**STEICOinternal is not present.** This is STEICO's published bundle, not the
complete material database bundled with WUFI.

## Use in Hamopy

After downloading, all 19 records can be used as coupled simulation materials:

```python
from hamopy.materials.steico import steico_therm, load_materials

materials = load_materials()  # keyed by the identifiers above
```

For XMLs outside the source checkout, pass their folder to `load_materials` or
set `HAMOPY_STEICO_XML_DIR` before using named imports. See the
[material guide](../../../docs/user-guide/materials.md#steico-materials-from-wufi-xml)
for liquid-transport selection, conversion equations and supported XML models.

## Import notes

- All 19 XML files parse successfully and have a matching PDF description.
- Numeric values use decimal commas, and source units are recorded in the XML.
- Each material contains a 12-point moisture-storage table. The importer keeps
  all points and uses monotone interpolation with an analytic moisture capacity.
- Liquid absorption and redistribution functions are marked `Generated=TRUE`
  and contain a single zero placeholder. Their actual curves must be generated
  from the associated parameters. The importer reconstructs these curves and
  converts diffusivity to pressure-based permeability. Redistribution is the
  default; suction can be selected explicitly.
- Moisture- and temperature-dependent conductivity functions also use generated
  definitions. The importer converts their reference temperature and moisture
  supplement. Vapour resistance has a single constant point.
- The XML `Info` fields retain the source's provenance and modelling limitations.
  These are retained in each material's `source` metadata. Download availability
  alone does not establish redistribution rights; no redistribution licence was
  identified in the archive.

WUFI's [XML import instructions](https://wufi.de/en/service/downloads/) are listed
under **Material files**.

## Verification

Run `python -m unittest discover -s tests -v` from the repository root. Synthetic
fixtures test the conversion without downloading manufacturer data. With all 19
local XMLs present, additional tests check every sorption table and run one-hour
coupled wetting and drying simulations for each material. These check numerical
integration and transport direction, not agreement with measured performance or
WUFI simulation results.
