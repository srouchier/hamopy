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

- STEICOinstall
- STEICOisorel
- STEICOphaltex
- STEICOprotect H
- STEICOprotect H dry
- STEICOprotect L dry
- STEICOprotect M
- STEICOroof dry
- STEICOsafe (120-160mm)
- STEICOsafe (40mm)
- STEICOsafe (60-100mm)
- STEICOspecial dry
- STEICOtherm
- STEICOtherm dry
- STEICOtherm SD
- STEICOtop
- STEICOuniversal
- STEICOuniversal dry (35mm)
- STEICOuniversal dry (60-100mm)

Names above are taken from the XML records. The last record's filename and
website listing say 40-100mm, while its internal name says 60-100mm; both are
preserved in the source files and index.

**STEICOinternal is not present.** This is STEICO's published bundle, not the
complete material database bundled with WUFI.

## Import notes

These are source files for a future importer; downloading them does not register
them as Hamopy simulation materials.

- All 19 XML files parse successfully and have a matching PDF description.
- Numeric values use decimal commas, and source units are recorded in the XML.
- Each material contains a 12-point moisture-storage table.
- Liquid absorption and redistribution functions are marked `Generated=TRUE`
  and contain a single zero placeholder. Their actual curves must be generated
  from the associated parameters, not read as zero liquid transport.
- Moisture- and temperature-dependent conductivity functions also use generated
  definitions. Vapour resistance has a single constant point.
- The XML `Info` fields retain the source's provenance and modelling limitations.
  Preserve these when converting the data. Download availability alone does not
  establish redistribution rights; no redistribution licence was identified in
  the archive.

WUFI's [XML import instructions](https://wufi.de/en/service/downloads/) are listed
under **Material files**.
