"""Supplement S7 scene list from the PlanetScope item metadata.

Reads <scene id>.json (Planet STAC item records) for the 14 acquisitions from PLANET_SCENES and writes
scene identifier, acquisition time (UTC), item type, instrument, sun elevation and azimuth and view angle
to SCENES_OUT (default r1/output/scenes.csv). The committed r1/scenes.csv holds the identifiers and times
only; imagery and metadata files are licensed from Planet Labs PBC and are not redistributed.

Run from the repository root:  python r1/scenes.py
"""
import csv
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCENES = [x['scene_id'] for x in csv.DictReader((HERE / 'scenes.csv').open(encoding='utf-8'))]


def main():
    src = os.environ.get('PLANET_SCENES')
    if not src:
        raise RuntimeError('PLANET_SCENES is not set: folder holding <scene id>.json item records')
    rows = []
    for sid in SCENES:
        j = json.loads((Path(src) / f'{sid}.json').read_text(encoding='utf-8')); p = j['properties']
        assert j['id'] == sid
        rows.append(dict(date=sid[:8], scene_id=sid, acquired_utc=p['datetime'], item_type=p['pl:item_type'], instrument=';'.join(p['instruments']),
                         sun_elevation_deg=p.get('view:sun_elevation'), sun_azimuth_deg=p.get('view:sun_azimuth'), off_nadir_deg=p.get('view:off_nadir')))
    out = Path(os.environ.get('SCENES_OUT') or HERE / 'output' / 'scenes.csv')
    out = out.resolve()
    if out.is_relative_to(Path(src).resolve()):
        raise PermissionError(f'{out} is inside the input folder {src}')
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f'{len(rows)} scenes written to {out}')


if __name__ == '__main__':
    main()
