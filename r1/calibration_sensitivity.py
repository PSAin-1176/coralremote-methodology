"""Supplement S1 calibration sensitivities of the pooled StableSand pair.

Computes, on the calibration pixels of r1/revision_runner.py:
  - the largest polygon's pixel share and the pair without it (0.203 / 0.209), and the range when any single
    polygon is omitted;
  - Lyzenga's variance-covariance estimator on the same pixels (0.393 / 0.386);
  - the deep-water subtraction sensitivity table of S1: per-date deep-water reflectance (median, mean or
    10th percentile over the large-scale deep-water mask, or median over the deepest subset) subtracted
    before the log regression; pooled pair, SD of the eight single-date slopes, mean 2024 minus mean 2023
    slope, dates passing the reliability rule, pair without the largest polygon, observations lost, the
    paired polygon-bootstrap change in the absolute 2024 minus 2023 contrast, deep-water red as a share of
    sand red, and the share of valid reef pixels at or below the deep-water level.
Ported from the sensitivity analysis of the revision workspace (its tide-correlation and intertidal variants
are not part of the paper and are omitted); the computations are unchanged. Only coefficients, statistics
and counts are written.

Environment (required unless marked):
  CORALREMOTE_ROOT, CORALREMOTE_A3_ROOT, REVISION_UDM, CORALREMOTE_STABLESAND   as for revision_runner.py
  DEEP_MASK_BM       large-scale deep-water mask GeoTIFF (Phase 2 bm_deep_mask.tif: atlas depth > 3 m or the
                     manual deep-water polygon)
  DEEP_MASK_MDEEP    deepest-subset mask GeoTIFF (Phase 2 mdeep_mask.tif)
  CALIB_OUT          (optional) output folder, default r1/output/calibration_sensitivity
  CALIB_REFERENCE    (optional) folder of a previous run of this analysis; outputs are compared with the
                     files of the same name there (numbers to 1e-9)

Run from the repository root (a few minutes):  python r1/calibration_sensitivity.py
"""
import csv
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.argv = sys.argv[:1]
sys.path.insert(0, str(HERE))
import revision_runner as rr                      # noqa: E402
import numpy as np                                # noqa: E402
import rasterio                                   # noqa: E402

OUTD = Path(os.environ.get('CALIB_OUT') or HERE / 'output' / 'calibration_sensitivity')
BASE, DATES = rr.BASE, rr.NEEDED
D2023 = [d for d in DATES if d.startswith('2023')]; D2024 = [d for d in DATES if d.startswith('2024')]
FLOOR = 1e-4
BANDS = {'B': 1, 'G1': 2, 'R': 5}
PAIR = ('BR', 'GR')
GLINT_PCTL = 95          # NIR outlier cut inside the deep-water mask (glint and whitecaps)


def _need(name):
    v = os.environ.get(name)
    if not v:
        raise RuntimeError(f'{name} is not set (see the docstring of r1/calibration_sensitivity.py)')
    return Path(v)


def sand_geometry():
    """As revision_runner.sand_geometry, reading the StableSand file from CORALREMOTE_STABLESAND."""
    from rasterio.features import rasterize
    from rasterio.warp import transform_geom
    from affine import Affine
    ref = rr.io.load_glint_stack().ref; geo = rr.load(rr.stablesand_path())
    shapes = [(transform_geom('EPSG:4326', 'EPSG:32647', f['geometry']), i + 1) for i, f in reversed(list(enumerate(geo['features'])))]
    own = rasterize(shapes, out_shape=(ref.height, ref.width), transform=Affine(*ref.transform), dtype='int16')
    _, labels = rr.rasterize_labels_stage1(rr.io.load_geojson(str(rr.paths.STAGE1_LABELS['test'])), ref, rr.stage1.POLY_CLASS_MAP)
    r_, c_ = np.indices(own.shape); tiles = r_ // 32 * 10000 + c_ // 32
    own[np.isin(tiles, np.unique(tiles[labels > 0]))] = 0
    return own


def lyzenga(M):
    """Lyzenga's variance-covariance estimator on within-date-centred pooled moments."""
    z = M.sum(1); n = z[:, 0]; good = n > 1; z = z[good]; n = n[good]
    vx = np.sum(z[:, 4] - z[:, 1] ** 2 / n); cov = np.sum(z[:, 5:7] - z[:, 1, None] * z[:, 2:4] / n[:, None], axis=0)
    vy = np.sum(z[:, 7:9] - z[:, 2:4] ** 2 / n[:, None], axis=0)
    a = (vy - vx) / (2 * cov)
    return a + np.sqrt(a * a + 1)


def drop_one(M, npoly):
    out = []
    for p in range(npoly):
        if M[..., p, 0].sum() == 0:
            continue
        w = np.ones(npoly); w[p] = 0
        out.append((p + 1, rr.estimate(M, w)[0]))
    ks = np.array([k for _, k in out]); pid = [p for p, _ in out]
    return ks.min(0), ks.max(0), [pid[i] for i in ks.argmin(0)], [pid[i] for i in ks.argmax(0)]


def ci(a):
    return np.nanquantile(a, [.025, .975], axis=0)


def main():
    rr.log('loading geometry and reflectance')
    own = sand_geometry(); NPOLY = int(own.max())
    GL = {d: rr.pp.glint_reflectance(d) for d in DATES}
    RAW = {d: rr.pp.raw_stack().get(d) for d in DATES}
    UDM = {d: rr.pp.udm2_clear(d) for d in DATES}
    raw = {d: rr.raw_calibration(d, own) for d in DATES}
    BM = rasterio.open(_need('DEEP_MASK_BM')).read(1) > 0
    MDEEP = rasterio.open(_need('DEEP_MASK_MDEEP')).read(1) > 0
    reef = rr.pp.reef_mask()

    def deep_set(d, mask):
        A, R = GL[d], RAW[d]; inmask = mask & (own == 0)
        ok = inmask & UDM[d] & np.isfinite(A).all(2) & (np.isfinite(R) & (R > 0)).all(2)
        cut = float(np.percentile(R[:, :, 7][ok], GLINT_PCTL))
        return ok & ~(R[:, :, 7] > cut)

    LINF, deep_rows = {}, []
    for d in DATES:
        for mname, mask in (('BM', BM), ('MDeep', MDEEP)):
            ok = deep_set(d, mask); row = dict(date=d, mask=mname, n_used=int(ok.sum()))
            for b, i in BANDS.items():
                v = GL[d][:, :, i][ok].astype(float)
                for st, val in (('median', np.median(v)), ('mean', v.mean()), ('p10', np.percentile(v, 10))):
                    LINF.setdefault((mname, st), {}).setdefault(d, {})[b] = float(val); row[f'Linf_{st}_{b}'] = float(val)
            if mname == 'BM':
                valid = raw[d][0]
                for b, i in BANDS.items():
                    row[f'sand_median_{b}'] = float(np.median(GL[d][:, :, i][valid].astype(float)))
                row['ratio_Linf_median_R_to_sand_median_R'] = row['Linf_median_R'] / row['sand_median_R']
            deep_rows.append(row)

    def variant_pixels(d, linf=None):
        valid, x0, y0, ids = raw[d]; keep = np.ones(len(ids), bool); lost = {}
        if linf is None:
            x, y = x0, y0
        else:
            A = GL[d]; diff = {b: A[:, :, i][valid].astype(float) - linf[d][b] for b, i in BANDS.items()}
            both = (diff['B'] > FLOOR) & (diff['G1'] > FLOOR) & (diff['R'] > FLOOR)
            lost = dict(lost_any=int((keep & ~both).sum()), n_before_subtraction=int(keep.sum())); keep &= both
            x = np.log(np.where(both, diff['R'], 1.0)); y = np.stack([np.log(np.where(both, diff[b], 1.0)) for b in ('B', 'G1')], 1)
        return x[keep], y[keep], ids[keep], lost

    VARIANTS = {'A_baseline': None, 'B_median_BM': ('BM', 'median'), 'B_mean_BM': ('BM', 'mean'),
                'B_p10_BM': ('BM', 'p10'), 'B_median_MDeep': ('MDeep', 'median')}
    results, per_date, BS = [], [], {}
    i23 = [DATES.index(d) for d in D2023]; i24 = [DATES.index(d) for d in D2024]
    for name, key in VARIANTS.items():
        rr.log(f'variant {name}')
        linf = LINF[key] if key else None
        pix = {d: variant_pixels(d, linf) for d in DATES}
        M = {d: rr.moments(*pix[d][:3], NPOLY) for d in DATES}; Mb = np.stack([M[d] for d in BASE])
        pool, bs_pool = rr.fit_summary(Mb); dmin, dmax, dargmin, dargmax = drop_one(Mb, NPOLY)
        largest = int(np.argmax(Mb[:, :, 0].sum(0)))
        w = np.ones(NPOLY); w[largest] = 0; nolargest = rr.estimate(Mb, w)[0]
        share = float(Mb[:, largest, 0].sum() / Mb[:, :, 0].sum())
        fits, bsd = {}, {}
        for d in DATES:
            fits[d], bsd[d] = rr.fit_summary(M[d][None])
        S = np.array([fits[d]['pair'] for d in DATES]); SB = np.stack([bsd[d] for d in DATES], 1)
        sd = S.std(0, ddof=1); con = S[i24].mean(0) - S[i23].mean(0)
        con_b = np.nanmean(SB[:, i24], 1) - np.nanmean(SB[:, i23], 1); BS[name] = con_b
        lost_base = {k: sum(pix[d][3].get(k, 0) for d in BASE) for k in ('lost_any', 'n_before_subtraction')}
        lyz = lyzenga(Mb)
        for j, band in enumerate(PAIR):
            results.append(dict(variant=name, band=band, pooled_slope=pool['pair'][j], pooled_lower95=pool['CI95'][0][j],
                                 pooled_upper95=pool['CI95'][1][j], pooled_R2=pool['R2'][j], pooled_n=pool['n'],
                                 lyzenga_vc_estimator=lyz[j], drop1poly_min=dmin[j], drop1poly_min_poly=dargmin[j],
                                 drop1poly_max=dmax[j], drop1poly_max_poly=dargmax[j], largest_poly=largest + 1,
                                 pooled_without_largest_poly=nolargest[j], largest_poly_pixel_share=share,
                                 perdate_sd=sd[j], contrast_2024_minus_2023=con[j], contrast_lower95=ci(con_b)[0][j],
                                 contrast_upper95=ci(con_b)[1][j], n_dates_eligible=sum(fits[d]['eligible'] for d in DATES),
                                 pixels_lost_nonpositive_baseline3=lost_base['lost_any'],
                                 pixels_before_subtraction_baseline3=lost_base['n_before_subtraction']))
            for d in DATES:
                f = fits[d]
                per_date.append(dict(variant=name, date=d, band=band, slope=f['pair'][j], lower95=f['CI95'][0][j],
                                     upper95=f['CI95'][1][j], R2=f['R2'][j], n=f['n'], polygons=f['occupied_polygons'],
                                     eligible=f['eligible'], failure_reasons=';'.join(f['failure_reasons'])))
    paired = []
    for name in VARIANTS:
        if name == 'A_baseline':
            continue
        for j, band in enumerate(PAIR):
            dlt = np.abs(BS[name])[:, j] - np.abs(BS['A_baseline'])[:, j]; q = np.nanquantile(dlt, [.025, .5, .975])
            paired.append(dict(variant=name, band=band, metric='abs_contrast_2024_minus_2023 (variant minus A)',
                               median=q[1], lower95=q[0], upper95=q[2]))
    reef_rows = []
    for d in DATES:
        F = rr.pp.phase5_features(d); ok = reef & np.isfinite(F).all(2); row = dict(date=d, reef_valid_px=int(ok.sum()))
        for b, i in BANDS.items():
            row[f'frac_le_Linf_median_BM_{b}'] = float(np.mean(GL[d][:, :, i][ok] - LINF[('BM', 'median')][d][b] <= FLOOR))
        reef_rows.append(row)
    rr.guard_out(OUTD); OUTD.mkdir(parents=True, exist_ok=True)
    rr.csvsave(OUTD / 'results.csv', results); rr.csvsave(OUTD / 'per_date.csv', per_date)
    rr.csvsave(OUTD / 'deepwater.csv', deep_rows); rr.csvsave(OUTD / 'paired_bootstrap.csv', paired)
    rr.csvsave(OUTD / 'reef_below_deepwater.csv', reef_rows)
    rr.save(OUTD / 'design.json', dict(dates=DATES, calibration_dates=BASE, raw_NIR_interval=[.005, .12], floor_reflectance=FLOOR,
                                       estimator='OLS, within-date centred pooled fit (revision_runner.estimate)',
                                       bootstrap='revision_runner.bootstrap: 1000 polygon resamples, seed 20261118; same draws for every date and variant',
                                       deep_pixel_filters=['inside the mask and not a StableSand pixel', 'UDM2 clear',
                                                           'glint product finite in all bands; raw reflectance finite and > 0 in all bands',
                                                           f'raw NIR at or below the date-specific {GLINT_PCTL}th percentile inside the mask'],
                                       inputs=[dict(name=n, sha256=rr.sha(p)) for n, p in (('DEEP_MASK_BM', _need('DEEP_MASK_BM')),
                                                                                      ('DEEP_MASK_MDEEP', _need('DEEP_MASK_MDEEP')),
                                                                                      ('StableSand', rr.stablesand_path()))]))
    a = {x['band']: x for x in results if x['variant'] == 'A_baseline'}
    rr.log(f"pooled {a['BR']['pooled_slope']:.4f}/{a['GR']['pooled_slope']:.4f}; Lyzenga {a['BR']['lyzenga_vc_estimator']:.4f}/"
           f"{a['GR']['lyzenga_vc_estimator']:.4f}; largest polygon {a['BR']['largest_poly']} share {a['BR']['largest_poly_pixel_share']:.3f}")


def compare(ref):
    ref = Path(ref); ok = True
    for f in ('results.csv', 'per_date.csv', 'deepwater.csv', 'paired_bootstrap.csv', 'reef_below_deepwater.csv'):
        new = list(csv.DictReader((OUTD / f).open(encoding='utf-8'))); old = list(csv.DictReader((ref / f).open(encoding='utf-8')))
        rename = {'pooled_without_largest_poly': 'pooled_without_poly36', 'largest_poly_pixel_share': 'poly36_pixel_share'}

        def match(n, o):
            for k, v in n.items():
                ko = rename.get(k, k)
                if ko not in o:
                    continue
                if v != o[ko]:
                    try:
                        if abs(float(v) - float(o[ko])) > 1e-9 * max(1, abs(float(o[ko]))):
                            return False
                    except ValueError:
                        return False
            return True
        keyf = [k for k in ('variant', 'date', 'band', 'mask', 'metric') if new and k in new[0]]
        index = {tuple(o.get(k) for k in keyf): o for o in old}
        diffs = [n for n in new if tuple(n.get(k) for k in keyf) not in index or not match(n, index[tuple(n.get(k) for k in keyf)])]
        print(f"{'match' if not diffs else 'DIFFERS':8s} {f} ({len(new)} rows)")
        ok &= not diffs
    return ok


if __name__ == '__main__':
    main()
    if os.environ.get('CALIB_REFERENCE'):
        sys.exit(0 if compare(os.environ['CALIB_REFERENCE']) else 1)
