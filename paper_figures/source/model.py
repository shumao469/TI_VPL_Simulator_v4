"""Reproducible v4 analytic-model data on the bundled Allen CCF masks.

The atlas is used for display/sampling only, not as a conductor boundary.
This module does not perform a finite-element solution or animal registration.
"""
from __future__ import annotations
import base64
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from ti_vpl.analytic import PreviewParameters, preview_fields
from ti_vpl.coordinates import crossed_montage
from ti_vpl.envelope import maximum_orientation_envelope

TARGET = np.array([-1.94, -1.75, -3.10])
DEFAULT = PreviewParameters()


class Atlas:
    def __init__(self, path=None):
        self.path = Path(path or ROOT / 'assets/allen_ccf_100um_web.json')
        self.payload = json.loads(self.path.read_text())
        self.shape = tuple(self.payload['shape_AP_DV_ML'])
        for key, name in [('brain_mask_bitpack_little_base64', 'brain'),
                          ('left_vpl_mask_bitpack_little_base64', 'vpl')]:
            values = np.frombuffer(base64.b64decode(self.payload[key]), dtype=np.uint8)
            setattr(self, name, np.unpackbits(values, bitorder='little')[:np.prod(self.shape)].reshape(self.shape).astype(bool))
        self.indices = np.argwhere(self.brain)
        self.points = self.to_world(self.indices)
        self.vpl_indices = np.argwhere(self.vpl)
        self.vpl_points = self.to_world(self.vpl_indices)
        self.vpl_centroid = self.vpl_points.mean(0)
        self.nearest = self.vpl_points[np.linalg.norm(self.vpl_points - TARGET, axis=1).argmin()]
        self.distance = float(np.linalg.norm(self.nearest - TARGET))
        self.slices = np.round(self.vpl_centroid, 1)
        self.fields = preview_fields(self.points, DEFAULT)
        self.reference = float(self.fields['mTI'].max())
        self.vpl_selection = self.vpl[self.brain]

    @staticmethod
    def to_world(i):
        return np.column_stack([5.4 - i[:, 0] * .1, i[:, 2] * .1 - 5.7, -i[:, 1] * .1])

    def sample_mask(self, points, name='brain'):
        shp = points.shape[:-1]
        p = np.asarray(points).reshape(-1, 3)
        ijk = np.rint(np.column_stack([(5.4-p[:, 0])/.1, -p[:, 2]/.1, (p[:, 1]+5.7)/.1])).astype(int)
        inside = np.all((ijk >= 0) & (ijk < self.shape), axis=1)
        result = np.zeros(len(p), bool)
        valid = ijk[inside]
        result[inside] = getattr(self, name)[valid[:, 0], valid[:, 1], valid[:, 2]]
        return result.reshape(shp)

    def plane(self, kind, coordinate=None, params=DEFAULT, nx=280):
        if kind == 'sagittal':
            coordinate = float(self.slices[1] if coordinate is None else coordinate)
            x, y = np.linspace(-7.8, 5.5, nx), np.linspace(-8, 0, int(nx*.62))
            xx, yy = np.meshgrid(x, y)
            xyz = np.stack([xx, np.full_like(xx, coordinate), yy], -1)
            labels, k = ('AP (mm)', 'DV (mm)'), 1
        elif kind == 'coronal':
            coordinate = float(self.slices[0] if coordinate is None else coordinate)
            x, y = np.linspace(-5.6, 5.6, nx), np.linspace(-8, 0, int(nx*.73))
            xx, yy = np.meshgrid(x, y)
            xyz = np.stack([np.full_like(xx, coordinate), xx, yy], -1)
            labels, k = ('ML (mm)', 'DV (mm)'), 0
        elif kind == 'axial':
            coordinate = float(self.slices[2] if coordinate is None else coordinate)
            x, y = np.linspace(-5.6, 5.6, nx), np.linspace(-7.8, 5.5, int(nx*1.18))
            xx, yy = np.meshgrid(x, y)
            xyz = np.stack([yy, xx, np.full_like(xx, coordinate)], -1)
            labels, k = ('ML (mm)', 'AP (mm)'), 2
        else:
            raise ValueError(kind)
        f = preview_fields(xyz, params)
        return dict(kind=kind, coordinate=coordinate, axis=k, x=x, y=y, xyz=xyz,
                    brain=self.sample_mask(xyz), vpl=self.sample_mask(xyz, 'vpl'),
                    mTI=f['mTI'], E1=np.linalg.norm(f['E1'], axis=-1),
                    E2=np.linalg.norm(f['E2'], axis=-1), labels=labels)

    def summary(self):
        field = self.fields['mTI']
        peak = self.points[field.argmax()]
        vpl_values = field[self.vpl_selection]
        high = field >= np.quantile(field, .95)
        centroid = np.average(self.points[high], weights=field[high], axis=0)
        target = float(preview_fields(TARGET[None], DEFAULT)['mTI'][0])
        return dict(scientific_status='ANALYTIC_PREVIEW_NOT_REGISTERED_FEM',
                    anatomy='Allen CCFv3 2017, 100 um; VPL id 718; left hemisphere',
                    atlas_asset_sha256=hashlib.sha256(self.path.read_bytes()).hexdigest(),
                    parameters=asdict(DEFAULT), requested_target_mm=TARGET.tolist(),
                    vpl_centroid_mm=self.vpl_centroid.tolist(),
                    nearest_vpl_voxel_center_mm=self.nearest.tolist(),
                    target_to_nearest_vpl_voxel_center_mm=self.distance,
                    brain_voxel_count=len(field), vpl_voxel_count=int(self.vpl.sum()),
                    voxel_volume_mm3=.001,
                    baseline_brain_max_preview_V_per_m=self.reference,
                    baseline_brain_peak_coordinate_mm=peak.tolist(),
                    baseline_target_preview_V_per_m=target,
                    baseline_vpl_mean_preview_V_per_m=float(vpl_values.mean()),
                    baseline_vpl_median_preview_V_per_m=float(np.median(vpl_values)),
                    baseline_brain_p95_preview_V_per_m=float(np.quantile(field, .95)),
                    vpl_median_over_brain_p95=float(np.median(vpl_values)/np.quantile(field, .95)),
                    top_5_percent_field_weighted_centroid_mm=centroid.tolist(),
                    normalization='Every spatial map is divided by the same baseline whole-brain maximum; no per-panel normalization.',
                    graphical_color_scale='Logarithmic, 0.01-1 for mTI; channel comparison uses 0.01-3.',
                    surface_rendering='8 iterations of 0.15 Laplacian smoothing for display only; field and ROI arrays unchanged.',
                    limitations=['Generic atlas mapping is not animal registration.',
                                 'Point sources are at the ideal DV +0.50 mm plane, not a reconstructed skull surface.',
                                 'The scalar factor 0.028 is inherited from v4 and is not a calibrated skull-transfer model.',
                                 'Brain anatomy masks the output but does not constrain the PDE.',
                                 'mTI is an orientation-maximized modulation depth, not activation probability.',
                                 '0.1 mm atlas sampling is not an FEM mesh-convergence study.'])


def generate_data(out, atlas):
    out.mkdir(parents=True, exist_ok=True)
    summary = atlas.summary()
    (out / 'model_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    volume = np.full(atlas.shape, np.nan, np.float32)
    volume[atlas.brain] = atlas.fields['mTI']
    np.savez_compressed(out / 'baseline_volume.npz', mTI_preview_V_per_m=volume,
                        brain=atlas.brain, left_vpl=atlas.vpl,
                        spacing_AP_DV_ML_mm=np.array([.1,.1,.1]),
                        generic_ccf_origin_AP_DV_ML_mm=np.array([5.4,0,5.7]),
                        reference_preview_V_per_m=atlas.reference)
    import csv
    rows=[]
    for i1 in np.linspace(.05, .95, 19):
        p=replace(DEFAULT,channel_1_current_mA=float(i1),channel_2_current_mA=float(1-i1))
        field=preview_fields(atlas.points,p)['mTI']
        rows.append(dict(I1_mA=i1,I2_mA=1-i1,
                         vpl_mean_relative=float(field[atlas.vpl_selection].mean()/atlas.reference),
                         target_relative=float(preview_fields(TARGET[None],p)['mTI'][0]/atlas.reference),
                         brain_p95_relative=float(np.quantile(field,.95)/atlas.reference)))
    with (out/'current_sweep.csv').open('w',newline='') as file:
        w=csv.DictWriter(file,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    slices={}
    for kind in ['sagittal','coronal','axial']:
        for label,coord in [('vpl',None),('target',TARGET[dict(sagittal=1,coronal=0,axial=2)[kind]])]:
            d=atlas.plane(kind,coord)
            for key in ['x','y','brain','vpl','mTI','E1','E2']:
                slices[f'{kind}_{label}_{key}']=d[key]
    np.savez_compressed(out/'orthogonal_slices.npz',**slices)
    return summary
