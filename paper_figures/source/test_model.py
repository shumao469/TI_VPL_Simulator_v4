"""Numerical checks for figure calculations, not biological validation."""
import unittest
from dataclasses import replace
import numpy as np
from model import Atlas, DEFAULT, TARGET, preview_fields, maximum_orientation_envelope


class FigureModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.atlas=Atlas()
        cls.points=cls.atlas.points[::791]

    def test_atlas_masks_and_coordinate_roundtrip(self):
        a=self.atlas
        self.assertEqual(int(a.vpl.sum()),441)
        self.assertTrue(np.all(a.brain[a.vpl]))
        self.assertTrue(np.all(a.sample_mask(a.vpl_points,'vpl')))
        self.assertAlmostEqual(a.distance,.6356886030125125,places=12)

    def test_field_linear_current_scaling(self):
        a=preview_fields(self.points,DEFAULT)
        p=replace(DEFAULT,channel_1_current_mA=1.,channel_2_current_mA=1.)
        b=preview_fields(self.points,p)
        for name in ['E1','E2','mTI']:
            np.testing.assert_allclose(b[name],2*a[name],rtol=1e-12,atol=1e-12)

    def test_nonzero_frequency_difference_does_not_change_space(self):
        ref=preview_fields(self.points,DEFAULT)['mTI']
        for f2 in [2001.,2040.,2100.]:
            got=preview_fields(self.points,replace(DEFAULT,f2_Hz=f2))['mTI']
            np.testing.assert_array_equal(got,ref)

    def test_equal_frequencies_and_single_channel_have_no_ti(self):
        for p in [replace(DEFAULT,f2_Hz=2000.),replace(DEFAULT,channel_2_current_mA=0.)]:
            np.testing.assert_array_equal(preview_fields(self.points,p)['mTI'],0.)

    def test_envelope_against_direction_search(self):
        rng=np.random.default_rng(73)
        x=rng.normal(size=(40,3));y=rng.normal(size=(40,3))
        n=rng.normal(size=(20000,3));n/=np.linalg.norm(n,axis=1)[:,None]
        sampled=(2*np.minimum(np.abs(x@n.T),np.abs(y@n.T))).max(1)
        exact=maximum_orientation_envelope(x,y)
        self.assertTrue(np.all(sampled<=exact+1e-12))
        self.assertLess(float(np.max((exact-sampled)/exact)),.025)

    def test_expected_orthogonal_and_collinear_limits(self):
        np.testing.assert_allclose(maximum_orientation_envelope(np.array([[1.,0,0]]),np.array([[0,1.,0]])),np.sqrt(2))
        np.testing.assert_allclose(maximum_orientation_envelope(np.array([[2.,0,0]]),np.array([[-1.,0,0]])),2.)
        np.testing.assert_allclose(maximum_orientation_envelope(np.zeros((1,3)),np.array([[1.,2,3]])),0.)

    def test_target_is_not_snapped_into_vpl(self):
        a=self.atlas
        self.assertFalse(a.sample_mask(TARGET[None],'vpl')[0])
        self.assertAlmostEqual(float(preview_fields(TARGET[None],DEFAULT)['mTI'][0]),.3316171396768962,places=12)
        peak=a.points[a.fields['mTI'].argmax()]
        self.assertGreater(np.linalg.norm(peak-TARGET),2.)

    def test_normalization_is_one_fixed_reference(self):
        a=self.atlas
        self.assertAlmostEqual(a.reference,3.0793786491321113,places=11)
        self.assertGreater(a.reference,float(a.fields['mTI'][a.vpl_selection].max()))
        for kind in ['axial','sagittal','coronal']:
            d=a.plane(kind)
            self.assertTrue(d['vpl'].any())
            self.assertTrue(np.all(np.isfinite(d['mTI'][d['brain']])))


if __name__=='__main__':unittest.main()
