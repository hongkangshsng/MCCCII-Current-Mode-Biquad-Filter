import unittest
import numpy as np
from scipy.signal import freqs
from src.filters import BASE,coefficients,metrics,response,nodal_response,sweep

class Filters(unittest.TestCase):
    def test_baseline(self):
        for k in BASE:
            w,f,q=metrics(k,BASE[k]); self.assertAlmostEqual(w,1000 if k!='BPF' else 1000/np.sqrt(3))
            self.assertAlmostEqual(q,.5 if k!='BPF' else np.sqrt(3)/5)
    def test_independent_nodal_solver(self):
        rng=np.random.default_rng(42)
        for k,p in BASE.items():
            for _ in range(8):
                varied={key:v*10**rng.uniform(-1,1) for key,v in p.items()}; f=np.logspace(-2,6,71)
                np.testing.assert_allclose(response(k,varied,f),nodal_response(k,varied,f),rtol=2e-9,atol=1e-11)
    def test_scipy(self):
        for k,p in BASE.items():
            f=np.logspace(-2,6,101); n,d=coefficients(k,p)
            np.testing.assert_allclose(response(k,p,f),freqs(n,d,2*np.pi*f)[1],rtol=1e-12)
    def test_limits_and_phase(self):
        self.assertAlmostEqual(abs(response('LPF',BASE['LPF'],[1e-8])[0]),1)
        self.assertAlmostEqual(abs(response('HPF',BASE['HPF'],[1e10])[0]),1)
        for k in BASE:
            f0=metrics(k,BASE[k])[1]; h=response(k,BASE[k],[f0])[0]
            self.assertAlmostEqual(np.angle(h,deg=True),{'LPF':-90,'HPF':90,'BPF':0}[k])
        self.assertAlmostEqual(abs(response('BPF',BASE['BPF'],[metrics('BPF',BASE['BPF'])[1]])[0]),.2)
    def test_stability_and_tuning(self):
        for k,p in BASE.items():
            for v in [.25,.5,1,2,4]:
                for mode in ['q','frequency']:
                    varied=sweep(k,mode,v); self.assertTrue(np.all(np.real(np.roots(coefficients(k,varied)[1]))<0))
                    w,_,q=metrics(k,varied); wbase,_,qbase=metrics(k,p)
                    self.assertAlmostEqual(w/wbase,v if mode=='frequency' else 1)
                    if mode=='frequency': self.assertAlmostEqual(q,qbase)
    def test_bpf_bandwidth(self):
        p=BASE['BPF']; w,_,q=metrics('BPF',p); bw=w/q
        lo=(np.sqrt(bw*bw+4*w*w)-bw)/2; hi=lo+bw
        h=abs(response('BPF',p,np.array([lo,hi])/(2*np.pi)))
        np.testing.assert_allclose(h,.2/np.sqrt(2),rtol=1e-12)
        self.assertAlmostEqual(w/(hi-lo),q)
    def test_db_and_invalid(self):
        self.assertAlmostEqual(20*np.log10(.5),-6.020599913279624)
        self.assertAlmostEqual(20*np.log10(1/np.sqrt(2)),-3.010299956639812)
        for k,p in BASE.items():
            for bad in [0,-1,np.nan,np.inf]:
                varied=p.copy(); varied[next(iter(p))]=bad
                with self.assertRaises(ValueError): coefficients(k,varied)
        with self.assertRaises(ValueError): response('LPF',BASE['LPF'],[0])
if __name__=='__main__': unittest.main()
