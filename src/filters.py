"""Ideal Rx=0 MCCCII transfer functions; SI units throughout."""
import numpy as np

BASE = {
    'LPF': dict(G2=1e-3, G4=1e-3, C1=1e-6, C3=1e-6),
    'HPF': dict(G1=1e-3, G3=1e-3, C2=1e-6, C4=1e-6),
    'BPF': dict(G1=1e-3, G3=1e-3, G4=1e-3, C1=1e-6, C2=1e-6, C3=1e-6),
}

def coefficients(kind, p):
    if set(p) != set(BASE[kind]) or any(not np.isfinite(v) or v <= 0 for v in p.values()):
        raise ValueError('Expected exactly the positive finite SI parameters for '+kind)
    if kind == 'LPF':
        a=p['C1']*p['C3']; b=p['G2']*(p['C1']+p['C3']); c=p['G2']*p['G4']; n=[0,0,c]
    elif kind == 'HPF':
        a=p['C2']*p['C4']; b=p['C2']*(p['G1']+p['G3']); c=p['G1']*p['G3']; n=[a,0,0]
    else:
        a=p['C1']*p['C2']+p['C1']*p['C3']+p['C2']*p['C3']
        b=p['C2']*p['G1']+p['C1']*p['G3']+p['C3']*p['G1']+p['C2']*p['G3']+p['C2']*p['G4']
        c=p['G1']*p['G3']; n=[0,p['C2']*p['G4'],0]
    return np.array(n)/a, np.array([a,b,c])/a

def metrics(kind,p):
    n,d=coefficients(kind,p); w0=np.sqrt(d[2]); q=w0/d[1]
    return w0, w0/(2*np.pi), q

def response(kind,p,f):
    f=np.asarray(f)
    if np.any(~np.isfinite(f)) or np.any(f<=0): raise ValueError('frequency must be positive finite Hz')
    n,d=coefficients(kind,p); s=2j*np.pi*f
    return np.polyval(n,s)/np.polyval(d,s)

def admittances(kind,p,s):
    if kind=='LPF': return s*p['C1'],p['G2'],s*p['C3'],p['G4']
    if kind=='HPF': return p['G1'],s*p['C2'],p['G3'],s*p['C4']
    return s*p['C1']+p['G1'],s*p['C2'],s*p['C3']+p['G3'],p['G4']

def nodal_response(kind,p,f):
    """Independent 3-unknown solve for Va, Vy and Ix, Iin=1 A."""
    out=[]
    for freq in np.atleast_1d(f):
        y1,y2,y3,y4=admittances(kind,p,2j*np.pi*freq)
        matrix=np.array([[y1+y3,-y3,1],[-y3,y2+y3,0],[-y4,y4,1]],complex)
        out.append(np.linalg.solve(matrix,[1,0,0])[2])
    return np.array(out)

def sweep(kind,mode,value):
    p=BASE[kind].copy()
    if mode=='frequency':
        # Scale every conductance; fixed capacitances => w0 scales, Q unchanged.
        for key in p:
            if key.startswith('G'): p[key]*=value
    elif mode=='q':
        # Fixed w0: LPF G2G4 fixed; HPF G1G3 fixed; BPF G4 changes damping only.
        if kind=='LPF': p['G4']*=value; p['G2']/=value
        elif kind=='HPF': p['G1']*=value; p['G3']/=value
        else: p['G4']*=value
    else: raise ValueError(mode)
    return p
