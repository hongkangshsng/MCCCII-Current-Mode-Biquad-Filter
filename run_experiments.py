"""Run from repository root: python run_experiments.py [--output results]."""
import argparse,csv,json,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.signal import freqs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.filters import BASE, coefficients, metrics, response, sweep

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',default='results'); args=parser.parse_args()
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    f=np.logspace(-1,5,2001); rows=[]
    for kind,base in BASE.items():
        for mode,values in [('baseline',[1]),('q',[0.25,1,4]),('frequency',[0.5,1,2])]:
            fig,axes=plt.subplots(2,1,figsize=(8,7),sharex=True)
            for value in values:
                p=base.copy() if mode=='baseline' else sweep(kind,mode,value)
                h=response(kind,p,f); w0,f0,q=metrics(kind,p)
                n,d=coefficients(kind,p); _,hs=freqs(n,d,2*np.pi*f)
                error=float(np.max(np.abs(h-hs)))
                label=f'factor={value:g}, f0={f0:.2f} Hz, Q={q:.4f}'
                axes[0].semilogx(f,20*np.log10(abs(h)),label=label)
                axes[1].semilogx(f,np.unwrap(np.angle(h))*180/np.pi,label=label)
                tag=f'{kind}_{mode}_{value:g}'
                np.savetxt(out/(tag+'.csv'),np.column_stack([f,abs(h),20*np.log10(abs(h)),np.unwrap(np.angle(h))*180/np.pi]),delimiter=',',header='frequency_hz,magnitude,gain_db,phase_deg',comments='')
                rows.append(dict(filter=kind,experiment=mode,factor=value,w0_rad_s=w0,f0_hz=f0,Q=q,parameters_si=json.dumps(p,sort_keys=True),scipy_max_abs_error=error))
            axes[0].set_ylabel('Current gain (dB)'); axes[1].set_ylabel('Phase (deg)'); axes[1].set_xlabel('Frequency (Hz)')
            for ax in axes: ax.grid(True,which='both',alpha=.25)
            axes[0].legend(fontsize=8); fig.suptitle(f'{kind}: {mode} | ideal transfer-function simulation (Rx=0)'); fig.tight_layout()
            fig.savefig(out/f'{kind}_{mode}.svg'); plt.close(fig)
    with (out/'summary.csv').open('w') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    lines=['# Generated numerical results / 新增理論模型數值結果','', '| Filter | Sweep | Factor | f0 (Hz) | Q |','|---|---|---:|---:|---:|']
    lines += [f"| {r['filter']} | {r['experiment']} | {r['factor']} | {r['f0_hz']:.6f} | {r['Q']:.6f} |" for r in rows]
    (out/'summary.md').write_text('\n'.join(lines)+'\n')
    (out/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,model='Ideal Rx=0; no measured data'),indent=2)+'\n')
    print(f'Generated {len(rows)} cases, 9 SVG plots and summary in {out}')
if __name__=='__main__': main()
