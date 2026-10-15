"""Read every integration point; never infer PASS from solver completion."""
from pathlib import Path
import re,json,sys,csv
import numpy as np
from generate import histories
NUM=r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?'
def native(p):
    txt=p.read_text(errors='replace')
    blocks=re.split(r'Output for time\s+('+NUM+')',txt)[1:]
    t=[];S=[];damage=[]
    for k in range(0,len(blocks),2):
        b=blocks[k+1]
        stresses=re.findall(r'stresses\s+((?:'+NUM+r'\s*){6})',b)
        dt=re.findall(r'damageTension\s+('+NUM+')',b);dc=re.findall(r'damageCompression\s+('+NUM+')',b)
        if len(stresses)!=8 or len(dt)!=8 or len(dc)!=8:raise ValueError('Incomplete OOFEM Gauss output')
        t.append(float(blocks[k]));S.append([[float(v) for v in s.split()] for s in stresses]);damage.append(np.array([dt,dc],float).T)
    return np.array(t),np.array(S)[:,:,[0,1,2,5,4,3]],np.array(damage)
def main():
    p=Path(sys.argv[1]);m=json.loads((p/'manifest.json').read_text());status=json.loads((p/'status.json').read_text());results=[];curves={}
    for g,gs in zip(m['groups'],status['groups']):
        f=p/g['id']/'abaqus.npz';a=np.load(f) if f.exists() else None
        for i,c in enumerate(m['cases']):
            row=dict(group=g['id'],case=c['id'],status='NOT_EVALUATED',native_completed=gs['native'][i]['completed'],abaqus_completed=gs.get('abaqus_completed',False))
            if a is None or not len(a['time']) or not row['native_completed']:results.append(row);continue
            t,s,d=native(p/g['id']/c['id']/'reference.out');nt=len(a['time'])
            # An aborted Abaqus step may append a duplicate-time terminal
            # diagnostic frame. It is not a new converged comparison point.
            while nt>1 and a['time'][nt-1]<=a['time'][nt-2]:nt-=1
            n=min(nt,len(t));inds=slice(8*i,8*i+8)
            ss=a['S'][:n,inds];le=a['LE'][:n,inds]
            # ODB frame times can be serialized in single precision even when
            # integration-point fields use double precision. Check against a
            # float32 serialization bound, NOT a relaxed material criterion.
            time_error=float(abs(t[:n]-a['time'][:n]).max())
            if time_error>8*np.finfo(np.float32).eps*m['T']:raise ValueError('Time grid mismatch: %s/%s max_error=%g'%(g['id'],c['id'],time_error))
            vs=np.concatenate([a['SDV59'][:n,inds],a['SDV60'][:n,inds]],axis=-1)
            _,eps,_=histories(c,m['n']);ev=eps[1:n+1][:,[0,1,2,0,0,1],[0,1,2,1,2,2]]
            # Abaqus strain field components use engineering shear convention;
            # ev is a true symmetric tensor, used unmodified for elasticity.
            le_expected=ev.copy();le_expected[:,3:]*=2
            scale=max(g['properties'][3],float(abs(s).max()));diff=ss-s[:n];crit=m['criteria']
            row.update(frames=n,stress_max=float(abs(diff).max()/scale),stress_rms=float(np.sqrt(np.mean(diff**2))/scale),damage_max=float(abs(vs-d[:n]).max()),strain_max=float(abs(le-le_expected[:,None,:]).max()),gp_spread=float(np.ptp(ss,axis=1).max()))
            lengths=a['SDV39'][:n,inds,0];active=lengths>0;expected_length=g['properties'][12]*np.abs(np.asarray(c['Q'])[:,0]).sum()
            row['length_max_abs_error']=float(abs(lengths[active]-expected_length).max()) if np.any(active) else 0.
            limit=crit['elastic_max_scaled_stress'] if c['category']=='elastic' else crit['nonlinear_max_scaled_stress']
            checks=[row['stress_max']<=limit,row['stress_rms']<=crit['nonlinear_rms_scaled_stress'],row['damage_max']<=crit['max_damage_difference'],row['strain_max']<=crit['max_strain_difference'],row['gp_spread']<=crit['max_gp_stress_spread']]
            row['status']=('PASS' if all(checks) else 'FAIL') if row['abaqus_completed'] and n==m['n'] else 'INCOMPLETE'
            # Independent elastic law gate, not just agreement of two codes.
            if c['category']=='elastic':
                E,nu=g['properties'][:2];G=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu));expected=2*G*ev;expected[:,:3]+=lam*ev[:,:3].sum(axis=1)[:,None]
                row['analytical_stress_max']=float(abs(ss-expected[:,None,:]).max()/scale)
                if row['analytical_stress_max']>crit['elastic_max_scaled_stress']:row['status']='FAIL'
            results.append(row);key=g['id']+'__'+c['id'];curves[key+'_S_umat']=ss.mean(axis=1);curves[key+'_S_oofem']=s[:n].mean(axis=1);curves[key+'_damage_umat']=vs.mean(axis=1);curves[key+'_damage_oofem']=d[:n].mean(axis=1);curves[key+'_strain']=ev
    (p/'results.json').write_text(json.dumps(results,indent=2))
    keys=list(dict.fromkeys(k for r in results for k in r))
    with (p/'results.csv').open('w',newline='') as f:w=csv.DictWriter(f,keys);w.writeheader();w.writerows(results)
    np.savez_compressed(p/'curves.npz',**curves)
    print({s:sum(r['status']==s for r in results) for s in sorted(set(r['status'] for r in results))})
if __name__=='__main__':main()
