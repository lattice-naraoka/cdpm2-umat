"""Deterministic, independent one-element CDPM2 comparisons.

Abaqus: F=Q diag(exp(e)) Q.T. OOFEM: epsilon=Q diag(e) Q.T.
Q is constant throughout each test; all stretches commute, so there is no
rotation-history ambiguity. Prescribed lateral strains are NOT stress-free BCs.
Each Abaqus batch contains disconnected, fully prescribed one-element specimens
with identical material constants. Native OOFEM runs each specimen separately.
"""
from pathlib import Path
import argparse, json, hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
XYZ=np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)],float)
BASE=[30000.,.2,40.,3.,.525,.3,.01,.85,.08,.003,2.,1e-6,1.,.9,.02,.15,.0001,15.,1.,1e-5,1e-5,1e-9,20000.,.001]
def write(p,s):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_bytes(s.encode('utf-8'))

def paths():
    q=np.eye(3); r=np.array([[1,1,0],[-1,1,0],[0,0,np.sqrt(2)]])/np.sqrt(2)
    out=[]
    def add(name,points,rotation=q,category='nonlinear'):
        out.append(dict(id=name,knots=np.asarray(points,float).tolist(),Q=rotation.tolist(),category=category))
    for name,v,rot in [('elastic_x',[1e-6,0,0],q),('elastic_hydro',[1e-6]*3,q),('elastic_dev',[1e-6,-1e-6,0],q),('elastic_xy',[1e-6,-1e-6,0],r)]:
        add(name,[[0,0,0],v],rot,'elastic')
    add('tension_constrained',[[0,0,0],[.003,0,0]])
    add('compression_constrained',[[0,0,0],[-.01,0,0]])
    add('pure_shear_xy',[[0,0,0],[.003,-.003,0]],r)
    add('reversal_x',[[0,0,0],[.0005,0,0],[-.003,0,0],[.001,0,0],[-.006,0,0]])
    add('tension_lateral_contraction',[[0,0,0],[.003,-.0006,-.0006]])
    add('compression_lateral_expansion',[[0,0,0],[-.01,.002,.002]])
    add('biaxial_tension',[[0,0,0],[.001,.001,0]])
    add('biaxial_compression',[[0,0,0],[-.006,-.006,0]])
    add('triaxial_compression',[[0,0,0],[-.01,-.002,-.002]])
    add('hydro_tension',[[0,0,0],[.0004]*3])
    add('hydro_compression',[[0,0,0],[-.01]*3])
    add('mixed_sign',[[0,0,0],[.002,-.003,0]])
    add('tension_cycles',[[0,0,0],[.0003,0,0],[0,0,0],[.001,0,0],[0,0,0],[.003,0,0],[0,0,0]])
    add('compression_cycles',[[0,0,0],[-.001,0,0],[0,0,0],[-.004,0,0],[0,0,0],[-.01,0,0],[0,0,0]])
    add('shear_reversal',[[0,0,0],[.001,-.001,0],[-.001,.001,0],[.003,-.003,0],[-.003,.003,0]],r)
    add('tension_hold',[[0,0,0],[.001,0,0],[.001,0,0],[.001,0,0],[0,0,0]])
    add('compression_hold',[[0,0,0],[-.006,0,0],[-.006,0,0],[-.006,0,0],[0,0,0]])
    add('nonproportional_diagonal',[[0,0,0],[.0005,0,0],[.0005,-.003,0],[0,-.003,.001],[-.004,0,0]])
    add('alternating_12cycles',[[0,0,0]]+[[((-1)**i)*(.0002+.0001*i),0,0] for i in range(1,25)])
    rng=np.random.RandomState(1242026)
    add('seeded_diagonal',np.vstack([np.zeros(3),rng.uniform(-.002,.0005,(8,3))]))
    return out

VARIANTS={'baseline':{},'E_low':{0:20000},'E_high':{0:40000},'nu_low':{1:.1},'fc_high':{2:60},'ft_high':{3:4.5,13:1.35},'Hp_05':{6:.5},'As_5':{17:5},'As_30':{17:30},'Df_075':{7:.75},'Gf_half':{14:.01,15:.075},'Efc_double':{16:.0002},'length_2':{12:2}}

def material(p,reference_length):
    keys=['E','n','fc','ft','ecc','kinit','hp','dilation','ahard','bhard','chard','dhard','helem','ft1','wf1','wf','efc','asoft','bsoft']
    v=p[:19].copy();v[12]=reference_length
    return 'con2dpm 1 d 0 talpha 0 '+' '.join('%s %.16g'%(k,x) for k,x in zip(keys,v))+' yieldtol 1e-6 newtoniter 100 sratetype 0'

def histories(case,n):
    a=np.asarray(case['knots']); ts=np.linspace(0,1,len(a)); t=np.linspace(0,1,n+1)
    eig=np.column_stack([np.interp(t,ts,a[:,i]) for i in range(3)])
    q=np.asarray(case['Q']); eps=np.einsum('ij,tj,kj->tik',q,eig,q)
    F=np.einsum('ij,tj,kj->tik',q,np.exp(eig),q)
    return t,eps,F

def native(case,p,n,T):
    t,e,F=histories(case,n);xyz=XYZ*p[12];u=np.einsum('tij,nj->tni',e,xyz)
    s=['reference.out','Independent homogeneous one-element CDPM2 strain-path test',
       'StaticStructural nsteps %d deltat %.16g rtolf 1e-8 maxiter 100'%(n,T/n),'domain 3d',
       'OutputManager tstep_all dofman_all element_all','ndofman 8 nelem 1 ncrosssect 1 nmat 1 nbc 24 nltf 24 nic 0']
    for i,x in enumerate(xyz,1):s.append('node %d coords 3 %s bc 3 %d %d %d'%(i,' '.join(map(str,x)),3*i-2,3*i-1,3*i))
    # OOFEM local node ordering differs from Abaqus C3D8.
    # For the cubic specimens and fixed Q in this catalogue, the projected
    # span at damage initiation is L on the axes and sqrt(2)*L for XY shear.
    # It is derived geometrically, not fitted to stress or damage results.
    q=np.asarray(case['Q']);reference_length=p[12]*np.abs(q[:,0]).sum()
    s+=['lspace 1 nodes 8 5 8 7 6 1 4 3 2 mat 1 crossSect 1','SimpleCS 1',material(p,reference_length)]
    for i in range(1,25):s.append('BoundaryCondition %d loadTimeFunction %d prescribedvalue 1'%(i,i))
    for i in range(8):
        for j in range(3):
            k=3*i+j+1
            s.append('PiecewiseLinFunction %d nPoints %d t %d %s f(t) %d %s'%(k,n+1,n+1,' '.join('%.16g'%x for x in t*T),n+1,' '.join('%.16g'%x for x in u[:,i,j])))
    return '\n'.join(s)+'\n'

def batch(cases,p,n,T):
    s=['*Heading','Disconnected one-element strain-path verification specimens','*File Format,ASCII','*Part,name=SPECIMENS','*Node']
    for k,c in enumerate(cases):
        for i,x in enumerate(XYZ*p[12],1):s.append('%d,%.16g,%.16g,%.16g'%(8*k+i,x[0]+3*k*p[12],x[1],x[2]))
    s+=['*Element,type=C3D8,elset=ALL']
    for k in range(len(cases)):s.append('%d,%s'%(k+1,','.join(str(8*k+i) for i in range(1,9))))
    s+=['*Solid Section,elset=ALL,material=CDPM2_VISCOUS_1049',',','*End Part','*Assembly,name=ASSEMBLY','*Instance,name=P,part=SPECIMENS','*End Instance','*End Assembly','*Material,name=CDPM2_VISCOUS_1049','*Depvar','89','*User Material,constants=24,unsymm']
    for k in range(0,24,8):s.append(','.join('%.16g'%x for x in p[k:k+8]))
    bcs=[]
    for k,c in enumerate(cases):
        t,e,F=histories(c,n);u=np.einsum('tij,nj->tni',F-np.eye(3),XYZ*p[12])
        known={}
        for i in range(8):
            for j in range(3):
                vals=u[:,i,j]; key=tuple(np.round(vals,15))
                if key not in known:
                    name='A%d_%d'%(k+1,len(known));known[key]=name
                    s.append('*Amplitude,name='+name+',time=TOTAL TIME')
                    pairs=['%.16g,%.16g'%(a,b) for a,b in zip(t*T,vals)]
                    for kk in range(0,len(pairs),4):s.append(','.join(pairs[kk:kk+4]))
                bcs+=['*Boundary,amplitude='+known[key],'P.%d,%d,%d,1.'%(8*k+i+1,j+1,j+1)]
    # STATIC,DIRECT is a different procedure identifier (JSTEP(2)=2) and is
    # outside the released UMAT's Static General interface (JSTEP(2)=1).
    # Keep the comparison grid fixed using equal automatic min/max increments.
    s+=['*Step,name=PATH,nlgeom=YES,inc=%d'%(n+5),'*Static','%.16g,%.16g,%.16g,%.16g'%(T/n,T,T/n,T/n)]+bcs
    s+=['*Output,field,frequency=1','*Element Output','S,LE,SDV14,SDV37,SDV38,SDV39,SDV59,SDV60',
        '*Node File,frequency=999999','U','*El File,frequency=999999','S','*End Step']
    return '\n'.join(s)+'\n'

def main():
    a=argparse.ArgumentParser();a.add_argument('--pilot',action='store_true');a.add_argument('--n',type=int,default=240);a.add_argument('--out');a.add_argument('--group',action='append');a.add_argument('--case',action='append');args=a.parse_args()
    cases=paths()[:8] if args.pilot else paths(); variants={'baseline':{}} if args.pilot else VARIANTS
    if args.group:variants={k:v for k,v in variants.items() if k in args.group}
    if args.case:cases=[c for c in cases if c['id'] in args.case]
    if not variants or not cases:raise ValueError('Empty test selection')
    out=ROOT/(args.out or ('pilot' if args.pilot else 'suite')); manifest=dict(schema=1,n=args.n,T=1000.,mu=.001,cases=cases,groups=[])
    for group,changes in variants.items():
        p=BASE.copy()
        for i,v in changes.items():p[i]=v
        write(out/group/'test.inp',batch(cases,p,args.n,1000.))
        for case in cases:write(out/group/case['id']/'reference.in',native(case,p,args.n,1000.))
        manifest['groups'].append(dict(id=group,properties=p))
    manifest['criteria']=dict(elastic_max_scaled_stress=1e-4,nonlinear_max_scaled_stress=.01,nonlinear_rms_scaled_stress=.002,max_damage_difference=.01,max_strain_difference=1e-9,max_gp_stress_spread=1e-7,stress_scale='max(ft,max_abs_reference_tensor_over_history)')
    manifest['scope']='Coaxial prescribed strain histories, not stress-free uniaxial tests. Native OOFEM small strains matched to Abaqus logarithmic strains at prescribed time knots. UMAT uses automatic geometric characteristic length; OOFEM receives its independently derived projected cubic span, L for axial paths and sqrt(2)*L for XY shear. The length_2 group scales physical cube size to 2. No global tangent/Newton verification.'
    write(out/'manifest.json',json.dumps(manifest,indent=2)+'\n')
    print(out, len(cases)*len(variants),'independent specimen comparisons')
if __name__=='__main__':main()
