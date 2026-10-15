from pathlib import Path
from odbAccess import openOdb
import numpy as np,json,time,sys
p=Path(__file__).resolve().parent
assert len(sys.argv)==2,'Usage: abaqus python extract_curve.py cdp|cdpm2'
job=sys.argv[1]
assert job in ['cdp','cdpm2']
o=openOdb(str(p/(job+'.odb')),readOnly=True)
try:
    s=o.steps['PULLOUT'];ins=o.rootAssembly.instances['PART-1-1']
    sets={n:sorted({x.label for x in ins.nodeSets[n].nodes}) for n in ['REB_TOP','REB_BOTTOM']}
    refs=[]
    for target in [(50.,0.,0.),(0.,50.,0.)]:
        labels=[n.label for n in ins.nodes if max(abs(float(n.coordinates[i])-target[i]) for i in range(3))<1e-5]
        assert len(labels)==1,(target,labels)
        refs.append(labels[0])
    def hist(n,v):return np.array(list(s.historyRegions['Node PART-1-1.%s'%n].historyOutputs[v].data),dtype=float)
    t=hist(sets['REB_TOP'][0],'U3')[:,0]
    def avg(ns,v):return sum(np.interp(t,hist(n,v)[:,0],hist(n,v)[:,1]) for n in ns)/len(ns)
    top=avg(sets['REB_TOP'],'U3');bar=avg(sets['REB_BOTTOM'],'U3');conc=avg(refs,'U3')
    force=4*sum(np.interp(t,hist(n,'RF3')[:,0],hist(n,'RF3')[:,1]) for n in sets['REB_TOP'])/1000
    assert np.allclose(top,2*t,atol=3e-6)
    rows={float(tt):list(map(float,vals)) for tt,vals in zip(t,zip(t,top,bar,bar-conc,force))}
    a=np.array([rows[k] for k in sorted(rows)])
finally:o.close()
np.savetxt(str(p/(job+'_curve.csv')),a,delimiter=',',header='time,loaded_U3_mm,unloaded_U3_mm,relative_slip_2point_mm,full_load_kN',comments='')
d=dict(extracted_epoch=time.time(),reference_node_labels=refs,latest=a[-1].tolist(),maximum_so_far=a[int(a[:,4].argmax())].tolist(),points=len(a))
(p/(job+'_curve.json')).write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps(d,indent=2))
