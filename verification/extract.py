from odbAccess import openOdb
from abaqusConstants import DOUBLE_PRECISION
from pathlib import Path
import numpy as np,sys,json
p=Path(sys.argv[1]);o=openOdb(str(p/'test.odb'),readOnly=True)
rows=[];times=[]; names=['S','LE','SDV14','SDV37','SDV38','SDV39','SDV59','SDV60']
def data(v):return np.atleast_1d(v.dataDouble if v.precision==DOUBLE_PRECISION else v.data).tolist()
for f in o.steps['PATH'].frames:
    if f.incrementNumber==0:continue
    result={}
    for name in names:
        vs=sorted(f.fieldOutputs[name].values,key=lambda v:(v.elementLabel,v.integrationPoint))
        result[name]=np.asarray([data(v) for v in vs])
    rows.append(result);times.append(f.frameValue)
o.close()
np.savez_compressed(str(p/'abaqus.npz'),time=np.asarray(times),**{name:np.asarray([r[name] for r in rows]) for name in names})
print('extracted',len(times),'frames')
