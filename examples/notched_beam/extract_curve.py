# SPDX-License-Identifier: LGPL-2.1-or-later
"""Run with Abaqus Python: extract_curve.py cdp_coarse.odb (etc.)."""
from odbAccess import openOdb
from abaqusConstants import DOUBLE_PRECISION
import os, sys
path = os.path.abspath(sys.argv[1])
job = os.path.splitext(os.path.basename(path))[0]
model, size = job.split('_')
assert model in ['cdp', 'cdpm2']
width, top = {'coarse': (.02, [127,1127]), 'medium': (.01,[493,1493]), 'fine': (.005,[1465,3465])}[size]
odb = openOdb(path=path, readOnly=True)
def component(value):
    return float((value.dataDouble if value.precision==DOUBLE_PRECISION else value.data)[1])
rows = []
for step in odb.steps.values():
    for frame in step.frames:
        u = {v.nodeLabel:component(v) for v in frame.fieldOutputs['U'].values if v.nodeLabel in top}
        rf = {v.nodeLabel:component(v) for v in frame.fieldOutputs['RF'].values if v.nodeLabel in top}
        assert len(u)==len(rf)==2
        rows.append((-sum(u.values())/2*1000, -2*sum(rf.values())*.05/width/1000, float(frame.frameValue)))
odb.close()
folder=os.path.join(os.path.dirname(path),'results')
if not os.path.isdir(folder):os.makedirs(folder)
target=os.path.join(folder,job+'_curve.txt')
if os.path.exists(target):raise RuntimeError('Existing curve preserved: '+target)
with open(target,'w') as stream:
    stream.write('# displacement_mm full_beam_load_50mm_kN Riks_step_time\n')
    for row in rows:stream.write('%.17g %.17g %.17g\n'%row)
print(target)
