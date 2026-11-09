from odbAccess import openOdb
import sys,json,csv
job=sys.argv[1]
o=openOdb(job+'.odb',readOnly=True)
try:
 ins=list(o.rootAssembly.instances.values())[0]
 names=['BAR_TOP','GRAPH_BAR_Z285','GRAPH_BAR_Z570','CONC_BOTTOM','CONC_TOP']
 labels={n:[v.label for v in ins.nodeSets[n].nodes] for n in names}
 rows=[];summary=[]
 for sn,step in o.steps.items():
  series={};common=None
  for name in names:
   for label in labels[name]:
    region_name='Node '+ins.name+'.'+str(label)
    region=step.historyRegions[region_name] if region_name in step.historyRegions else None
    for var in (['U3','RF3'] if name=='BAR_TOP' else ['U3']):
     key=(label,var)
     if key in series:continue
     h=region.historyOutputs[var] if region is not None and var in region.historyOutputs else None
     d=dict((h.data or []) if h is not None else []);series[key]=d
     common=set(d) if common is None else common.intersection(d)
  for t in sorted(common or []):
   means={n:sum(series[(l,'U3')][t] for l in labels[n])/len(labels[n]) for n in names}
   lo=means['GRAPH_BAR_Z285']-means['CONC_BOTTOM'];hi=means['GRAPH_BAR_Z570']-means['CONC_TOP']
   load=sum(series[(l,'RF3')][t] for l in labels['BAR_TOP'])/1000.
   rows.append([step.totalTime+t,lo,hi,(lo+hi)/2,load,means['BAR_TOP'],sn,t])
  summary.append(dict(step=sn,complete_points=len(common or [])))
 if not rows:
  print('No synchronized complete history samples yet; no CSV overwritten')
  sys.exit(0)
 with open(job+'_bond_slip_history.csv','w') as f:
  w=csv.writer(f);w.writerow(['time','lower_slip_mm','upper_slip_mm','mean_slip_mm','load_kN','loaded_displacement_mm','step','step_time']);w.writerows(rows)
 print(json.dumps(dict(counts={n:len(labels[n]) for n in names},steps=summary,last=rows[-1])))
finally:o.close()
