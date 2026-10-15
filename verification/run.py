"""Run generated independent one-element verification specimens.

Requires installed Abaqus/Standard, a working Fortran compiler and native OOFEM.
The UMAT source is never modified. Existing results are never overwritten.
CPU count is one; manage any site-specific license scheduling externally.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,time

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('suite',type=Path)
    ap.add_argument('--oofem',required=True,help='Native OOFEM executable')
    ap.add_argument('--reference-source',required=True,type=Path,help='concretedpm2.C used to build OOFEM')
    ap.add_argument('--abaqus',default='abaqus')
    ap.add_argument('--umat',type=Path,default=Path(__file__).resolve().parents[1]/'src/cdpm2.f90')
    args=ap.parse_args();suite=args.suite.resolve();root=Path(__file__).resolve().parent
    m=json.loads((suite/'manifest.json').read_text());state=dict(phase='starting',groups=[])
    status=suite/'status.json'
    if status.exists():raise RuntimeError('Existing run: use a new generated suite directory')
    def save():state['updated']=time.time();status.write_text(json.dumps(state,indent=2))
    def call(cmd,d,log):
        with (d/log).open('w') as f:return subprocess.call(list(map(str,cmd)),cwd=str(d),stdout=f,stderr=subprocess.STDOUT)
    native=Path(args.oofem).resolve()
    if not native.is_file():raise FileNotFoundError(native)
    save();build=suite/'umat_library';build.mkdir()
    shutil.copy2(str(args.umat),str(build/'cdpm2.f90'))
    (build/'abaqus_v6.env').write_text("compile_fortran=compile_fortran+['-free']\n")
    if call([args.abaqus,'make','library=cdpm2.f90'],build,'build.log')!=0:raise RuntimeError('UMAT build failed; see build.log')
    libs=list(build.glob('*standardU*'))
    provenance=dict(umat_sha256=sha(args.umat),native_oofem_executable_sha256=sha(native),native_oofem_source_sha256=sha(args.reference_source),libraries={p.name:sha(p) for p in libs if p.is_file()})
    (suite/'provenance.json').write_text(json.dumps(provenance,indent=2))
    for g in m['groups']:
        d=suite/g['id'];entry=dict(id=g['id'],native=[]);state['groups'].append(entry);state['phase']='native_'+g['id'];save()
        for c in m['cases']:
            p=d/c['id'];code=call([native,'-f','reference.in'],p,'native.log')
            entry['native'].append(dict(id=c['id'],exit_code=code,completed=code==0 and 'ANALYSIS FINISHED' in (p/'native.log').read_text(errors='replace')));save()
        state['phase']='abaqus_'+g['id'];save()
        (d/'abaqus_v6.env').write_text("import os\nusub_lib_dir=%r\nfor k in ['CDPM2_GEOMETRY_TABLE','CDPM2_FIL_FILE']:os.environ.pop(k,None)\nos.environ['CDPM2_VISCO_LOG_LIMIT']='0'\n"%str(build))
        if (d/'test.odb').exists():raise RuntimeError('Refusing to overwrite '+str(d/'test.odb'))
        call([args.abaqus,'job=test','cpus=1','mp_mode=threads','output_precision=full','interactive'],d,'solver.log')
        sta=d/'test.sta';entry['abaqus_completed']=sta.exists() and 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' in sta.read_text(errors='replace');save()
        if (d/'test.odb').exists():entry['extract_exit']=call([args.abaqus,'python',root/'extract.py',d],d,'extract.log')
        save()
    state['phase']='finished';save()
if __name__=='__main__':main()
