from pathlib import Path
import subprocess,json
import numpy as np
import argparse,tempfile
ap=argparse.ArgumentParser();ap.add_argument('--compiler',default='gfortran');args=ap.parse_args()
p=Path(__file__).resolve().parent;d=Path(tempfile.mkdtemp(prefix='cdpm2-spectral-'));print('Test directory:',d)
s=(p.parent/'src/cdpm2.f90').read_text().replace('public::endpoint_update,plastic_sensitivity','public::endpoint_update,plastic_sensitivity,split_sensitivity')
(d/'exposed.f90').write_text(s)
driver='''program probe
use endpoint_1076,only:split_sensitivity
use iso_fortran_env,only:real64
implicit none
real(real64)::s(6),ds(6,6),sp(6),sm(6),dp(6,6),dm(6,6)
integer::ios,rc
do
read(*,*,iostat=ios)s,ds
if(ios<0)exit
if(ios/=0)stop 1
call split_sensitivity(s,ds,sp,sm,dp,dm,rc,1e-9_real64)
if(rc/=0)stop 2
write(*,*)dp
enddo
end program
'''
(d/'driver.f90').write_text(driver);exe=d/'probe.exe'
r=subprocess.run([args.compiler,'-O2','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-ffree-line-length-none',str(d/'exposed.f90'),str(p/'test_stubs.f90'),str(p/'extra_stubs.f90'),str(d/'driver.f90'),'-o',str(exe)],cwd=d,capture_output=True,text=True);assert r.returncode==0,r.stderr
def pack(a):return a[[0,1,2,0,0,1],[0,1,2,1,2,2]]
def unpack(v):return np.array([[v[0],v[3],v[4]],[v[3],v[1],v[5]],[v[4],v[5],v[2]]])
rng=np.random.default_rng(20261106);cases=[];expected=[]
spectra=[[3,2,1],[-3,-2,-1],[3,.2,-1],[3,1e-10,-1],[3,-1e-10,-1],[3,0,-1],[2,2,-1],[0,0,0]]
for ev in spectra:
 for _ in range(20):
  q,_=np.linalg.qr(rng.normal(size=(3,3)));a=(q*np.array(ev))@q.T;v=pack(a);ds=rng.normal(size=(6,6));cases.append(np.r_[v,ds.flatten(order='F')])
  w,Q=np.linalg.eigh(a);tol=64*np.finfo(float).eps*max(abs(v).max(),np.finfo(float).tiny);slopes=np.where(w>tol,1.,np.where(w< -tol,0.,.5));F=np.zeros((3,3))
  for i in range(3):
   for j in range(3):F[i,j]=(max(w[i],0)-max(w[j],0))/(w[i]-w[j]) if abs(w[i]-w[j])>tol else .5*(slopes[i]+slopes[j])
  expected.append(np.column_stack([pack(Q@(F*(Q.T@unpack(ds[:,j])@Q))@Q.T) for j in range(6)]).flatten(order='F'))
r=subprocess.run([str(exe)],input='\n'.join(' '.join(format(x,'.17g') for x in row) for row in cases)+'\n',cwd=d,capture_output=True,text=True);assert r.returncode==0,r.stderr
(d/'output.txt').write_text(r.stdout);actual=np.loadtxt(d/'output.txt');expected=np.array(expected);err=np.linalg.norm(actual-expected,axis=1)/np.maximum(1,np.linalg.norm(expected,axis=1))
out={'cases':len(cases),'max_relative_error':float(err.max()),'threshold':1e-10,'pass':bool((err<1e-10).all()),'note':'Independent numpy eigendecomposition divided-difference check, including symmetric generalized slope at roundoff-zero; not a classical derivative claim at exact zero or an Abaqus test.'}
(d/'summary.json').write_text(json.dumps(out,indent=2));print(out)

assert out['pass'],out
