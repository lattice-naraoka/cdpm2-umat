! SPDX-License-Identifier: LGPL-2.1-or-later
! OOFEM original portions: Copyright (C) 1993 - 2013 Borek Patzak
! Modified for a Fortran Abaqus UMAT, 2026-09-30 to 2026-10-12.
! Distributed without warranty. See ../LICENSE and ../THIRD_PARTY_NOTICES.md.
! Derived from OOFEM (LGPL-2.1-or-later), commit
! 1c1dbd2b4b8c95dcc26878bf5009d44ed535e1a6, floatmatrixf.h eig and
! structuralmaterial.C computePrincipalValDir. Independent diagnostic port.
module oofem_direction_154
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine principal(e,val,v,rc)
  real(real64),intent(in)::e(6) ! 11,22,33,engineering12,13,23
  real(real64),intent(out)::val(3),v(3,3)
  integer,intent(out)::rc
  real(real64)::m(3,3),total,ssum,tol,aa,si,co,tt,aij,aji,col(3)
  integer::i,j,k,ite,ii,jj
  rc=2;val=0;v=0
  if(.not.all(ieee_is_finite(e)))return
  do i=1,3;v(i,i)=1;end do
  m=0;do i=1,3;m(i,i)=e(i);end do
  m(1,2)=e(4)/2;m(2,1)=m(1,2)
  m(1,3)=e(5)/2;m(3,1)=m(1,3)
  m(2,3)=e(6)/2;m(3,2)=m(2,3)
  val=e(1:3);total=sum(abs(m))
  if(.not.ieee_is_finite(total))return
  if(total<=0)then
   rc=0;return
  end if
  tol=.1_real64**10;ite=0
  do
   ssum=0
   do j=2,3
   do i=1,j-1
    if(abs(m(i,j))/total<=tol)cycle
    ssum=ssum+abs(m(i,j))
    aa=atan2(2*m(i,j),val(i)-val(j))/2;si=sin(aa);co=cos(aa)
    do k=1,i-1
     tt=m(k,i);m(k,i)=co*tt+si*m(k,j);m(k,j)=-si*tt+co*m(k,j)
     tt=v(k,i);v(k,i)=co*tt+si*v(k,j);v(k,j)=-si*tt+co*v(k,j)
    end do
    tt=val(i);val(i)=co*tt+si*m(i,j);aij=-si*tt+co*m(i,j)
    tt=v(i,i);v(i,i)=co*tt+si*v(i,j);v(i,j)=-si*tt+co*v(i,j)
    do k=i+1,j-1
     tt=m(i,k);m(i,k)=co*tt+si*m(k,j);m(k,j)=-si*tt+co*m(k,j)
     tt=v(k,i);v(k,i)=co*tt+si*v(k,j);v(k,j)=-si*tt+co*v(k,j)
    end do
    tt=m(i,j);aji=co*tt+si*val(j);val(j)=-si*tt+co*val(j)
    tt=v(j,i);v(j,i)=co*tt+si*v(j,j);v(j,j)=-si*tt+co*v(j,j)
    do k=j+1,3
     tt=m(i,k);m(i,k)=co*tt+si*m(j,k);m(j,k)=-si*tt+co*m(j,k)
     tt=v(k,i);v(k,i)=co*tt+si*v(k,j);v(k,j)=-si*tt+co*v(k,j)
    end do
    val(i)=co*val(i)+si*aji;val(j)=-si*aij+co*val(j);m(i,j)=0
   end do
   end do
   ite=ite+1
   if(ite>50)return
   if(abs(ssum)/total<=tol)exit
  end do
  do ii=1,2
  do jj=1,2
   if(val(jj+1)>val(jj))then
    tt=val(jj);val(jj)=val(jj+1);val(jj+1)=tt
    col=v(:,jj);v(:,jj)=v(:,jj+1);v(:,jj+1)=col
   end if
  end do
  end do
  if(.not.all(ieee_is_finite(val)).or..not.all(ieee_is_finite(v)))return
  rc=0
 end subroutine
end module

module strain_direction_146
 use oofem_direction_154,only:principal,real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine maximum_strain_direction(e,n,lambda,gap,rc,nodes)
 real(real64),intent(in)::e(6)
 real(real64),intent(in),optional::nodes(3,8)
 real(real64),intent(out)::n(3),lambda,gap
 integer,intent(out)::rc
 real(real64)::v(3,3),ev(3),b(3,3),a(3,3),el(6),edge(3),scale,len
 integer::i
 n=0;lambda=0;gap=0;rc=2
 if(.not.all(ieee_is_finite(e)))return
 b=0;do i=1,3;b(i,i)=1;end do
 if(present(nodes))then
  if(.not.all(ieee_is_finite(nodes)))return
  edge=nodes(:,2)-nodes(:,1);scale=sqrt(sum(edge**2))
  if(.not.ieee_is_finite(scale).or.scale<=tiny(scale))return
  b(:,1)=edge/scale
  edge=nodes(:,4)-nodes(:,1);len=sqrt(sum(edge**2))
  if(.not.ieee_is_finite(len).or.len<=tiny(len))return
  edge=edge/len;edge=edge-dot_product(edge,b(:,1))*b(:,1)
  len=sqrt(sum(edge**2));if(len<=1e-12_real64)return
  b(:,2)=edge/len
  b(:,3)=[b(2,1)*b(3,2)-b(3,1)*b(2,2),b(3,1)*b(1,2)-b(1,1)*b(3,2),b(1,1)*b(2,2)-b(2,1)*b(1,2)]
 end if
 a=0;do i=1,3;a(i,i)=e(i);end do
 a(1,2)=e(4)/2;a(2,1)=a(1,2);a(1,3)=e(5)/2;a(3,1)=a(1,3);a(2,3)=e(6)/2;a(3,2)=a(2,3)
 a=matmul(transpose(b),matmul(a,b))
 el=[a(1,1),a(2,2),a(3,3),2*a(1,2),2*a(1,3),2*a(2,3)]
 call principal(el,ev,v,rc);if(rc/=0)return
 n=matmul(b,v(:,1));lambda=ev(1);gap=ev(1)-ev(2)
 end subroutine
end module

module geometry_length_143
  use iso_fortran_env, only: real64
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  subroutine projected_span(x,n,span,rc)
    real(real64),intent(in)::x(:,:),n(3)
    real(real64),intent(out)::span
    integer,intent(out)::rc
    real(real64)::p(size(x,2))
    span=0;rc=1
    if(size(x,1)/=3.or.size(x,2)<1) return
    if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(n))) return
    p=matmul(n,x)
    if(.not.all(ieee_is_finite(p))) return
    span=maxval(p)-minval(p)
    if(.not.ieee_is_finite(span)) return
    rc=0
  end subroutine
end module

module length_state_145
  use iso_fortran_env, only: real64
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  pure subroutine select_length(oldle,old_damage,tensile_active,kappa,e0,tol,helem,span,mean_size,trialle,rc)
    real(real64),intent(in)::oldle,old_damage(2),kappa,e0,tol,helem,span,mean_size
    logical,intent(in)::tensile_active
    real(real64),intent(out)::trialle
    integer,intent(out)::rc
    real(real64)::candidate
    trialle=oldle;rc=1
    if(.not.all(ieee_is_finite([oldle,old_damage,kappa,e0,tol,helem]))) return
    if(oldle<0.or.any(old_damage<0).or.any(old_damage>1).or.e0<=0.or.tol<0.or.tol>=1) return
    rc=0
    if(.not.tensile_active) return
    if(kappa<=e0*(1-tol)) return
    candidate=oldle
    if(helem>0) then
      candidate=helem
    else if(all(old_damage==0)) then
      candidate=span
      if(candidate==0) candidate=mean_size
    else if(oldle==0) then
      candidate=mean_size
    end if
    ! Deliberate fail-safe vs upstream invalid geometry: no invalid state commit.
    if(.not.ieee_is_finite(candidate)) then
      rc=2;return
    end if
    if(candidate<=0) then
      rc=2;return
    end if
    trialle=candidate
  end subroutine
end module

module cdpm2_functions
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: hardening_one, hardening_two, yield_value
  public :: flow_gradient, ductility, hardening_rate, regular_residual, regular_jacobian
contains
  pure function flow_gradient(sig,rho,k,fc,ft,ecc,q0,hp,dilation) result(dg)
    ! Literal computeDGDInv translation; nonassociated potential, not df/dsigma.
    real(real64), intent(in) :: sig,rho,k,fc,ft,ecc,q0,hp,dilation
    real(real64) :: dg(2),q1,q2,m,ag,bg,mq,a,b,logmq
    q1=hardening_one(k,q0,hp)
    q2=hardening_two(k,hp)
    m=3*(fc**2-ft**2)/(fc*ft)*ecc/(ecc+1)
    ag=3*ft*q1*q2/fc+m*q1**2/2
    bg=q1*q2*(1+ft/fc)/3 / &
       (log(ag)+log(dilation+1)-log(2*dilation-1)-log(3*q1*q2+m*q1**2/2))
    logmq=log(ag)+(sig-ft*q1*q2/3)/(fc*bg)
    ! Reject numerically unusable exponential trials; never clamp the flow.
    ! Reserve headroom for squared flow norm and subsequent residual products.
    if (.not.ieee_is_finite(logmq)) then
      dg=ieee_value(0.0_real64,ieee_quiet_nan)
      return
    end if
    if (logmq>log(huge(1.0_real64))/4) then
      dg=ieee_value(0.0_real64,ieee_quiet_nan)
      return
    end if
    mq=exp(logmq)
    b=sig/fc+rho/(fc*sqrt(6.0_real64))
    a=(1-q1)*b**2+sqrt(1.5_real64)*rho/fc
    dg(1)=4*(1-q1)/fc*a*b+q1**2*mq/fc
    dg(2)=a/(sqrt(6.0_real64)*fc)*(4*(1-q1)*b+6)+m*q1**2/(sqrt(6.0_real64)*fc)
  end function

  pure real(real64) function ductility(sig,theta,fc,ah,bh,ch,dh) result(d)
    real(real64), intent(in) :: sig,theta,fc,ah,bh,ch,dh
    real(real64) :: x,eh,fh
    x=-(sig+fc/3)/fc
    if (x<0) then
      eh=bh-dh
      fh=(bh-dh)*ch/(ah-bh)
      d=(eh*exp(x/fh)+dh)/(2*cos(theta))**2
    else
      d=(ah+(bh-ah)*exp(-x/ch))/(2*cos(theta))**2
    end if
  end function

  pure real(real64) function hardening_rate(dg,d) result(rate)
    real(real64), intent(in) :: dg(2),d
    rate=sqrt(dg(1)**2/3+dg(2)**2)/d
  end function

  pure function regular_residual(x,trial,kold,bulk,shear,theta,fc,ft,ecc,q0,hp,dilation,ah,bh,ch,dh) result(r)
    ! Unknowns: mean stress, deviator norm, kappa, plastic multiplier.
    ! Residual definition only. No Newton iteration or vertex return yet.
    real(real64), intent(in) :: x(4),trial(2),kold,bulk,shear,theta,fc,ft,ecc,q0,hp,dilation,ah,bh,ch,dh
    real(real64) :: r(4),dg(2),d
    dg=flow_gradient(x(1),x(2),x(3),fc,ft,ecc,q0,hp,dilation)
    if (.not.all(ieee_is_finite(dg))) then
      r=ieee_value(0.0_real64,ieee_quiet_nan)
      return
    end if
    d=ductility(x(1),theta,fc,ah,bh,ch,dh)
    r(1)=x(1)-trial(1)+bulk*x(4)*dg(1)
    r(2)=x(2)-trial(2)+2*shear*x(4)*dg(2)
    r(3)=-x(3)+kold+x(4)*hardening_rate(dg,d)
    r(4)=yield_value(x(1),x(2),theta,x(3),fc,ft,ecc,q0,hp)
  end function

  pure real(real64) function hardening_one(k, q0, hp) result(q)
    real(real64), intent(in) :: k, q0, hp
    if (k <= 0.0_real64) then
      q = q0
    else if (k < 1.0_real64) then
      q = (1-q0-hp)*k**3 - (3*(1-q0)-3*hp)*k**2 &
          + (3*(1-q0)-2*hp)*k + q0
    else
      q = 1.0_real64
    end if
  end function

  pure real(real64) function hardening_two(k, hp) result(q)
    real(real64), intent(in) :: k, hp
    q = 1.0_real64
    if (k >= 1.0_real64) q = 1 + (k-1)*hp
  end function

  pure real(real64) function yield_value(sig, rho, theta, k, fc, ft, ecc, q0, hp) result(f)
    ! OOFEM invariant convention: sig=mean stress, rho=deviator norm,
    ! theta in [0,pi/3]. Valid material parameters required by caller.
    real(real64), intent(in) :: sig, rho, theta, k, fc, ft, ecc, q0, hp
    real(real64) :: q1, q2, r, a, m, c
    q1 = hardening_one(k,q0,hp)
    q2 = hardening_two(k,hp)
    c = cos(theta)
    r = (4*(1-ecc**2)*c**2+(2*ecc-1)**2) / &
        (2*(1-ecc**2)*c+(2*ecc-1)*sqrt(4*(1-ecc**2)*c**2+5*ecc**2-4*ecc))
    m = 3*(fc**2-ft**2)/(fc*ft)*ecc/(ecc+1)
    a = (1-q1)*(sig/fc+rho/(sqrt(6.0_real64)*fc))**2 + sqrt(1.5_real64)*rho/fc
    f = a**2+q1**2*q2*m*(sig/fc+rho*r/(sqrt(6.0_real64)*fc))-q1**2*q2**2
  end function
  pure subroutine regular_jacobian(x,bulk,shear,theta,fc,ft,ecc,q0,hp,dilation,ah,bh,ch,dh,jac,ok)
    ! Derivative of existing regular_residual w.r.t. (mean,rho,kappa,lambda).
    ! Trial stress, old kappa and trial theta are fixed. NOT Abaqus DDSDDE.
    ! At kappa=0/1 use admissible right-branch derivatives; do not average laws.
    real(real64),intent(in)::x(4),bulk,shear,theta,fc,ft,ecc,q0,hp,dilation,ah,bh,ch,dh
    real(real64),intent(out)::jac(4,4)
    logical,intent(out)::ok
    real(real64)::q,v,qp,vp,m,b,a,ag,bg,top,bot,z,num,mq,logmq,term
    real(real64)::dq(3),dv(3),db(3),da(3),dag(3),dbg(3),dz(3),dtop(3),dlog(3),dnum(3)
    real(real64)::dg(2),ddg(2,3),dd(3),dn(3),drate(3),df(3),dt(3),dterm(3)
    real(real64)::ddilat,n,rate,xx,eh,fh,c,shape,t
    ok=.false.;jac=0
    dg=flow_gradient(x(1),x(2),x(3),fc,ft,ecc,q0,hp,dilation)
    if(.not.all(ieee_is_finite(dg)))return
    q=hardening_one(x(3),q0,hp);v=hardening_two(x(3),hp);qp=0;vp=0
    if(x(3)>=0.and.x(3)<1)qp=3*(1-q0-hp)*x(3)**2-2*(3*(1-q0)-3*hp)*x(3)+(3*(1-q0)-2*hp)
    if(x(3)>=1)vp=hp
    dq=[0._real64,0._real64,qp];dv=[0._real64,0._real64,vp]
    m=3*(fc**2-ft**2)/(fc*ft)*ecc/(ecc+1)
    b=x(1)/fc+x(2)/(sqrt(6._real64)*fc)
    db=[1/fc,1/(sqrt(6._real64)*fc),0._real64]
    a=(1-q)*b**2+sqrt(1.5_real64)*x(2)/fc
    da=2*(1-q)*b*db-b**2*dq;da(2)=da(2)+sqrt(1.5_real64)/fc
    ag=3*ft*q*v/fc+m*q**2/2
    dag=3*ft/fc*(v*dq+q*dv)+m*q*dq
    z=3*q*v+m*q**2/2;dz=3*(v*dq+q*dv)+m*q*dq
    top=q*v*(1+ft/fc)/3;dtop=(v*dq+q*dv)*(1+ft/fc)/3
    bot=log(ag)+log(dilation+1)-log(2*dilation-1)-log(z)
    bg=top/bot;dbg=(dtop-bg*(dag/ag-dz/z))/bot
    num=x(1)-ft*q*v/3;dnum=-ft*(v*dq+q*dv)/3;dnum(1)=dnum(1)+1
    logmq=log(ag)+num/(fc*bg);mq=exp(logmq)
    dlog=dag/ag+dnum/(fc*bg)-num*dbg/(fc*bg**2)
    ddg(1,:)=4/fc*((1-q)*(b*da+a*db)-a*b*dq)+(2*q*mq*dq+q**2*mq*dlog)/fc
    term=4*(1-q)*b+6;dterm=4*((1-q)*db-b*dq)
    ddg(2,:)=(term*da+a*dterm+2*m*q*dq)/(sqrt(6._real64)*fc)
    ddilat=ductility(x(1),theta,fc,ah,bh,ch,dh);dd=0;xx=-(x(1)+fc/3)/fc
    if(xx<0)then
      eh=bh-dh;fh=(bh-dh)*ch/(ah-bh)
      dd(1)=-eh*exp(xx/fh)/(fc*fh*(2*cos(theta))**2)
    else
      dd(1)=(bh-ah)*exp(-xx/ch)/(fc*ch*(2*cos(theta))**2)
    end if
    n=sqrt(dg(1)**2/3+dg(2)**2)
    if(n<=tiny(1._real64).or.ddilat<=0)return
    dn=(dg(1)*ddg(1,:)/3+dg(2)*ddg(2,:))/n
    rate=n/ddilat;drate=dn/ddilat-n*dd/ddilat**2
    c=cos(theta)
    shape=(4*(1-ecc**2)*c**2+(2*ecc-1)**2)/ &
      (2*(1-ecc**2)*c+(2*ecc-1)*sqrt(4*(1-ecc**2)*c**2+5*ecc**2-4*ecc))
    t=x(1)/fc+x(2)*shape/(sqrt(6._real64)*fc)
    dt=[1/fc,shape/(sqrt(6._real64)*fc),0._real64]
    df=2*a*da+m*(2*q*v*t*dq+q**2*t*dv+q**2*v*dt)-2*q*v**2*dq-2*q**2*v*dv
    jac(1,1:3)=bulk*x(4)*ddg(1,:);jac(1,1)=jac(1,1)+1;jac(1,4)=bulk*dg(1)
    jac(2,1:3)=2*shear*x(4)*ddg(2,:);jac(2,2)=jac(2,2)+1;jac(2,4)=2*shear*dg(2)
    jac(3,1:3)=x(4)*drate;jac(3,3)=jac(3,3)-1;jac(3,4)=rate
    jac(4,1:3)=df
    ok=all(ieee_is_finite(jac))
  end subroutine

end module

module cdpm2_return
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  use cdpm2_functions
  implicit none
  private
  public :: material_parameters, regular_return, regular_return_seed
  type material_parameters
    real(real64) :: bulk=19940.1883333333_real64, shear=14955.14125_real64
    real(real64) :: fc=35.892339_real64,ft=2.71644205_real64,ecc=.525_real64
    real(real64) :: q0=.3_real64,hp=.5_real64,dilation=.85_real64
    real(real64) :: ah=.08_real64,bh=.003_real64,ch=2.0_real64,dh=1e-6_real64
  end type
contains
  subroutine regular_return(p,trial,theta,kold,state,rc,iterations,maxiter)
    ! Independent1098: same endpoint equations and two bounded Newton attempts.
    ! The fallback kappa is estimated from the trial yield surface with q1=1.
    ! This is only an initial guess; old history, trial stress and time stay fixed.
    ! Successful first attempts are unchanged. No strain/time subdivisions.
    type(material_parameters),intent(in)::p
    real(real64),intent(in)::trial(2),theta,kold
    real(real64),intent(inout)::state(4)
    integer,intent(out)::rc,iterations
    integer,intent(in),optional::maxiter
    integer::second_iterations,limit
    real(real64)::seed(4),f1,qa,qb,disc,qstar,guess
    limit=60
    if(present(maxiter))limit=maxiter
    call regular_return_seed(p,trial,theta,kold,state,rc,iterations,limit)
    !1119: rc4 is a failed Newton trial after validated entry, not invalid input.
    if(rc==1.or.rc==3.or.rc==4)then
      guess=1._real64
      f1=yield_value(trial(1),trial(2),theta,1._real64,p%fc,p%ft,p%ecc,p%q0,p%hp)
      if(ieee_is_finite(f1).and.f1>0.and.p%hp>0)then
        ! At q1=1: F=qa+qb*q2-q2**2, q2=1+(kappa-1)*hp.
        ! qa=(sqrt(3/2)*rho/fc)**2; qb=F(q2=1)-qa+1.
        qa=1.5_real64*(trial(2)/p%fc)**2
        qb=f1-qa+1._real64
        disc=norm2([qb,2*sqrt(qa)])
        if(qb>=0)then
          qstar=(qb+disc)/2
        else
          qstar=2*qa/(disc-qb)
        endif
        if(ieee_is_finite(qstar).and.qstar>1)then
          guess=1+(qstar-1)/p%hp
          if(.not.ieee_is_finite(guess))guess=1._real64
        endif
      endif
      ! The initial guess must remain admissible after pre-peak hardening.
      seed=[trial,max(kold,guess),0._real64]
      call regular_return_seed(p,trial,theta,kold,state,rc,second_iterations,limit,seed)
      iterations=iterations+second_iterations
    endif
  end subroutine

  subroutine regular_return_seed(p,trial,theta,kold,state,rc,iterations,maxiter,seed)
    type(material_parameters), intent(in) :: p
    real(real64), intent(in) :: trial(2),theta,kold
    real(real64), intent(inout) :: state(4)
    integer, intent(out) :: rc,iterations
    integer, intent(in), optional :: maxiter
    real(real64), intent(in), optional :: seed(4)
    real(real64) :: y(4),yp(4),ym(4),r(4),rp(4),rm(4),jac(4,4),dy(4),scale(4),h,a,n0
    real(real64)::colscale(4),jj(4,4),normal(4,4),mat(4,4),rhs(4),vec(4),step(4)
    real(real64)::yield_scale
    logical::active(4),found
    integer :: damp,aset,limit,j,ls
    logical :: ok
    ! state is committed only on rc=0. rc1=iteration/line search failure,
    ! rc2=invalid input/nonfinite arithmetic, rc3=singular Jacobian.
    rc=2
    iterations=0
    if (.not.all(ieee_is_finite([trial,theta,kold,p%bulk,p%shear,p%fc,p%ft,p%ecc, &
       p%q0,p%hp,p%dilation,p%ah,p%bh,p%ch,p%dh]))) return
    if (min(p%bulk,p%shear,p%ft,p%q0,p%dh,p%ch)<=0 .or. p%fc<=p%ft) return
    if (p%ecc<=.5_real64.or.p%ecc>=1.or.p%dilation<=.5_real64.or.p%dilation>=1) return
    if (p%ah<=p%bh.or.p%bh<=p%dh.or.p%q0>1.or.p%hp<0) return
    if (trial(2)<0.or.kold<0.or.theta<0.or.theta>acos(-1.0_real64)/3) return
    limit=60
    if (present(maxiter)) limit=maxiter
    if (yield_value(trial(1),trial(2),theta,kold,p%fc,p%ft,p%ecc,p%q0,p%hp)<=1e-10_real64) then
      state=[trial,kold,0.0_real64]
      rc=0
      return
    end if
    !1098: normalize ONLY yield convergence test by committed hardening scale.
    ! Keep residual, Jacobian, Newton direction and line-search merit unchanged.
    ! This explicitly changes the yield stopping norm, not the material equations.
    yield_scale=max(1._real64,(hardening_one(kold,p%q0,p%hp)*hardening_two(kold,p%hp))**2)
    scale=[p%fc,p%fc,1.0_real64,p%fc**2/p%bulk]
    y=[trial,kold,0.0_real64]/scale
    if(present(seed))y=seed/scale
    rc=1
    do iterations=0,limit
      call residual(y,r,ok)
      if (.not.ok) then
        !1103: valid entry, but unusable Newton trial (e.g. exponential guard).
        ! Distinct from invalid input/state rc2; endpoint maps rc4 to cutback.
        ! Do not clamp constitutive flow or accept this failed trial.
        rc=4
        return
      end if
      n0=maxval(abs(r))
      !861: tighter scaled regular residual; equations/activation/backtracking unchanged.
      if (maxval(abs(r(1:3)))<1e-12_real64.and.abs(r(4))/yield_scale<1e-12_real64) then
        state=y*scale
        rc=0
        return
      end if
      if (iterations==limit) return
      call regular_jacobian(y*scale,p%bulk,p%shear,theta,p%fc,p%ft,p%ecc,p%q0,p%hp, &
        p%dilation,p%ah,p%bh,p%ch,p%dh,jac,ok)
      if(.not.ok)then
        rc=2;return
      end if
      do j=1,4
        jac(:,j)=jac(:,j)*scale(j)
      end do
      jac(1:2,:)=jac(1:2,:)/p%fc
      call solve4(jac,-r,dy,ok)
      if (.not.ok) then
        rc=3
        return
      end if
      a=1
      do ls=0,24
        yp=y+a*dy
        if (yp(2)>=0.and.yp(3)>=kold.and.yp(4)>=0) then
          call residual(yp,rp,ok)
          if (ok) then
            if (maxval(abs(rp))<(1-1e-4_real64*a)*n0) exit
          end if
        end if
        a=a/2
      end do
      if (ls>24) then
        ! 892: independent constrained least-squares search direction, used
        ! ONLY after the unchanged Newton line search has failed. The material
        ! residual, admissible final state and convergence tolerance are unchanged.
        ! This is our numerical strategy, not a copied CDPM2 constitutive formula.
        colscale=max(sqrt(sum(jac**2,dim=1)),1.0_real64)
        do j=1,4
          jj(:,j)=jac(:,j)/colscale(j)
        enddo
        found=.false.
        do damp=-6,6,2
          normal=matmul(transpose(jj),jj)
          rhs=-matmul(transpose(jj),r)
          do j=1,4
            normal(j,j)=normal(j,j)+10.0_real64**damp
          enddo
          active=.false.
          do aset=1,4
            mat=normal;vec=rhs
            do j=1,4
              if(active(j))then
                mat(j,:)=0;mat(:,j)=0;mat(j,j)=1;vec(j)=0
              endif
            enddo
            call solve4(mat,vec,step,ok)
            if(.not.ok)exit
            dy=step/colscale
            if(y(2)<=0.and.dy(2)<0.and..not.active(2))then
              active(2)=.true.;cycle
            endif
            if(y(3)<=kold.and.dy(3)<0.and..not.active(3))then
              active(3)=.true.;cycle
            endif
            if(y(4)<=0.and.dy(4)<0.and..not.active(4))then
              active(4)=.true.;cycle
            endif
            exit
          enddo
          if(.not.ok)cycle
          a=1
          do ls=0,24
            yp=y+a*dy
            if(yp(2)>=0.and.yp(3)>=kold.and.yp(4)>=0)then
              call residual(yp,rp,ok)
              if(ok)then
                if(maxval(abs(rp))<(1-1e-4_real64*a)*n0)then
                  found=.true.;exit
                endif
              endif
            endif
            a=a/2
          enddo
          if(found)exit
        enddo
        if(.not.found)return
      endif
      y=yp
    end do
  contains
    subroutine residual(z,f,finite)
      real(real64), intent(in) :: z(4)
      real(real64), intent(out) :: f(4)
      logical, intent(out) :: finite
      f=regular_residual(z*scale,trial,kold,p%bulk,p%shear,theta,p%fc,p%ft,p%ecc, &
         p%q0,p%hp,p%dilation,p%ah,p%bh,p%ch,p%dh)
      f(1:2)=f(1:2)/p%fc
      finite=all(ieee_is_finite(f))
    end subroutine
  end subroutine

  subroutine solve4(a,b,x,ok)
    real(real64), intent(in) :: a(4,4),b(4)
    real(real64), intent(out) :: x(4)
    logical, intent(out) :: ok
    real(real64) :: aug(4,5),row(5),factor
    integer :: i,j,k
    aug(:,1:4)=a
    aug(:,5)=b
    ok=.false.
    x=0
    do i=1,4
      k=i-1+maxloc(abs(aug(i:4,i)),dim=1)
      if (abs(aug(k,i))<1e-14_real64*max(1.0_real64,maxval(abs(a)))) return
      row=aug(i,:)
      aug(i,:)=aug(k,:)
      aug(k,:)=row
      do j=i+1,4
        factor=aug(j,i)/aug(i,i)
        aug(j,i:5)=aug(j,i:5)-factor*aug(i,i:5)
      end do
    end do
    do i=4,1,-1
      x(i)=(aug(i,5)-sum(aug(i,i+1:4)*x(i+1:4)))/aug(i,i)
    end do
    ok=all(ieee_is_finite(x))
  end subroutine
end module

module cdpm2_vertex
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  use cdpm2_functions
  use cdpm2_return, only: material_parameters
  implicit none
  private
  public :: vertex_kappa, vertex_return
contains
  pure real(real64) function vertex_kappa(p,trial,kold,s) result(k)
    type(material_parameters), intent(in) :: p
    real(real64), intent(in) :: trial(2),kold,s
    real(real64) :: dep,d
    dep=sqrt(((trial(1)-s)/p%bulk)**2/9+(trial(2)/(2*p%shear))**2)
    ! Retains OOFEM's pi/3 vertex convention, including tensile vertex.
    d=ductility(s,acos(-1.0_real64)/3,p%fc,p%ah,p%bh,p%ch,p%dh)
    k=kold+dep/d
  end function

  subroutine vertex_return(p,trial,kold,state,rc,iterations,maxiter)
    ! Requires validated material parameters (same domain as regular_return).
    ! state=(mean,rho,kappa); no invented plastic multiplier at the apex.
    ! rc0=accepted, rc1=regular branch required, rc2=invalid/nonfinite,
    ! rc3=iteration budget exhausted. No output-state change on failure.
    type(material_parameters), intent(in) :: p
    real(real64), intent(in) :: trial(2),kold
    real(real64), intent(inout) :: state(3)
    integer, intent(out) :: rc,iterations
    integer, intent(in), optional :: maxiter
    real(real64) :: outside,inside,s,f,k,fo,fi,dg(2),ratio_p,ratio_t
    integer :: limit
    rc=2
    iterations=0
    if (.not.all(ieee_is_finite([trial,kold]))) return
    if (trial(2)<0.or.kold<0) return
    rc=1
    if (trial(1)==0) return
    if (trial(1)<0.and.kold>=1) return
    outside=trial(1)
    inside=0
    call evaluate(outside,fo,k)
    call evaluate(inside,fi,k)
    if (.not.all(ieee_is_finite([fo,fi,k]))) then
      rc=2
      return
    end if
    ! Bracket must start outside with origin inside; do not force a bracket.
    if (fo<=0.or.fi>=0) return
    limit=250
    if (present(maxiter)) limit=maxiter
    rc=3
    do iterations=1,limit
      s=inside+(outside-inside)/2
      call evaluate(s,f,k)
      if (.not.all(ieee_is_finite([f,k]))) then
        rc=2
        return
      end if
      if (f<=0) then
        inside=s
      else
        outside=s
      end if
      ! Independent599: resolve the scalar bracket and select its inside end.
      ! The rounded midpoint can stick at the outside end (f>0), so requiring
      ! that midpoint to be inside can falsely exhaust all250 iterations.
      if (abs(outside-inside)<=16*epsilon(s)*max(p%ft,abs(s))) then
        s=inside
        call evaluate(s,f,k)
        if (.not.all(ieee_is_finite([f,k]))) then
          rc=2
          return
        end if
      end if
      if (f<=0.and.abs(f)<1e-10_real64.and. &
          abs(outside-inside)<=16*epsilon(s)*max(p%ft,abs(s))) then
        dg=flow_gradient(s,0.0_real64,k,p%fc,p%ft,p%ecc,p%q0,p%hp,p%dilation)
        if (.not.all(ieee_is_finite(dg))) then
          rc=2
          return
        end if
        rc=1
        if (abs(dg(1))<tiny(1.0_real64).or.trial(1)==s) return
        ! 2G/K equals OOFEM 3(1-2nu)/(1+nu).
        ratio_p=dg(2)/dg(1)*(2*p%shear/p%bulk)
        ratio_t=trial(2)/(trial(1)-s)
        if ((trial(1)>0.and.ratio_p>=ratio_t).or.(trial(1)<0.and.ratio_p<=ratio_t)) then
          state=[s,0.0_real64,k]
          rc=0
        end if
        return
      end if
    end do
  contains
    subroutine evaluate(sig,value,kap)
      real(real64), intent(in) :: sig
      real(real64), intent(out) :: value,kap
      kap=vertex_kappa(p,trial,kold,sig)
      value=yield_value(sig,0.0_real64,0.0_real64,kap,p%fc,p%ft,p%ecc,p%q0,p%hp)
    end subroutine
  end subroutine
end module

module cdpm2_plastic
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  use cdpm2_functions
  use cdpm2_return
  use cdpm2_vertex
  implicit none
  private
  public :: elastic_stress, elastic_strain, invariants, plastic_update
contains
  pure function elastic_stress(p,e) result(s)
    type(material_parameters), intent(in) :: p
    real(real64), intent(in) :: e(6)
    real(real64) :: s(6),v
    ! Abaqus order 11,22,33,12,13,23; strain shears are engineering.
    v=sum(e(1:3))
    s(1:3)=2*p%shear*e(1:3)+(p%bulk-2*p%shear/3)*v
    s(4:6)=p%shear*e(4:6)
  end function
  pure function elastic_strain(p,s) result(e)
    type(material_parameters), intent(in) :: p
    real(real64), intent(in) :: s(6)
    real(real64) :: e(6),mean
    mean=sum(s(1:3))/3
    e(1:3)=(s(1:3)-mean)/(2*p%shear)+mean/(3*p%bulk)
    e(4:6)=s(4:6)/p%shear
  end function
  pure subroutine invariants(s,mean,rho,theta,dev)
    real(real64), intent(in) :: s(6)
    real(real64), intent(out) :: mean,rho,theta,dev(6)
    real(real64) :: a(6),det,c,m(3,3),r(3,3),ev(3),phi,cs,sn,tmp
    integer :: sweep,i,j,ip,iq
    mean=sum(s(1:3))/3
    dev=s
    dev(1:3)=s(1:3)-mean
    rho=sqrt(sum(dev(1:3)**2)+2*sum(dev(4:6)**2))
    theta=0
    if (rho<=tiny(1.0_real64)) return
    a=dev/rho
    det=a(1)*a(2)*a(3)+2*a(4)*a(5)*a(6)-a(1)*a(6)**2-a(2)*a(5)**2-a(3)*a(4)**2
    c=max(-1.0_real64,min(1.0_real64,3*sqrt(6.0_real64)*det))
    theta=acos(c)/3
    ! Independent candidate: same Lode angle from ordered deviatoric eigenvalues.
    ! Avoid acos(cos(3 theta)) endpoint amplification near repeated eigenvalues.
    ! Retain original formula only if the bounded Jacobi iteration does not finish.
    m(1,:)=[a(1),a(4),a(5)];m(2,:)=[a(4),a(2),a(6)];m(3,:)=[a(5),a(6),a(3)]
    do sweep=1,60
      ip=1;iq=2
      do i=1,2;do j=i+1,3
        if(abs(m(i,j))>abs(m(ip,iq)))then;ip=i;iq=j;end if
      end do;end do
      if(abs(m(ip,iq))<8*epsilon(1._real64))exit
      phi=.5_real64*atan2(2*m(ip,iq),m(iq,iq)-m(ip,ip))
      cs=cos(phi);sn=sin(phi);r=0
      do i=1,3;r(i,i)=1;end do
      r(ip,ip)=cs;r(iq,iq)=cs;r(ip,iq)=sn;r(iq,ip)=-sn
      m=matmul(transpose(r),matmul(m,r))
    end do
    if(sweep<=60)then
      ev=[m(1,1),m(2,2),m(3,3)]
      do i=1,2;do j=i+1,3
        if(ev(j)>ev(i))then;tmp=ev(i);ev(i)=ev(j);ev(j)=tmp;end if
      end do;end do
      theta=min(acos(-1._real64)/3,atan2(sqrt(3._real64)*(ev(2)-ev(3)),2*ev(1)-ev(2)-ev(3)))
    end if
  end subroutine
  subroutine plastic_update(p,e,epold,kold,s,ep,knew,rc,branch,committed_effective)
    type(material_parameters), intent(in) :: p
    real(real64), intent(in) :: e(6),epold(6),kold
    real(real64), intent(inout) :: s(6),ep(6),knew
    integer, intent(out) :: rc,branch
    real(real64),optional,intent(in)::committed_effective(6)
    real(real64)::seed(4),om,orr,ot,od(6)
    real(real64)::alt(4),xs(4)
    integer::ar,ai
    real(real64) :: trial(6),dev(6),mean,rho,theta,x(4),v(3),sn(6),epn(6),kn,f
    integer :: it
    ! branch0 elastic,1 regular,2 vertex; -1 means not committed.
    rc=2
    branch=-1
    if (.not.all(ieee_is_finite([e,epold,kold]))) return
    ! Reuse regular entry's parameter validation using an elastic origin probe.
    x=0
    call regular_return(p,[0.0_real64,0.0_real64],0.0_real64,kold,x,rc,it)
    if (rc/=0) return
    trial=elastic_stress(p,e-epold)
    call invariants(trial,mean,rho,theta,dev)
    f=yield_value(mean,rho,theta,kold,p%fc,p%ft,p%ecc,p%q0,p%hp)
    kn=kold
    sn=trial
    if (f<=1e-10_real64) then
      branch=0
    else
      call vertex_return(p,[mean,rho],kold,v,rc,it)
      if (rc==0) then
        sn=0
        sn(1:3)=v(1)
        kn=v(3)
        branch=2
      else if (rc==1) then
        !1131: preserve all existing successful returns. One extra initial
        !guess only after existing bounded Newton attempts fail.
        call regular_return(p,[mean,rho],theta,kold,x,rc,it)
        if(present(committed_effective))then
          call invariants(committed_effective,om,orr,ot,od)
          seed=[om,orr,kold,0._real64]
          if(rc==1.or.rc==3.or.rc==4)then
            call regular_return_seed(p,[mean,rho],theta,kold,x,rc,it,60,seed)
          else if(rc==0)then
            !1134: one independent committed-state seed, not root enumeration.
            ! Two converged but distinct endpoint solutions are ambiguous:
            ! return failure for external cutback; never select by small kappa.
            alt=0
            call regular_return_seed(p,[mean,rho],theta,kold,alt,ar,ai,60,seed)
            if(ar==0)then
              xs=[p%fc,p%fc,max(1._real64,abs(x(3)),abs(alt(3))), &
                max(p%fc**2/p%bulk,abs(x(4)),abs(alt(4)))]
              if(maxval(abs(x-alt)/xs)>1e-3_real64)rc=1
            endif
          endif
        endif
        if (rc/=0) return
        if (rho<=tiny(1.0_real64)) then
          rc=1
          return
        end if
        ! Fixed Lode angle/directions in OOFEM regular return: radial dev scaling
        ! is equivalent to spectral reconstruction and needs no eigenvectors.
        sn=dev*(x(2)/rho)
        sn(1:3)=sn(1:3)+x(1)
        kn=x(3)
        branch=1
      else
        return
      end if
    end if
    epn=e-elastic_strain(p,sn)
    if (branch==0) epn=epold
    if (.not.all(ieee_is_finite([sn,epn,kn]))) then
      rc=2
      branch=-1
      return
    end if
    s=sn
    ep=epn
    knew=kn
    rc=0
  end subroutine
end module

module cdpm2_damage
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: equivalent_strain, compression_damage, tension_damage
contains
  subroutine tension_damage(eq,k1,k2,le,old,young,ft,ft1,wf1,wf,tol,omega,rc)
    ! Rate-independent bilinear law, with irreversibility applied on EVERY exit.
    ! Analytical intersection of (1-d)Eeq with traction(le*(k1+d*k2)).
    real(real64), intent(in) :: eq,k1,k2,le,old,young,ft,ft1,wf1,wf,tol
    real(real64), intent(inout) :: omega
    integer, intent(out) :: rc
    real(real64) :: intercept,slope,den,d,w,wtol,ftemp
    integer :: branch
    rc=2
    if (.not.all(ieee_is_finite([eq,k1,k2,le,old,young,ft,ft1,wf1,wf,tol]))) return
    if (min(k1,old)<0.or.old>1.or.min(le,young,ft,ft1,wf1)<=0.or.wf<=wf1) return
    if (tol<0.or.tol>=1) return
    ftemp=ft*(1-tol)
    if (ft1>=ftemp) return
    ! Same floating-point activation expression as the history and length rules.
    if (eq<=(ft/young)*(1-tol)) then
      ! OOFEM histories can be signed below initiation during unloading.
      ! No traction intersection is evaluated here; retain history, not clamp it.
      omega=old
      rc=0
      return
    end if
    if (k2<0) return
    wtol=64*epsilon(1.0_real64)*wf
    do branch=1,2
      if (branch==1) then
        intercept=ftemp
        slope=(ft1-ftemp)/wf1
      else
        slope=-ft1/(wf-wf1)
        intercept=ft1-slope*wf1
      end if
      den=young*eq+slope*le*k2
      if (abs(den)<=epsilon(1.0_real64)*max(young*eq,abs(slope*le*k2))) cycle
      d=(young*eq-intercept-slope*le*k1)/den
      w=le*(k1+d*k2)
      if ((branch==1.and.w>=-wtol.and.w<=wf1+wtol).or. &
          (branch==2.and.w>=wf1-wtol.and.w<=wf+wtol)) then
        omega=max(old,min(1.0_real64,d))
        rc=0
        return
      end if
    end do
    if (le*(k1+k2)>=wf-wtol) then
      omega=1
      rc=0
      return
    end if
    rc=1
  end subroutine
  pure real(real64) function equivalent_strain(sig,rho,theta,fc,ft,ecc,young) result(eq)
    real(real64), intent(in) :: sig,rho,theta,fc,ft,ecc,young
    real(real64) :: r,m,ph,qh,c
    c=cos(theta)
    r=(4*(1-ecc**2)*c**2+(2*ecc-1)**2)/ &
      (2*(1-ecc**2)*c+(2*ecc-1)*sqrt(4*(1-ecc**2)*c**2+5*ecc**2-4*ecc))
    m=3*(fc**2-ft**2)/(fc*ft)*ecc/(ecc+1)
    ph=-m*(rho*r/(sqrt(6.0_real64)*fc)+sig/fc)
    qh=-1.5_real64*(rho/fc)**2
    eq=max(0.0_real64,-ph/2+sqrt(ph**2/4-qh))*ft/young
  end function
  subroutine compression_damage(eq,k1,k2,old,young,ft,ef,tol,omega,rc)
    real(real64), intent(in) :: eq,k1,k2,old,young,ft,ef,tol
    real(real64), intent(inout) :: omega
    integer, intent(out) :: rc
    real(real64) :: lo,hi,mid,f,ftemp
    integer :: i
    rc=2
    if (.not.all(ieee_is_finite([eq,k1,k2,old,young,ft,ef,tol]))) return
    if (min(k1,old)<0.or.old>1.or.min(young,ft,ef)<=0.or.tol<0.or.tol>=1) return
    ftemp=ft*(1-tol)
    ! Same floating-point activation expression as the history and length rules.
    if (eq<=(ft/young)*(1-tol)) then
      ! Same inactive-domain rule as OOFEM; signed auxiliary histories are unused.
      omega=old
      rc=0
      return
    end if
    lo=old
    hi=1
    ! OOFEM's signed threshold increments can leave kappaTwo negative even
    ! above initiation. For k2<0, F'(omega)=-E*eq+(k2/ef)*traction <0:
    ! the existing [old,1] bracket therefore has at most one root.
    ! Reject unrepresentable exponential trials, NOT signed histories per se.
    if (-(k1+min(0.0_real64,k2))/ef>log(huge(1.0_real64))/4) return
    f=(1-lo)*young*eq-ftemp*exp(-(k1+lo*k2)/ef)
    if (f<=0) then
      omega=old
      rc=0
      return
    end if
    ! 760: residual-only bisection gives a quantized damage map, contaminating
    ! strain derivatives even in real128. Also resolve the root bracket to
    ! representational precision. The residual criterion is NOT relaxed.
    ! real64 keeps the original100-iteration bound; real128 diagnostics need
    ! additional bits. Constitutive/activation parameters remain unchanged.
    do i=1,max(100,digits(1.0_real64)+8)
      mid=lo+(hi-lo)/2
      f=(1-mid)*young*eq-ftemp*exp(-(k1+mid*k2)/ef)
      if (f>0) then
        lo=mid
      else
        hi=mid
      end if
      if (hi-lo<=16*epsilon(mid)*max(1.0_real64,abs(mid))) then
        ! Upper damage bracket: admissible side, never above the physical1.
        mid=hi
        f=(1-mid)*young*eq-ftemp*exp(-(k1+mid*k2)/ef)
        if (abs(f)/ft<1e-10_real64) then
          omega=mid
          rc=0
          return
        end if
      end if
    end do
    rc=1
  end subroutine
end module

module cdpm2_split
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: split_stress, damaged_stress
contains
  subroutine split_stress(s,plus,minus,alpha,rc)
    real(real64), intent(in) :: s(6)
    real(real64), intent(inout) :: plus(6),minus(6),alpha
    integer, intent(out) :: rc
    real(real64) :: a(3,3),v(3,3),r(3,3),ap(3,3),am(3,3),scale,phi,c,t,ev(3),den
    integer :: i,j,k,p,q,sweep
    rc=2
    if (.not.all(ieee_is_finite(s))) return
    scale=maxval(abs(s))
    if (scale<=tiny(1.0_real64)) then
      plus=0
      minus=0
      alpha=1 ! OOFEM zero-stress convention.
      rc=0
      return
    end if
    a(1,:)=[s(1),s(4),s(5)]/scale
    a(2,:)=[s(4),s(2),s(6)]/scale
    a(3,:)=[s(5),s(6),s(3)]/scale
    v=0
    do i=1,3
      v(i,i)=1
    end do
    rc=1
    do sweep=1,60
      p=1
      q=2
      do i=1,2
        do j=i+1,3
          if (abs(a(i,j))>abs(a(p,q))) then
            p=i
            q=j
          end if
        end do
      end do
      if (abs(a(p,q))<8*epsilon(1.0_real64)) exit
      phi=.5_real64*atan2(2*a(p,q),a(q,q)-a(p,p))
      c=cos(phi)
      t=sin(phi)
      r=0
      do i=1,3
        r(i,i)=1
      end do
      r(p,p)=c
      r(q,q)=c
      r(p,q)=t
      r(q,p)=-t
      a=matmul(transpose(r),matmul(a,r))
      v=matmul(v,r)
    end do
    if (sweep>60) return
    ev=[a(1,1),a(2,2),a(3,3)]
    ap=0
    am=0
    do k=1,3
      do i=1,3
        do j=1,3
          ap(i,j)=ap(i,j)+max(ev(k),0.0_real64)*v(i,k)*v(j,k)
          am(i,j)=am(i,j)+min(ev(k),0.0_real64)*v(i,k)*v(j,k)
        end do
      end do
    end do
    plus=[ap(1,1),ap(2,2),ap(3,3),ap(1,2),ap(1,3),ap(2,3)]*scale
    minus=[am(1,1),am(2,2),am(3,3),am(1,2),am(1,3),am(2,3)]*scale
    den=sum(ev**2)
    alpha=max(0.0_real64,min(1.0_real64,sum(min(ev,0.0_real64)**2)/den))
    rc=0
  end subroutine
  subroutine damaged_stress(effective,dt,dc,nominal,rc)
    real(real64), intent(in) :: effective(6),dt,dc
    real(real64), intent(inout) :: nominal(6)
    integer, intent(out) :: rc
    real(real64) :: sp(6),sm(6),alpha
    rc=2
    if (.not.all(ieee_is_finite([dt,dc]))) return
    if (min(dt,dc)<0.or.max(dt,dc)>1) return
    call split_stress(effective,sp,sm,alpha,rc)
    if (rc/=0) return
    nominal=(1-dt)*sp+(1-dc)*sm
  end subroutine
end module

module cdpm2_damage_history
  use geometry_length_143
  use length_state_145
  use strain_direction_146
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use cdpm2_return
  use cdpm2_functions, only: hardening_two
  use cdpm2_plastic, only: invariants
  use cdpm2_damage
  use cdpm2_split
  implicit none
  private
  public :: damage_history, damage_parameters, update_damage
  type damage_history
    real(real64) :: selected_le=0
    real(real64) :: equivalent=0, accumulated(2)=0, alpha=0
    real(real64) :: threshold(2)=0, one(2)=0, two(2)=0, damage(2)=0
  end type
  type damage_parameters
    logical :: automatic_length=.false.
    real(real64) :: reference_nodes(3,8)=0,mean_size=0
    real(real64) :: young=35892.339_real64,le=20.0_real64
    real(real64) :: ft1=.814932615_real64,wf1=.021974493916716565_real64
    real(real64) :: wf=.14649662611144376_real64,efc=1e-4_real64
    real(real64) :: asoft=15.0_real64,bsoft=1.0_real64
    ! Pinned OOFEM initializeFrom: yieldTolDamage=10*yieldTol, default yieldTol=1e-6.
    ! Independent of e0=ft/E and of our local plastic Newton stopping tolerance.
    real(real64) :: damage_tolerance=1e-5_real64
  end type
contains
  subroutine update_damage(p,d,old,sold,snew,dep,kp,new,nominal,rc,diagnostic,total_strain)
    type(material_parameters), intent(in) :: p
    type(damage_parameters), intent(in) :: d
    type(damage_history), intent(in) :: old
    type(damage_history), intent(inout) :: new
    real(real64), intent(in) :: sold(6),snew(6),dep(6),kp
    real(real64), intent(inout) :: nominal(6)
    integer, intent(out) :: rc
    real(real64), optional, intent(out) :: diagnostic(5)
    real(real64), optional, intent(in) :: total_strain(6)
    real(real64) :: normal(3),span,lambda,gap,used_le
    type(damage_history) :: trial
    real(real64) :: sp(6),sm(6),alpha,eq,emin,eqplus,eqminus,mid,sig,rho,theta,dev(6)
    real(real64) :: e0,tol,rs,duct,pnorm,factor,extra,delta,omega,sn(6)
    integer :: i,j
    logical :: reversal
    rc=2
    if (present(diagnostic)) diagnostic=0
    if (d%young<=0.or.(.not.d%automatic_length.and.d%le<=0).or.d%asoft<1.or.d%bsoft<=0) return
    if (.not.ieee_is_finite(d%damage_tolerance)) return
    if (d%damage_tolerance<0.or.d%damage_tolerance>=1) return
    e0=p%ft/d%young
    tol=d%damage_tolerance
    call split_stress(snew,sp,sm,alpha,rc)
    if (rc/=0) return
    eq=equiv(snew)
    eqplus=equiv(sold+.01_real64*(snew-sold))
    eqminus=equiv(sold+.99_real64*(snew-sold))
    emin=old%equivalent
    reversal=old%equivalent>eqplus.and.eq>eqminus.and. &
      abs(eqplus-old%equivalent)>tol/100.and.abs(eqminus-eq)>tol/100
    if (reversal) then
      do i=1,100
        mid=equiv(sold+real(i,real64)/100*(snew-sold))
        if (mid>emin) exit
        emin=mid
      end do
    end if
    trial=old
    trial%equivalent=eq
    trial%alpha=alpha
    trial%accumulated(1)=old%accumulated(1)+eq-old%equivalent
    if (reversal) then
      trial%accumulated(2)=old%accumulated(2)+old%alpha*(emin-old%equivalent)+alpha*(eq-emin)
    else
      trial%accumulated(2)=old%accumulated(2)+alpha*(eq-old%equivalent)
    end if
    call invariants(snew,sig,rho,theta,dev)
    rs=0
    if (sig<0) rs=-sig/(.40824829_real64*max(rho,1e-16_real64))
    duct=1+(d%asoft-1)*rs**d%bsoft
    do j=1,2
      delta=trial%accumulated(j)-old%threshold(j)
      if (delta/e0 < -tol) cycle
      trial%threshold(j)=trial%accumulated(j)
      ! Correct engineering shear norm (OOFEM raw Euclidean norm differs).
      pnorm=norm2([dep(1:3),dep(4:6)/sqrt(2.0_real64)])
      if (j==2) pnorm=alpha*pnorm
      ! At activation itself there is no post-activation plastic portion.
      if (trial%threshold(j)<=e0*(1-tol)) then
        pnorm=0
      else if (trial%threshold(j)>e0*(1-tol).and.old%threshold(j)<e0*(1-tol)) then
        ! DIAGNOSTIC ONLY: align interpolation with the existing activation level.
        ! The OOFEM source uses e0 here but e0*(1-tol) for activation, which can
        ! give a negative fraction and a positive norm jump on first activation.
        factor=1-(e0*(1-tol)-old%threshold(j))/(trial%threshold(j)-old%threshold(j))
        pnorm=abs(factor)*pnorm ! Retains original norm of the scaled increment.
      end if
      if (j==2) then
        extra=p%ft*hardening_two(kp,p%hp)*sqrt(2.0_real64/3)/ &
          (max(rho,1e-16_real64)*sqrt(1+2*p%dilation**2))
        pnorm=pnorm*extra
      end if
      trial%one(j)=old%one(j)+pnorm/duct
      trial%two(j)=old%two(j)+delta/duct
      omega=old%damage(j)
      if (present(diagnostic)) diagnostic=[real(j,real64),trial%threshold(j),trial%one(j),trial%two(j),omega]
      if (j==1) then
        used_le=d%le
        if(d%automatic_length) then
          ! Preserve positive placeholder below activation: tension_damage validates le.
          used_le=1.0_real64
          if(trial%threshold(1)>e0*(1-tol)) then
            span=0
            if(all(old%damage==0)) then
              if(.not.present(total_strain)) then
                rc=2;return
              end if
              call maximum_strain_direction(total_strain,normal,lambda,gap,rc,d%reference_nodes)
              if(rc/=0) then
                rc=2;return ! Explicit unsupported degeneracy, not cutback tuning.
              end if
              call projected_span(d%reference_nodes,normal,span,rc)
              if(rc/=0) then
                rc=2;return
              end if
            end if
            call select_length(old%selected_le,old%damage,.true.,trial%threshold(1), &
              e0,tol,0.0_real64,span,d%mean_size,trial%selected_le,rc)
            if(rc/=0) then
              rc=2;return
            end if
            used_le=trial%selected_le
          end if
        end if
        call tension_damage(trial%threshold(j),trial%one(j),trial%two(j),used_le,omega, &
          d%young,p%ft,d%ft1,d%wf1,d%wf,tol,trial%damage(j),rc)
      else
        call compression_damage(trial%threshold(j),trial%one(j),trial%two(j),omega, &
          d%young,p%ft,d%efc,tol,trial%damage(j),rc)
      end if
      if (rc/=0) return
    end do
    call damaged_stress(snew,trial%damage(1),trial%damage(2),sn,rc)
    if (rc/=0) return
    nominal=sn
    new=trial
  contains
    real(real64) function equiv(s) result(value)
      real(real64), intent(in) :: s(6)
      real(real64) :: mean,r,t,v(6)
      call invariants(s,mean,r,t,v)
      value=equivalent_strain(mean,r,t,p%fc,p%ft,p%ecc,d%young)
    end function
  end subroutine
end module

module cdpm2_model
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  use cdpm2_damage_history,only:damage_history
  implicit none
  private
  public :: material_state
  type material_state
    real(real64) :: strain(6)=0,plastic(6)=0,kappa=0,effective(6)=0,stress(6)=0
    type(damage_history) :: history
  end type
end module

module state_geometry_150
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 use cdpm2_model,only:material_state
 implicit none
 integer,parameter::geometry_state_size=39
contains
 subroutine unpack_geometry(v,s,rc)
  real(real64),intent(in)::v(:)
  type(material_state),intent(inout)::s
  integer,intent(out)::rc
  type(material_state)::t
  rc=2
  if(size(v)<39)return
  if(.not.all(ieee_is_finite(v(1:39))))return
  if(all(v(1:39)==0))then
   s=t;rc=0;return
  end if
  if(v(1)/=2)return ! Never silently migrate version1.
  if(v(14)<0.or.v(27)<0.or.v(28)<0.or.v(28)>1)return
  if(any(v(33:34)<0).or.any(v(37:38)<0).or.any(v(37:38)>1).or.v(39)<0)return
  t%strain=v(2:7);t%plastic=v(8:13);t%kappa=v(14)
  t%effective=v(15:20);t%stress=v(21:26)
  t%history%equivalent=v(27);t%history%alpha=v(28)
  t%history%accumulated=v(29:30);t%history%threshold=v(31:32)
  t%history%one=v(33:34);t%history%two=v(35:36);t%history%damage=v(37:38)
  t%history%selected_le=v(39)
  s=t;rc=0
 end subroutine
 subroutine pack_geometry(s,v,rc)
  type(material_state),intent(in)::s
  real(real64),intent(inout)::v(:)
  integer,intent(out)::rc
  real(real64)::a(39)
  type(material_state)::checked
  rc=2;if(size(v)<39)return
  a=[2.0_real64,s%strain,s%plastic,s%kappa,s%effective,s%stress, &
     s%history%equivalent,s%history%alpha,s%history%accumulated, &
     s%history%threshold,s%history%one,s%history%two,s%history%damage,s%history%selected_le]
  call unpack_geometry(a,checked,rc)
  if(rc/=0)return
  v(1:39)=a
 end subroutine
end module

module geometry_payload_159
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine decode_geometry(payload,nlgeom,frame,nodes,mean_size,rc)
  real(real64),intent(in)::payload(:)
  integer,intent(in)::nlgeom,frame ! only0: small displacement/global frame
  real(real64),intent(inout)::nodes(3,8),mean_size
  integer,intent(out)::rc
  real(real64)::x(3,8),expected(3,8),edge,scale,tol
  integer::i
  rc=2
  if(size(payload)/=27.or.nlgeom/=0.or.frame/=0)return
  if(.not.all(ieee_is_finite(payload)))return
  if(payload(1)/=1.or.payload(2)/=1)return
  edge=payload(3);if(edge<=0)return
  x=reshape(payload(4:27),[3,8])
  expected=reshape([-1, -1, -1,  1, -1, -1,  1, 1, -1, -1, 1, -1, &
                    -1, -1,  1,  1, -1,  1,  1, 1,  1, -1, 1,  1],[3,8])*edge/2
  scale=max(edge,maxval(abs(x)))
  tol=64*epsilon(edge)*scale ! geometry roundoff only, no material tolerance
  if(maxval(abs(x-expected))>tol)return
  ! Avoid degenerate floating-point edge coordinates at extreme magnitudes.
  do i=1,3
   if(maxval(x(i,:))-minval(x(i,:))<=0)return
  end do
  nodes=x;mean_size=edge;rc=0
 end subroutine
end module

module properties_161
 use iso_fortran_env,only:real64
 use cdpm2_return,only:material_parameters
 use cdpm2_damage_history,only:damage_parameters
 use cdpm2_return,only:regular_return
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
  subroutine decode_properties(props,p,d,cap,h,budget,rc)
    real(real64), intent(in) :: props(:)
    type(material_parameters), intent(out) :: p
    type(damage_parameters), intent(out) :: d
    real(real64), intent(out) :: cap,h
    integer, intent(out) :: budget,rc
    real(real64) :: state(4),zero(2),nu
    integer :: iterations
    rc=2
    if (size(props)/=23) return
    if (.not.all(ieee_is_finite(props))) return
    if (props(1)<=0.or.props(2)<=-1.or.props(2)>=.5_real64) return
    if (minval(props(13:19))<=0.or.props(18)<1) return
    if (props(16)<=props(15).or.props(20)<0.or.props(20)>=1) return
    if (props(14)>=props(4)*(1-props(20))) return
    if (minval(props(21:22))<=0.or.props(23)<1.or.props(23)>1000000) return
    if (props(23)/=aint(props(23))) return
    nu=props(2)
    p%bulk=props(1)/(3*(1-2*nu)); p%shear=props(1)/(2*(1+nu))
    p%fc=props(3); p%ft=props(4); p%ecc=props(5); p%q0=props(6)
    p%hp=props(7); p%dilation=props(8)
    p%ah=props(9); p%bh=props(10); p%ch=props(11); p%dh=props(12)
    d%young=props(1); d%le=props(13); d%ft1=props(14)
    d%wf1=props(15); d%wf=props(16); d%efc=props(17)
    d%asoft=props(18); d%bsoft=props(19); d%damage_tolerance=props(20)
    cap=props(21); h=props(22); budget=int(props(23))
    zero=0; state=0
    call regular_return(p,zero,0.0_real64,0.0_real64,state,rc,iterations)
  end subroutine

end module

module rotation_281
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine rotate_state(q,v,out,rc)
 real(real64),intent(in)::q(3,3),v(39)
 real(real64),intent(inout)::out(39)
 integer,intent(out)::rc
 real(real64)::eye(3,3),det,t(39)
 integer::i
 rc=2
 if(.not.all(ieee_is_finite(q)).or..not.all(ieee_is_finite(v)))return
 if(v(1)/=2.and..not.all(v==0))return
 eye=0;do i=1,3;eye(i,i)=1;end do
 if(maxval(abs(matmul(transpose(q),q)-eye))>1e-12_real64)return
 det=q(1,1)*(q(2,2)*q(3,3)-q(2,3)*q(3,2))-q(1,2)*(q(2,1)*q(3,3)-q(2,3)*q(3,1)) &
     +q(1,3)*(q(2,1)*q(3,2)-q(2,2)*q(3,1))
 if(abs(det-1)>1e-12_real64)return
 t=v
 call tensor(q,v(2:7),2._real64,t(2:7))
 call tensor(q,v(8:13),2._real64,t(8:13))
 call tensor(q,v(15:20),1._real64,t(15:20))
 call tensor(q,v(21:26),1._real64,t(21:26))
 if(.not.all(ieee_is_finite(t)))return
 out=t;rc=0
 end subroutine
 subroutine tensor(q,v,f,out)
 real(real64),intent(in)::q(3,3),v(6),f
 real(real64),intent(out)::out(6)
 real(real64)::a(3,3),b(3,3)
 a(1,:)=[v(1),v(4)/f,v(5)/f]
 a(2,:)=[v(4)/f,v(2),v(6)/f]
 a(3,:)=[v(5)/f,v(6)/f,v(3)]
 b=matmul(matmul(q,a),transpose(q))
 out=[b(1,1),b(2,2),b(3,3),f*b(1,2),f*b(1,3),f*b(2,3)]
 end subroutine
end module

module rotation_mock_282
 use rotation_281,only:rotate_state,real64
 use state_geometry_150,only:unpack_geometry,material_state
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine align_state(q,old,external_s,external_e,out,rc)
 real(real64),intent(in)::q(3,3),old(39),external_s(6),external_e(6)
 real(real64),intent(inout)::out(39)
 integer,intent(out)::rc
 real(real64)::trial(39)
 type(material_state)::checked
 rc=2
 if(.not.all(ieee_is_finite(external_s)).or..not.all(ieee_is_finite(external_e)))return
 call unpack_geometry(old,checked,rc)
 if(rc/=0)return
 call rotate_state(q,old,trial,rc)
 if(rc/=0)return
 rc=2
 if(maxval(abs(external_s-trial(21:26)))>1e-10_real64*max(1._real64,maxval(abs(external_s))))return
 if(maxval(abs(external_e-trial(2:7)))>1e-10_real64*max(1._real64,maxval(abs(external_e))))return
 call unpack_geometry(trial,checked,rc)
 if(rc/=0)return
 out=trial
 end subroutine
end module

module frame_transport_378
 use rotation_281,only:rotate_state,real64
 use rotation_mock_282,only:align_state
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine stage_frame(drot,qold,x0,v,s,e,qnew,xnew,vnew,rc)
 real(real64),intent(in)::drot(3,3),qold(3,3),x0(3,8),v(39),s(6),e(6)
 real(real64),intent(inout)::qnew(3,3),xnew(3,8),vnew(39)
 integer,intent(out)::rc
 real(real64)::zero(39),check(39),qt(3,3),xt(3,8),vt(39)
 rc=2
 if(.not.all(ieee_is_finite(x0)))return
 zero=0;check=0
 call rotate_state(qold,zero,check,rc)
 if(rc/=0)return
 call align_state(drot,v,s,e,vt,rc)
 if(rc/=0)return
 qt=matmul(drot,qold)
 call rotate_state(qt,zero,check,rc)
 if(rc/=0)return
 xt=matmul(qt,x0)
 rc=2
 if(.not.all(ieee_is_finite(xt)))return
 qnew=qt;xnew=xt;vnew=vt;rc=0
 end subroutine
end module

module volume_map_376
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine volume_jacobian(s,c,out,rc)
 real(real64),intent(in)::s(6),c(6,6)
 real(real64),intent(inout)::out(6,6)
 integer,intent(out)::rc
 real(real64)::t(6,6)
 integer::j
 rc=2
 if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(c)))return
 t=c
 do j=1,3
  t(:,j)=t(:,j)+s
 end do
 if(.not.all(ieee_is_finite(t)))return
 out=t;rc=0
 end subroutine
end module

module geometry_general_404
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 pure real(real64) function determinant(a) result(d)
 real(real64),intent(in)::a(3,3)
 d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
  +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
 end function
 pure subroutine derivatives(nn,r,s,t,dn)
 integer,intent(in)::nn
 real(real64),intent(in)::r,s,t
 real(real64),intent(out)::dn(8,3)
 integer,parameter::sx(8)=[-1,1,1,-1,-1,1,1,-1],sy(8)=[-1,-1,1,1,-1,-1,1,1],sz(8)=[-1,-1,-1,-1,1,1,1,1]
 real(real64)::b(3),br(3),bs(3),f
 integer::i,j,k
 dn=0
 if(nn==8)then
  do i=1,8
   dn(i,1)=sx(i)*(1+sy(i)*s)*(1+sz(i)*t)/8._real64
   dn(i,2)=sy(i)*(1+sx(i)*r)*(1+sz(i)*t)/8._real64
   dn(i,3)=sz(i)*(1+sx(i)*r)*(1+sy(i)*s)/8._real64
  end do
 else
  b=[1-r-s,r,s];br=[-1._real64,1._real64,0._real64];bs=[-1._real64,0._real64,1._real64]
  do k=0,1
   f=2*k-1
   do j=1,3
    i=j+3*k;dn(i,1)=br(j)*(1+f*t)/2;dn(i,2)=bs(j)*(1+f*t)/2;dn(i,3)=b(j)*f/2
   end do
  end do
 end if
 end subroutine
 subroutine prepare_geometry(x,nodes,volume,mean_size,min_j,rc)
 real(real64),intent(in)::x(:,:)
 real(real64),intent(inout)::nodes(3,8),volume,mean_size,min_j
 integer,intent(out)::rc
 real(real64)::xc(3,8),center(3),dn(8,3),jac,vol,mj,rr,ss,tt,g,w,l
 real(real64),parameter::tri(2,3)=reshape([1._real64/6,1._real64/6,2._real64/3,1._real64/6,1._real64/6,2._real64/3],[2,3])
 real(real64),parameter::corn(2,3)=reshape([0._real64,0._real64,1._real64,0._real64,0._real64,1._real64],[2,3])
 integer::nn,i,j,k,pass
 rc=2;nn=size(x,2)
 if(size(x,1)/=3.or.(nn/=6.and.nn/=8))return
 if(.not.all(ieee_is_finite(x)))return
 center=sum(x,dim=2)/nn;xc=0
 do i=1,nn
  xc(:,i)=x(:,i)-center
 end do
 vol=0;mj=huge(mj);g=1/sqrt(3._real64)
 do pass=1,2
  do k=-1,1,2
   tt=k;if(pass==1)tt=tt*g
   if(nn==8)then
    do j=-1,1,2
     ss=j;if(pass==1)ss=ss*g
     do i=-1,1,2
      rr=i;if(pass==1)rr=rr*g
      call derivatives(nn,rr,ss,tt,dn);jac=determinant(matmul(xc,dn))
      if(.not.ieee_is_finite(jac).or.jac<=0)return
      mj=min(mj,jac);if(pass==1)vol=vol+jac
     end do
    end do
   else
    do i=1,3
     rr=corn(1,i);ss=corn(2,i);w=0
     if(pass==1)then
      rr=tri(1,i);ss=tri(2,i);w=1._real64/6
     end if
     call derivatives(nn,rr,ss,tt,dn);jac=determinant(matmul(xc,dn))
     if(.not.ieee_is_finite(jac).or.jac<=0)return
     mj=min(mj,jac);vol=vol+w*jac
    end do
   end if
  end do
 end do
 if(.not.ieee_is_finite(vol).or.vol<=0)return
 l=vol**(1._real64/3)
 if(.not.ieee_is_finite(l).or.l<=0)return
 ! Repeated nodes preserve projection extrema exactly; not wedge shape functions.
 if(nn==6)then
  xc(:,7)=xc(:,6);xc(:,8)=xc(:,6)
 end if
 nodes=xc;volume=vol;mean_size=l;min_j=mj;rc=0
 end subroutine
end module

module geometry_payload_406
 use geometry_payload_159,only:decode_geometry
 use geometry_general_404,only:prepare_geometry,real64,ieee_is_finite
 implicit none
contains
 subroutine decode_initial_geometry(payload,nodes,mean_size,rc)
 real(real64),intent(in)::payload(:)
 real(real64),intent(inout)::nodes(3,8),mean_size
 integer,intent(out)::rc
 real(real64)::x(3,8),xt(3,8),v,l,j
 integer::nn
 rc=2
 if(size(payload)/=27)return
 if(.not.all(ieee_is_finite(payload)))return
 if(payload(1)==1)then
  call decode_geometry(payload,0,0,nodes,mean_size,rc);return
 end if
 if(payload(1)/=2.or.payload(3)/=0)return
 if(payload(2)/=6.and.payload(2)/=8)return
 nn=int(payload(2));x=reshape(payload(4:27),[3,8])
 if(nn==6.and.any(x(:,7:8)/=0))return
 xt=0;v=0;l=0;j=0
 call prepare_geometry(x(:,:nn),xt,v,l,j,rc)
 if(rc/=0)return
 nodes=xt;mean_size=l
 end subroutine
end module

module corot_core_406
 use geometry_payload_406,only:decode_initial_geometry
 use properties_161
 use frame_transport_378,only:stage_frame
 use volume_map_376,only:volume_jacobian
 implicit none
end module

module geometry_provider_421
 use iso_fortran_env,only:real64
 use geometry_payload_406,only:decode_initial_geometry
 implicit none
 private
 public::load_geometry_table,get_element_geometry,clear_geometry_table,geometry_count,install_geometry
 integer,allocatable,save::labels(:)
 real(real64),allocatable,save::geometry(:,:)
contains
 subroutine clear_geometry_table()
 if(allocated(labels))deallocate(labels)
 if(allocated(geometry))deallocate(geometry)
 end subroutine
 integer function geometry_count()
 geometry_count=0
 if(allocated(labels))geometry_count=size(labels)
 end function
 subroutine load_geometry_table(path,rc)
 character(*),intent(in)::path
 integer,intent(out)::rc
 integer::u,ios,n,i,nn,stat,id,original,back_label,back_rc
 integer,allocatable::lt(:)
 real(real64),allocatable::gt(:,:)
 real(real64)::x(24),nodes(3,8),length
 character(80)::magic,instance,back_instance
 character(2048)::record
 character(80)::extra_token
 logical::instance_table
 call clear_geometry_table();rc=2
 open(newunit=u,file=path,status='old',action='read',iostat=ios)
 if(ios/=0)return
 read(u,'(a)',iostat=ios)magic
 if(ios/=0)goto 900
 instance_table=trim(magic)=='CDPM2_GEOMETRY_V2'
 if(.not.instance_table.and.trim(magic)/='CDPM2_GEOMETRY_V1')goto 900
 if(instance_table)then
  read(u,'(a)',iostat=ios)magic
  if(ios/=0)goto 900
  if(trim(magic)/='GLOBAL_REFERENCE_NO_ORIENTATION')goto 900
 endif
 read(u,*,iostat=ios)n
 if(ios/=0)goto 900
 if(n<1.or.n>1000000)goto 900
 allocate(lt(n),gt(27,n),stat=stat)
 if(stat/=0)goto 900
 do i=1,n
  if(instance_table)then
   read(u,'(a)',iostat=ios)record
   if(ios/=0)goto 900
   read(record,*,iostat=ios)instance,original,nn,x
   if(ios/=0.or.original<=0)goto 900
   read(record,*,iostat=ios)instance,original,nn,x,extra_token
   if(ios>=0)goto 900
   ! Empty quoted name denotes an assembly-level element (official convention).
   call GETINTERNAL(instance,original,1,id,stat)
   if(stat/=0.or.id<=0)goto 900
   call GETPARTINFO(id,1,back_instance,back_label,back_rc)
   if(back_rc/=0.or.back_label/=original)goto 900
   if(trim(back_instance)/=trim(instance))goto 900
  else
   ! V1 is the historical flat-model/internal-ID table. Do not reinterpret it.
   read(u,*,iostat=ios)id,nn,x
   if(ios/=0)goto 900
   if(id<=0)goto 900
   if(i>1)then
    if(id<=lt(i-1))goto 900
   endif
  end if
  gt(:,i)=[2._real64,real(nn,real64),0._real64,x]
  nodes=0;length=0
  call decode_initial_geometry(gt(:,i),nodes,length,stat)
  if(stat/=0)goto 900
  lt(i)=id
 end do
 ! Any nonblank trailing record is rejected; no silently ignored second table.
 do
  read(u,'(a)',iostat=ios)magic
  if(ios<0)exit
  if(ios/=0)goto 900
  if(len_trim(magic)>0)goto 900
 end do
 if(instance_table)then
  call sort_geometry(lt,gt,1,n)
  do i=2,n
   if(lt(i)<=lt(i-1))goto 900
  enddo
 endif
 call move_alloc(lt,labels);call move_alloc(gt,geometry);rc=0
900 close(u)
 end subroutine
 subroutine install_geometry(ids,data,rc)
 integer,intent(in)::ids(:)
 real(real64),intent(in)::data(:,:)
 integer,intent(out)::rc
 integer,allocatable::lt(:)
 real(real64),allocatable::gt(:,:)
 real(real64)::x(3,8),length
 integer::i,stat
 call clear_geometry_table();rc=2
 if(size(ids)<1.or.size(data,1)/=27.or.size(data,2)/=size(ids))return
 allocate(lt(size(ids)),gt(27,size(ids)),stat=stat);if(stat/=0)return
 lt=ids;gt=data
 do i=1,size(ids)
  if(lt(i)<=0)return
  x=0;length=0;call decode_initial_geometry(gt(:,i),x,length,stat)
  if(stat/=0)return
 enddo
 call sort_geometry(lt,gt,1,size(lt))
 do i=2,size(lt)
  if(lt(i)<=lt(i-1))return
 enddo
 call move_alloc(lt,labels);call move_alloc(gt,geometry);rc=0
 end subroutine
 recursive subroutine sort_geometry(ids,data,lo,hi)
 integer,intent(inout)::ids(:)
 real(real64),intent(inout)::data(:,:)
 integer,intent(in)::lo,hi
 integer::i,j,pivot,tmp
 real(real64)::row(27)
 if(lo>=hi)return
 i=lo;j=hi;pivot=ids(lo+(hi-lo)/2)
 do
  do while(ids(i)<pivot)
   i=i+1
  enddo
  do while(ids(j)>pivot)
   j=j-1
  enddo
  if(i<=j)then
   tmp=ids(i);ids(i)=ids(j);ids(j)=tmp
   row=data(:,i);data(:,i)=data(:,j);data(:,j)=row
   i=i+1;j=j-1
  endif
  if(i>j)exit
 enddo
 if(lo<j)call sort_geometry(ids,data,lo,j)
 if(i<hi)call sort_geometry(ids,data,i,hi)
 end subroutine
 subroutine get_element_geometry(noel,payload,rc)
 integer,intent(in)::noel
 real(real64),intent(inout)::payload(27)
 integer,intent(out)::rc
 integer::lo,hi,mid
 rc=2
 if(.not.allocated(labels))return
 lo=1;hi=size(labels)
 do while(lo<=hi)
  mid=lo+(hi-lo)/2
  if(labels(mid)==noel)then
   payload=geometry(:,mid);rc=0;return
  else if(labels(mid)<noel)then
   lo=mid+1
  else
   hi=mid-1
  end if
 end do
 end subroutine
end module

module tangent_observations
 use iso_fortran_env,only:int64
 implicit none
 private
 public::reset_tangent_observations,record_tangent_observations,report_tangent_observations
 integer(int64),save::calls=0,counts(8)=0
 character(40),parameter::names(8)=[character(40):: &
  'vertex_centered_fd','vertex_one_sided_fd','vertex_elastic_derivative', &
  'regular_elastic_derivative','damage_derivative_frozen','spectral_local_fd', &
  'zero_increment_frozen_query','selected_length_derivative_frozen']
contains
 subroutine reset_tangent_observations()
  calls=0;counts=0
 end subroutine
 subroutine record_tangent_observations(diag,noel,npt,kstep,kinc)
  ! Caller holds logger Mutex1. No extra material-call mutex is introduced.
  integer,intent(in)::diag(8),noel,npt,kstep,kinc
  integer::i,ios,rank
  calls=calls+1
  do i=1,8
   if(diag(i)<=0)cycle
   if(i>=3.and.i<=5.and.counts(i)==0)then
    call GETRANK(rank)
    write(6,*,iostat=ios)'CDPM2_TANGENT_WARNING trial-call approximation:', &
     trim(names(i)),' rank/NOEL/NPT/STEP/INC=',rank,noel,npt,kstep,kinc
    ! Diagnostic I/O failure must not change material acceptance.
   endif
   counts(i)=counts(i)+int(diag(i),int64)
  enddo
 end subroutine
 subroutine report_tangent_observations()
  integer::i,ios,rank
  call GETRANK(rank)
  write(6,*,iostat=ios)'CDPM2_TANGENT_SUMMARY rank/successful_trial_calls=',rank,calls
  write(6,*,iostat=ios)'Counts are session-local trial events; not converged increments or restart state.'
  do i=1,8
   write(6,*,iostat=ios)'CDPM2_TANGENT_EVENT ',trim(names(i)),counts(i)
  enddo
 end subroutine
end module

module trial_log_391
 use tangent_observations,only:record_tangent_observations
 use iso_fortran_env,only:real64
 implicit none
 private
 public::write_corot_trial
 logical,save::initialized=.false.,failures_only=.false.,mapping_done=.false.
 integer,save::limit=0,written=0,config_error=0,target_noel=0
 character(2048),save::path=''
contains
 subroutine write_corot_trial(noel,npt,js,ki,rc,t,dt,mu,e,de,dq,f0,f1,b,a,ierr,diagnostics)
  integer,optional,intent(in)::diagnostics(8)
 integer,intent(in)::noel,npt,js(4),ki,rc
 real(real64),intent(in)::t(2),dt,mu,e(6),de(6),dq(3,3),f0(3,3),f1(3,3),b(105),a(105)
 integer,intent(out)::ierr
 ! Wrapper ensures every internal early return still releases the mutex.
 call MutexLock(1)
  if(present(diagnostics))then
   if(rc==0)call record_tangent_observations(diagnostics,noel,npt,js(1),ki)
  endif
 call write_corot_trial_unlocked(noel,npt,js,ki,rc,t,dt,mu,e,de,dq,f0,f1,b,a,ierr)
 call MutexUnlock(1)
 end subroutine
 subroutine write_corot_trial_unlocked(noel,npt,js,ki,rc,t,dt,mu,e,de,dq,f0,f1,b,a,ierr)
 integer,intent(in)::noel,npt,js(4),ki,rc
 real(real64),intent(in)::t(2),dt,mu,e(6),de(6),dq(3,3),f0(3,3),f1(3,3),b(105),a(105)
 integer,intent(out)::ierr
 integer::u,n,st,cl,j,process_count,process_rank
 character(15)::rank_suffix
 character(80)::mapped_part
 integer::mapped_local,mapped_back,mapped_rc
 character(32)::value
 if(.not.initialized)then
  initialized=.true.
  call get_environment_variable('CDPM2_VISCO_LOG_LIMIT',value,length=n,status=st)
  if(st==1)then
   limit=0
  else if(st/=0.or.n<1.or.n>len(value))then
   config_error=1
  else
   do j=1,n
    if(value(j:j)<'0'.or.value(j:j)>'9')config_error=1
   end do
   if(config_error==0)then
    read(value(:n),*,iostat=st)limit
    if(st/=0)config_error=1
    if(limit<0.or.limit>10000)config_error=1
   end if
  end if
  if(config_error==0.and.limit>0)then
   call get_environment_variable('CDPM2_VISCO_LOG_NOEL',value,length=n,status=st)
   if(st/=1)then
    if(st/=0.or.n<1.or.n>len(value))then
     config_error=1
    else
     do j=1,n
      if(value(j:j)<'0'.or.value(j:j)>'9')config_error=1
     enddo
     if(config_error==0)then
      read(value(:n),*,iostat=st)target_noel
      if(st/=0)then
       config_error=1
      else if(target_noel<1)then
       config_error=1
      endif
     endif
    endif
   endif
  endif
  if(config_error==0.and.limit>0)then
   call get_environment_variable('CDPM2_VISCO_LOG_MODE',value,length=n,status=st)
   if(st==1)then
    failures_only=.false.
   else if(st/=0.or.n<1.or.n>len(value))then
    config_error=1
   else
    select case(value(:n))
    case('ALL')
     failures_only=.false.
    case('FAILURES')
     failures_only=.true.
    case default
     config_error=1
    end select
   endif
  endif
  if(config_error==0.and.limit>0)then
   call get_environment_variable('CDPM2_VISCO_LOG_PATH',path,length=n,status=st)
   if(st/=0.or.n<1.or.n>len(path))config_error=1
   if(config_error==0)then
    ! Abaqus2026 GETNUMCPUS means MPI processes, not OpenMP threads.
    call GETNUMCPUS(process_count)
    call GETRANK(process_rank)
    if(process_count<1.or.process_rank<0.or.process_rank>=process_count)then
     config_error=1
    else if(process_count>1)then
     write(rank_suffix,'(a,i10.10)',iostat=st)'.rank',process_rank
     if(st/=0.or.n+len(rank_suffix)>len(path))then
      config_error=1
     else
      path=path(:n)//rank_suffix
     endif
    endif
   endif
  end if
 end if
 ierr=config_error
 if(ierr/=0.or.limit==0.or.written>=limit)return
 ! Unselected elements do not consume the per-rank record budget.
 if(target_noel>0.and.noel/=target_noel)return
 if(target_noel>0.and..not.mapping_done)then
  mapping_done=.true.
  call cdpm2_map_element(noel,target_noel,mapped_part,mapped_local,mapped_back,mapped_rc)
  write(6,*)'CDPM2_MAP:',noel,mapped_local,mapped_back,mapped_rc,trim(mapped_part)
  if(mapped_rc/=0)then
   config_error=1;ierr=1;return
  endif
 endif
 ! Successful calls must not consume the bounded failure record budget.
 if(failures_only.and.rc==0.and.a(105)>=1.0_real64)return
 !1034: identify every recorded all-element failure, not every successful call.
 if(target_noel==0)then
  call GETPARTINFO(noel,1,mapped_part,mapped_local,mapped_rc)
  mapped_back=-1
  if(mapped_rc==0)then
   call GETINTERNAL(mapped_part,mapped_local,1,mapped_back,mapped_rc)
   if(mapped_rc==0.and.mapped_back/=noel)mapped_rc=4
  endif
  write(6,*)'CDPM2_FAILURE_MAP:',written+1,noel,npt,mapped_local,mapped_back,mapped_rc,trim(mapped_part)
  ! Keep numeric failure even if lookup fails; interpretation must flag mapping.
 endif
 open(newunit=u,file=trim(path),status='unknown',position='append',action='write',iostat=ierr)
 if(ierr/=0)return
 ! Explicit round-trip precision; schema1051 fixed262numeric fields, one physical record.
 write(u,'(9(i12,1x),253(es26.17e3,1x))',iostat=ierr)1051,noel,npt,js,ki,rc,t,dt,mu,e,de,dq,f0,f1,b,a
 close(u,iostat=cl)
 if(cl/=0)ierr=cl
 if(ierr==0)written=written+1
 end subroutine
end module

module work_diagnostic_435
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine update_work(old_work,sold,snew,de,total,plastic,new_work,rc)
 real(real64),intent(in)::old_work(3),sold(6),snew(6),de(6),total(6),plastic(6)
 real(real64),intent(inout)::new_work(3)
 integer,intent(out)::rc
 real(real64)::w(3)
 rc=2
 if(.not.all(ieee_is_finite(old_work)).or..not.all(ieee_is_finite(sold)).or. &
 .not.all(ieee_is_finite(snew)).or..not.all(ieee_is_finite(de)).or. &
 .not.all(ieee_is_finite(total)).or..not.all(ieee_is_finite(plastic)))return
 if(all(de==0))then
  new_work=old_work;rc=0;return
 end if
 w(1)=old_work(1)+dot_product((sold+snew)/2,de)
 w(2)=dot_product(snew,total-plastic)/2
 w(3)=w(1)-w(2)
 if(.not.all(ieee_is_finite(w)))return
 new_work=w;rc=0
 end subroutine
end module

module viscosity_1049
 use cdpm2_return,only:material_parameters
 use cdpm2_split,only:damaged_stress
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::relaxation,viscous_stress
contains
 subroutine relaxation(mu,dt,old,target,new,rc)
 real(real64),intent(in)::mu,dt,old(8),target(8)
 real(real64),intent(inout)::new(8)
 integer,intent(out)::rc
 real(real64)::a,b,r,tmp(8)
 rc=2
 if(.not.all(ieee_is_finite([mu,dt,old,target])))return
 if(mu<0.or.dt<0)return
 if(any(old(7:8)<0).or.any(old(7:8)>1))return
 if(any(target(7:8)<old(7:8)).or.any(target(7:8)>1))return
 if(mu==0)then
  tmp=target
 else if(dt==0)then
  tmp=old
 else
  if(dt>=mu)then
   r=mu/dt;a=1/(1+r);b=r/(1+r)
  else
   r=dt/mu;a=r/(1+r);b=1/(1+r)
  endif
  tmp=b*old+a*target
  ! Enforce exact irreversibility against floating-point convex-sum roundoff.
  tmp(7:8)=max(old(7:8),min(target(7:8),tmp(7:8)))
 endif
 if(.not.all(ieee_is_finite(tmp)))return
 new=tmp;rc=0
 end subroutine

 subroutine viscous_stress(p,e,v,s,rc)
 type(material_parameters),intent(in)::p
 real(real64),intent(in)::e(6),v(8)
 real(real64),intent(inout)::s(6)
 integer,intent(out)::rc
 real(real64)::ee(6),eff(6),out(6),lam
 rc=2
 if(.not.all(ieee_is_finite([e,v,p%bulk,p%shear])))return
 if(min(p%bulk,p%shear)<=0.or.any(v(7:8)<0).or.any(v(7:8)>1))return
 ee=e-v(1:6);lam=p%bulk-2*p%shear/3
 eff(1:3)=2*p%shear*ee(1:3)+lam*sum(ee(1:3))
 eff(4:6)=p%shear*ee(4:6)
 out=0
 call damaged_stress(eff,v(7),v(8),out,rc)
 if(rc==0)s=out
 end subroutine

end module

module endpoint_1076
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 use cdpm2_model,only:material_state
 use cdpm2_return,only:material_parameters
 use cdpm2_plastic,only:plastic_update,elastic_stress,elastic_strain,invariants
 use cdpm2_functions,only:regular_jacobian,regular_residual,flow_gradient,hardening_two
 use cdpm2_damage_history,only:damage_parameters,damage_history,update_damage
 use cdpm2_damage,only:equivalent_strain
 use cdpm2_split,only:split_stress
 use viscosity_1049,only:relaxation,viscous_stress
 use state_geometry_150,only:pack_geometry,unpack_geometry
 implicit none
 private
 public::endpoint_update,plastic_sensitivity
 type sensitivity
  real(real64)::ep(6,6)=0,k(6)=0,eff(6,6)=0,hist(13,6)=0
 end type
contains
 subroutine solve46(a,b,x,ok)
 real(real64),intent(in)::a(4,4),b(4,6)
 real(real64),intent(out)::x(4,6)
 logical,intent(out)::ok
 real(real64)::z(4,10),row(10),s
 integer::i,j,k
 z(:,1:4)=a;z(:,5:10)=b;x=0;ok=.false.
 do i=1,4
  s=maxval(abs(z(i,1:4)));if(s<=tiny(s))return
  z(i,:)=z(i,:)/s
 enddo
 do i=1,4
  k=i-1+maxloc(abs(z(i:4,i)),dim=1)
  if(abs(z(k,i))<1e-14_real64)return
  row=z(i,:);z(i,:)=z(k,:);z(k,:)=row
  do j=i+1,4
   z(j,i:10)=z(j,i:10)-z(j,i)/z(i,i)*z(i,i:10)
  enddo
 enddo
 do i=4,1,-1
  x(i,:)=(z(i,5:10)-matmul(z(i,i+1:4),x(i+1:4,:)))/z(i,i)
 enddo
 ok=all(ieee_is_finite(x))
 end subroutine

 subroutine inv_sensitivity(s,ds,mean,rho,theta,dev,dm,dr,dt)
 real(real64),intent(in)::s(6),ds(6,6)
 real(real64),intent(out)::mean,rho,theta,dev(6),dm(6),dr(6),dt(6)
 real(real64)::j3,c,gdet(6),dd(6,6),dc(6)
 integer::j
 call invariants(s,mean,rho,theta,dev)
 dm=sum(ds(1:3,:),dim=1)/3;dr=0;dt=0
 if(rho>1e-12_real64*max(1._real64,maxval(abs(s))))then
  dr=matmul(dev(1:3),ds(1:3,:))/rho+2*matmul(dev(4:6),ds(4:6,:))/rho
  j3=dev(1)*dev(2)*dev(3)+2*dev(4)*dev(5)*dev(6)-dev(1)*dev(6)**2-dev(2)*dev(5)**2-dev(3)*dev(4)**2
  c=max(-1._real64,min(1._real64,3*sqrt(6._real64)*j3/rho**3))
  if(1-c*c>1e-12_real64)then
   gdet=[dev(2)*dev(3)-dev(6)**2,dev(1)*dev(3)-dev(5)**2,dev(1)*dev(2)-dev(4)**2, &
    2*(dev(5)*dev(6)-dev(3)*dev(4)),2*(dev(4)*dev(6)-dev(2)*dev(5)),2*(dev(4)*dev(5)-dev(1)*dev(6))]
   dd=ds;do j=1,3;dd(j,:)=dd(j,:)-dm;enddo
   dc=3*sqrt(6._real64)/rho**3*(matmul(gdet,dd)-3*j3/rho*dr)
   dt=-dc/(3*sqrt(1-c*c))
  endif
 endif
 end subroutine

 subroutine plastic_sensitivity(p,e,epold,kold,s,ep,k,branch,de,depold,dkold,ds,dep,dk,fallback,details)
 integer,optional,intent(out)::details(3)
 type(material_parameters),intent(in)::p
 real(real64),intent(in)::e(6),epold(6),kold,s(6),ep(6),k,de(6,6),depold(6,6),dkold(6)
 integer,intent(in)::branch
 real(real64),intent(out)::ds(6,6),dep(6,6),dk(6)
 integer,intent(out)::fallback
 real(real64)::tr(6),dtr(6,6),dv(6),p0,r0,th,pn,rn,tn,v(6),dp0(6),dr0(6),dth(6)
 real(real64)::x(4),jac(4,4),rhs(4,6),dx(4,6),dg(2),rplus(4),rminus(4),rt(4),h,a,b
 real(real64)::sp(6),sm(6),epp(6),epm(6),kp,km,hs,ee(6),eoldp(6),eoldm(6)
 integer::j,rp,rm,bp,bm,events(3)
 logical::ok
 events=0;if(present(details))details=0
 fallback=0;tr=elastic_stress(p,e-epold)
 do j=1,6;dtr(:,j)=elastic_stress(p,de(:,j)-depold(:,j));enddo
 ds=dtr;dep=depold;dk=dkold
 if(branch==0)return
 if(branch==2)then
  ! Rare apex branch: differentiate only this one local return, not full history.
  do j=1,6
   hs=1e-9_real64/max(1._real64,maxval(abs(de(:,j))),maxval(abs(depold(:,j))),abs(dkold(j))*1e-4_real64)
   sp=s;sm=s;epp=ep;epm=ep;kp=k;km=k
   eoldp=epold+hs*depold(:,j);eoldm=epold-hs*depold(:,j)
   call plastic_update(p,e+hs*de(:,j),eoldp,max(0._real64,kold+hs*dkold(j)),sp,epp,kp,rp,bp)
   call plastic_update(p,e-hs*de(:,j),eoldm,max(0._real64,kold-hs*dkold(j)),sm,epm,km,rm,bm)
   if(rp==0.and.rm==0.and.bp==2.and.bm==2)then
    ds(:,j)=(sp-sm)/(2*hs);dep(:,j)=(epp-epm)/(2*hs);dk(j)=(kp-km)/(2*hs)
    events(1)=events(1)+1
   else if(rp==0.and.bp==2)then
    ds(:,j)=(sp-s)/hs;dep(:,j)=(epp-ep)/hs;dk(j)=(kp-k)/hs;fallback=fallback+1
    events(2)=events(2)+1
   else if(rm==0.and.bm==2)then
    ds(:,j)=(s-sm)/hs;dep(:,j)=(ep-epm)/hs;dk(j)=(k-km)/hs;fallback=fallback+1
    events(2)=events(2)+1
   else
    fallback=fallback+1
     events(3)=events(3)+1
   endif
  enddo
   if(present(details))details=events
   return
  endif
  call inv_sensitivity(tr,dtr,p0,r0,th,dv,dp0,dr0,dth)
 call invariants(s,pn,rn,tn,v)
 if(r0<=tiny(r0))then;fallback=1;return;endif
 dg=flow_gradient(pn,rn,k,p%fc,p%ft,p%ecc,p%q0,p%hp,p%dilation)
 if(abs(dg(2))<=tiny(r0))then;fallback=1;return;endif
 x=[pn,rn,k,(r0-rn)/(2*p%shear*dg(2))]
 call regular_jacobian(x,p%bulk,p%shear,th,p%fc,p%ft,p%ecc,p%q0,p%hp,p%dilation, &
  p%ah,p%bh,p%ch,p%dh,jac,ok)
 if(.not.ok)then;fallback=1;return;endif
 h=1e-6_real64;a=min(acos(-1._real64)/3,th+h);b=max(0._real64,th-h)
 rplus=regular_residual(x,[p0,r0],kold,p%bulk,p%shear,a,p%fc,p%ft,p%ecc,p%q0,p%hp, &
  p%dilation,p%ah,p%bh,p%ch,p%dh)
 rminus=regular_residual(x,[p0,r0],kold,p%bulk,p%shear,b,p%fc,p%ft,p%ecc,p%q0,p%hp, &
  p%dilation,p%ah,p%bh,p%ch,p%dh)
 rt=(rplus-rminus)/(a-b)
 do j=1,6;rhs(:,j)=[dp0(j),dr0(j),-dkold(j),0._real64]-rt*dth(j);enddo
 call solve46(jac,rhs,dx,ok)
 if(.not.ok)then;fallback=1;return;endif
 do j=1,6
  v=dtr(:,j);v(1:3)=v(1:3)-dp0(j)
  ds(:,j)=rn/r0*v+dv*(dx(2,j)/r0-rn*dr0(j)/r0**2)
  ds(1:3,j)=ds(1:3,j)+dx(1,j)
  dep(:,j)=de(:,j)-elastic_strain(p,ds(:,j))
 enddo
 dk=dx(3,:)
 end subroutine

 subroutine equiv_alpha(p,d,s,eq,alpha,rc)
 type(material_parameters),intent(in)::p
 type(damage_parameters),intent(in)::d
 real(real64),intent(in)::s(6)
 real(real64),intent(out)::eq,alpha
 integer,intent(out)::rc
 real(real64)::m,r,t,v(6),sp(6),sm(6)
 call invariants(s,m,r,t,v)
 eq=equivalent_strain(m,r,t,p%fc,p%ft,p%ecc,d%young)
 alpha=0;call split_stress(s,sp,sm,alpha,rc)
 end subroutine

 subroutine equiv_sensitivity(p,d,s,ds,eq,alpha,dq,da)
 type(material_parameters),intent(in)::p
 type(damage_parameters),intent(in)::d
 real(real64),intent(in)::s(6),ds(6,6)
 real(real64),intent(out)::eq,alpha,dq(6),da(6)
 real(real64)::qp,qm,ap,am,h,m,r,t,dev(6),dm(6),dr(6),dt(6),sp(6),sm(6),den
 real(real64)::mm,shape,shape_p,shape_m,ta,tb,bb,db(6),root
 integer::j,rc
 call inv_sensitivity(s,ds,m,r,t,dev,dm,dr,dt)
 eq=equivalent_strain(m,r,t,p%fc,p%ft,p%ecc,d%young)
 call split_stress(s,sp,sm,alpha,rc);dq=0;da=0
 den=sum(s(1:3)**2)+2*sum(s(4:6)**2)
 if(den>tiny(den))da=2/den*(matmul(sm(1:3)-alpha*s(1:3),ds(1:3,:)) &
  +2*matmul(sm(4:6)-alpha*s(4:6),ds(4:6,:)))
 mm=3*(p%fc**2-p%ft**2)/(p%fc*p%ft)*p%ecc/(p%ecc+1)
 shape=section_shape(t,p%ecc);ta=min(acos(-1._real64)/3,t+1e-6_real64);tb=max(0._real64,t-1e-6_real64)
 shape_p=section_shape(ta,p%ecc);shape_m=section_shape(tb,p%ecc)
 bb=mm*(r*shape/(sqrt(6._real64)*p%fc)+m/p%fc)
 db=mm*((shape*dr+r*(shape_p-shape_m)/(ta-tb)*dt)/(sqrt(6._real64)*p%fc)+dm/p%fc)
 root=sqrt(bb**2/4+1.5_real64*(r/p%fc)**2)
 if(root>1e-14_real64)then
  dq=p%ft/d%young*((.5_real64+bb/(4*root))*db+1.5_real64*r/p%fc**2/root*dr)
 else
  ! At the zero-stress corner use bounded local directional differences.
  do j=1,6
   if(maxval(abs(ds(:,j)))==0)cycle
   h=1e-6_real64/maxval(abs(ds(:,j)))
   call equiv_alpha(p,d,s+h*ds(:,j),qp,ap,rc)
   call equiv_alpha(p,d,s-h*ds(:,j),qm,am,rc)
   dq(j)=(qp-qm)/(2*h)
  enddo
 endif
 end subroutine

 pure real(real64) function section_shape(t,ecc) result(r)
 real(real64),intent(in)::t,ecc
 real(real64)::c
 c=cos(t)
 r=(4*(1-ecc**2)*c**2+(2*ecc-1)**2)/(2*(1-ecc**2)*c+(2*ecc-1)*sqrt(4*(1-ecc**2)*c**2+5*ecc**2-4*ecc))
 end function

 subroutine damage_sensitivity(p,d,old,new,go,gn,ge,fallback)
 type(material_parameters),intent(in)::p
 type(damage_parameters),intent(in)::d
 type(material_state),intent(in)::old,new
 type(sensitivity),intent(in)::go
 type(sensitivity),intent(inout)::gn
 real(real64),intent(in)::ge(6,6)
 integer,intent(out)::fallback
 real(real64)::eq,alpha,dq(6),da(6),qp,qm,ap,dummy(6),qmgrad(6),emin,demin(6),mid,dm(6)
 real(real64)::pmean,rho,theta,dev(6),dsig(6),drho(6),dtheta(6),rs,drs(6),duct,dd(6)
 real(real64)::pn,dpn(6),pn0,dpn0(6),deplast(6),gdplast(6,6),delta,ddelta(6),e0,tol
 real(real64)::fac,dfac(6),ex,dex(6),q2,dq2(6),w,wo,den,traction,slope,L,kap1,kap2,xx(6)
 integer::i,j,rc
 logical::reversal
 fallback=0;gn%hist=go%hist;e0=p%ft/d%young;tol=d%damage_tolerance
 call equiv_sensitivity(p,d,new%effective,gn%eff,eq,alpha,dq,da)
 call equiv_alpha(p,d,old%effective+.01_real64*(new%effective-old%effective),qp,ap,rc)
 call equiv_alpha(p,d,old%effective+.99_real64*(new%effective-old%effective),qm,ap,rc)
 reversal=old%history%equivalent>qp.and.eq>qm.and.abs(qp-old%history%equivalent)>tol/100.and.abs(qm-eq)>tol/100
 emin=old%history%equivalent;demin=go%hist(2,:)
 if(reversal)then
  do i=1,100
   fac=real(i,real64)/100
   call equiv_sensitivity(p,d,old%effective+fac*(new%effective-old%effective), &
    go%eff+fac*(gn%eff-go%eff),mid,ap,dm,dummy)
   if(mid>emin)exit
   emin=mid;demin=dm
  enddo
 endif
 gn%hist(2,:)=dq;gn%hist(5,:)=da
 gn%hist(3,:)=go%hist(3,:)+dq-go%hist(2,:)
 if(reversal)then
  gn%hist(4,:)=go%hist(4,:)+go%hist(5,:)*(emin-old%history%equivalent) &
   +old%history%alpha*(demin-go%hist(2,:))+da*(eq-emin)+alpha*(dq-demin)
 else
  gn%hist(4,:)=go%hist(4,:)+da*(eq-old%history%equivalent)+alpha*(dq-go%hist(2,:))
 endif
 call inv_sensitivity(new%effective,gn%eff,pmean,rho,theta,dev,dsig,drho,dtheta)
 rs=0;drs=0
 if(pmean<0)then
  rs=-pmean/(.40824829_real64*max(rho,1e-16_real64))
  drs=-dsig/(.40824829_real64*max(rho,1e-16_real64))
  if(rho>1e-16_real64)drs=drs-rs*drho/rho
 endif
 duct=1+(d%asoft-1)*rs**d%bsoft;dd=0
 if(rs>0)dd=(d%asoft-1)*d%bsoft*rs**(d%bsoft-1)*drs
 deplast=new%plastic-old%plastic;gdplast=gn%ep-go%ep
 pn0=norm2([deplast(1:3),deplast(4:6)/sqrt(2._real64)]);dpn0=0
 if(pn0>tiny(pn0))dpn0=(matmul(deplast(1:3),gdplast(1:3,:))+.5_real64*matmul(deplast(4:6),gdplast(4:6,:)))/pn0
 do j=1,2
  delta=new%history%accumulated(j)-old%history%threshold(j)
  if(delta/e0 < -tol)cycle
  ddelta=gn%hist(2+j,:)-go%hist(5+j,:)
  gn%hist(5+j,:)=gn%hist(2+j,:)
  pn=pn0;dpn=dpn0
  if(j==2)then;pn=alpha*pn0;dpn=alpha*dpn0+pn0*da;endif
  if(new%history%threshold(j)<=e0*(1-tol))then
   pn=0;dpn=0
  else if(old%history%threshold(j)<e0*(1-tol))then
   fac=1-(e0*(1-tol)-old%history%threshold(j))/(new%history%threshold(j)-old%history%threshold(j))
   dfac=go%hist(5+j,:)/delta+(e0*(1-tol)-old%history%threshold(j))*ddelta/delta**2
   dpn=abs(fac)*dpn+sign(1._real64,fac)*pn*dfac;pn=abs(fac)*pn
  endif
  if(j==2)then
   q2=hardening_two(new%kappa,p%hp);dq2=0;if(new%kappa>=1)dq2=p%hp*gn%k
   ex=p%ft*q2*sqrt(2._real64/3)/(max(rho,1e-16_real64)*sqrt(1+2*p%dilation**2))
   dex=ex*dq2/q2;if(rho>1e-16_real64)dex=dex-ex*drho/rho
   dpn=ex*dpn+pn*dex;pn=pn*ex
  endif
  gn%hist(7+j,:)=go%hist(7+j,:)+dpn/duct-pn*dd/duct**2
  gn%hist(9+j,:)=go%hist(9+j,:)+ddelta/duct-delta*dd/duct**2
  w=new%history%damage(j);wo=old%history%damage(j)
  if(w==wo)cycle
  if(w==1)then;gn%hist(11+j,:)=0;cycle;endif
  kap1=new%history%one(j);kap2=new%history%two(j)
  if(j==1)then
   L=d%le;if(d%automatic_length)L=new%history%selected_le
   slope=(d%ft1-p%ft*(1-tol))/d%wf1
   if(L*(kap1+w*kap2)>d%wf1+64*epsilon(1._real64)*d%wf)slope=-d%ft1/(d%wf-d%wf1)
   den=d%young*new%history%threshold(j)+slope*L*kap2
   xx=(1-w)*d%young*gn%hist(5+j,:)-slope*L*(gn%hist(7+j,:)+w*gn%hist(9+j,:))
  else
   traction=p%ft*(1-tol)*exp(-(kap1+w*kap2)/d%efc)
   den=d%young*new%history%threshold(j)-traction*kap2/d%efc
   xx=(1-w)*d%young*gn%hist(5+j,:)+traction/d%efc*(gn%hist(7+j,:)+w*gn%hist(9+j,:))
  endif
  if(abs(den)>1e-12_real64*p%ft)then
   gn%hist(11+j,:)=xx/den
  else
   fallback=fallback+1
  endif
 enddo
 ! Selected length is carried as state, but its tangent is intentionally frozen.
 gn%hist(1,:)=0
 end subroutine

 subroutine split_sensitivity(s,ds,sp,sm,dsp,dsm,rc,hfd,local_fd)
 logical,optional,intent(out)::local_fd
 real(real64),intent(in)::s(6),ds(6,6)
 real(real64),intent(out)::sp(6),sm(6),dsp(6,6),dsm(6,6)
 integer,intent(out)::rc
 real(real64),intent(in)::hfd
 real(real64)::alpha,a(3,3),v(3,3),rot(3,3),b(3,3),f(3,3),lam(3),slopes(3),scale,phi,cs,sn,tol
 real(real64)::pp(6),pm(6),mp(6),mm(6),band
 integer::j,i,k,l,ip,iq,sweep
  if(present(local_fd))local_fd=.false.
  call split_stress(s,sp,sm,alpha,rc);if(rc/=0)return
 scale=maxval(abs(s));v=0
 do i=1,3;v(i,i)=1;enddo
 if(scale<=tiny(scale))then
  dsp=.5_real64*ds;dsm=ds-dsp;return
 endif
 a(1,:)=[s(1),s(4),s(5)]/scale;a(2,:)=[s(4),s(2),s(6)]/scale;a(3,:)=[s(5),s(6),s(3)]/scale
 do sweep=1,60
  ip=1;iq=2
  do i=1,2;do k=i+1,3
   if(abs(a(i,k))>abs(a(ip,iq)))then;ip=i;iq=k;endif
  enddo;enddo
  if(abs(a(ip,iq))<8*epsilon(1._real64))exit
  phi=.5_real64*atan2(2*a(ip,iq),a(iq,iq)-a(ip,ip));cs=cos(phi);sn=sin(phi)
  rot=0;do i=1,3;rot(i,i)=1;enddo
  rot(ip,ip)=cs;rot(iq,iq)=cs;rot(ip,iq)=sn;rot(iq,ip)=-sn
  a=matmul(transpose(rot),matmul(a,rot));v=matmul(v,rot)
 enddo
 if(sweep>60)then;rc=1;return;endif
 lam=[a(1,1),a(2,2),a(3,3)];slopes=.5_real64;tol=1e-12_real64
 ! Zero-principal-stress kink: centered local spectral derivative in the same
 ! strain probe width as the prior algorithm. Never reintegrates plasticity.
 ! This avoids selecting a side from insignificant lateral stress noise.
 band=hfd*maxval(abs(ds))/scale
 if(minval(abs(lam))<=4*band)then
   if(present(local_fd))local_fd=.true.
  do j=1,6
   call split_stress(s+hfd*ds(:,j),pp,mp,alpha,rc);if(rc/=0)return
   call split_stress(s-hfd*ds(:,j),pm,mm,alpha,rc);if(rc/=0)return
   dsp(:,j)=(pp-pm)/(2*hfd)
  enddo
  dsm=ds-dsp;return
 endif
 where(lam>tol)slopes=1
 where(lam< -tol)slopes=0
 do i=1,3;do k=1,3
  if(abs(lam(i)-lam(k))>tol)then
   f(i,k)=(max(lam(i),0._real64)-max(lam(k),0._real64))/(lam(i)-lam(k))
  else
   f(i,k)=.5_real64*(slopes(i)+slopes(k))
  endif
 enddo;enddo
 do j=1,6
  b(1,:)=[ds(1,j),ds(4,j),ds(5,j)];b(2,:)=[ds(4,j),ds(2,j),ds(6,j)];b(3,:)=[ds(5,j),ds(6,j),ds(3,j)]
  b=matmul(transpose(v),matmul(b,v));b=f*b;b=matmul(v,matmul(b,transpose(v)))
  dsp(:,j)=[b(1,1),b(2,2),b(3,3),b(1,2),b(1,3),b(2,3)]
 enddo
 dsm=ds-dsp
 end subroutine

 subroutine endpoint_update(p,d,old,oldv,e,mu,dt,cap,budget,hfd,new,newv,s,c,rc,stats,diagnostics)
 integer,optional,intent(out)::diagnostics(8)
 type(material_parameters),intent(in)::p
 type(damage_parameters),intent(in)::d
 type(material_state),intent(in)::old
 real(real64),intent(in)::oldv(8),e(6),mu,dt,cap,hfd
 integer,intent(in)::budget
 type(material_state),intent(inout)::new
 real(real64),intent(inout)::newv(8),s(6),c(6,6)
 integer,intent(out)::rc
 integer,optional,intent(out)::stats(5)
 type(material_state)::cur,tr,checked
 type(sensitivity)::g,gt
 real(real64)::delta(6),length,ge(6,6),id(6,6),v(8),b,dv(8,6)
 real(real64)::sv(6),effective(6),ceff(6,6),sp(6),sm(6),dsp(6,6),dsm(6,6),ct(6,6)
 real(real64)::packed(39)
 integer::br,lr,j,nacc,nrej,npf,ndf,nvertex,fb,diag(8),pd(3)
 logical::local_fd
 rc=2;if(present(stats))stats=0
 diag=0;if(present(diagnostics))diagnostics=0
 if(.not.all(ieee_is_finite([e,old%strain,oldv,mu,dt,cap,hfd])))return
 if(min(mu,dt)<0.or.cap<=0.or.budget<1.or.hfd<=0)return
 if(any(oldv(7:8)<0).or.any(oldv(7:8)>old%history%damage))return
 call pack_geometry(old,packed,rc);if(rc/=0)return
 call unpack_geometry(packed,checked,rc);if(rc/=0)return
 rc=2;cur=checked;g=sensitivity();delta=e-old%strain
 length=norm2([delta(1:3),delta(4:6)/sqrt(2._real64)])
 id=0;do j=1,6;id(j,j)=1;enddo
 nacc=0;nrej=0;npf=0;ndf=0;nvertex=0
 if((dt>0.or.mu==0).and.length>tiny(length))then
  tr=cur;tr%strain=e
  call plastic_update(p,e,cur%plastic,cur%kappa,tr%effective,tr%plastic,tr%kappa,lr,br,cur%effective)
  if(lr==0)call update_damage(p,d,cur%history,cur%effective,tr%effective,tr%plastic-cur%plastic, &
   tr%kappa,tr%history,tr%stress,lr,total_strain=e)
  if(lr/=0)then
   rc=1;if(lr==2)rc=2
   if(present(stats))stats=[0,1,0,0,0]
   return
  endif
  !1133: numerical increment-control limit, not a damage ceiling.
  ! Reject an abrupt trial; all state stays committed until a smaller increment
  ! is accepted. Constitutive damage may still reach1 over accepted increments.
  if(maxval(tr%history%damage-cur%history%damage)>.25_real64)then
   rc=1
   if(present(stats))stats=[0,1,0,0,0]
   return
  endif
  ! Committed old state is held fixed. Derivative of endpoint strain is identity.
  gt=g;ge=id
  call plastic_sensitivity(p,e,cur%plastic,cur%kappa,tr%effective,tr%plastic,tr%kappa,br, &
   ge,g%ep,g%k,gt%eff,gt%ep,gt%k,fb,pd)
  diag(1:3)=pd
  if(br/=2)diag(4)=fb
  npf=fb;if(br==2)nvertex=1
  call damage_sensitivity(p,d,cur,tr,g,gt,ge,fb);ndf=fb
  diag(5)=fb
  if(d%automatic_length)diag(8)=1
  cur=tr;g=gt;nacc=1
 else if(dt>0)then
  ! Zero increment: practical frozen-internal-variable query tangent.
  g=sensitivity()
  diag(7)=1
 endif
 if(mu==0)then
  b=1
 else if(dt==0)then
  b=0
 else if(dt>=mu)then
  b=1/(1+mu/dt)
 else
  b=(dt/mu)/(1+dt/mu)
 endif
 v=oldv
 call relaxation(mu,dt,oldv,[cur%plastic,cur%history%damage],v,rc);if(rc/=0)return
 dv(1:6,:)=b*g%ep;dv(7:8,:)=b*g%hist(12:13,:)
 effective=elastic_stress(p,e-v(1:6))
 do j=1,6;ceff(:,j)=elastic_stress(p,id(:,j)-dv(1:6,j));enddo
 call split_sensitivity(effective,ceff,sp,sm,dsp,dsm,rc,hfd,local_fd);if(rc/=0)return
 if(local_fd)diag(6)=1
 call viscous_stress(p,e,v,sv,rc);if(rc/=0)return
 do j=1,6
  ct(:,j)=(1-v(7))*dsp(:,j)+(1-v(8))*dsm(:,j)-sp*dv(7,j)-sm*dv(8,j)
 enddo
 rc=2;if(.not.all(ieee_is_finite([sv,reshape(ct,[36])])) )return
 new=cur;newv=v;s=sv;c=ct;rc=0
 if(present(stats))stats=[nacc,nrej,nvertex,npf,ndf]
 if(present(diagnostics))diagnostics=diag
 end subroutine
end module

module endpoint_payload_1079
 use endpoint_1076,only:endpoint_update
 use viscosity_1049,only:viscous_stress
 use corot_core_406,only:decode_properties,decode_initial_geometry,stage_frame,volume_jacobian
 use rotation_281,only:tensor
 use cdpm2_plastic,only:elastic_stress
 use cdpm2_model,only:material_state
 use cdpm2_return,only:material_parameters
 use cdpm2_damage_history,only:damage_parameters
 use state_geometry_150,only:pack_geometry,unpack_geometry
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 implicit none
contains
 subroutine endpoint_payload(props,mu,dt,e,de,drot,jdet,s,v,c,pnewdt,rc,diagnostics)
 integer,optional,intent(out)::diagnostics(8)
 real(real64),intent(in)::props(50),mu,dt,e(6),de(6),drot(3,3),jdet
 real(real64),intent(inout)::s(6),v(59),c(6,6),pnewdt
 integer,intent(out)::rc
 type(material_parameters)::p
 type(damage_parameters)::d
 type(material_state)::old,new
 real(real64)::out(59),st(6),ct(6,6),cj(6,6),q(3,3),qt(3,3),x0(3,8),xt(3,8)
 real(real64)::vt(39),iv(8),nv(8),bs(6),expected(6),cap,h,unit(6)
 integer::budget,i,diag(8)
 rc=2
 diag=0;if(present(diagnostics))diagnostics=0
 if(.not.all(ieee_is_finite([mu,dt,e,de,s,v,reshape(c,[36]),pnewdt,jdet])))return
 if(mu<0.or.dt<0.or.jdet<=0.or.pnewdt<=0)return
 if(any(v/=0))then
  if(v(58)/=1049.or.v(59)/=mu.or.v(40)/=386)return
 endif
 out=v;st=s;ct=c
 call decode_properties(props(1:23),p,d,cap,h,budget,rc);if(rc/=0)return
 call decode_initial_geometry(props(24:50),x0,d%mean_size,rc);if(rc/=0)return
 q=0;qt=0;xt=0;vt=0
 if(all(v==0))then
  do i=1,3
   q(i,i)=1
  enddo
 else
  q=reshape(v(41:49),[3,3])
 endif
 call tensor(drot,v(21:26),1._real64,bs)
 call stage_frame(drot,q,x0,v(1:39),bs,e,qt,xt,vt,rc);if(rc/=0)return
 iv=v(50:57);call tensor(drot,v(50:55),2._real64,iv(1:6))
 expected=0;call viscous_stress(p,e,iv,expected,rc);if(rc/=0)return
 rc=2
 if(maxval(abs(expected-s))>1e-10_real64*max(1._real64,maxval(abs(s))))return
 call unpack_geometry(vt,old,rc);if(rc/=0)return
 d%automatic_length=.true.;d%reference_nodes=xt;new=old;nv=iv
 call endpoint_update(p,d,old,iv,e+de,mu,dt,cap,budget,h,new,nv,st,ct,rc,diagnostics=diag)
 if(rc/=0)then
  if(rc==1.or.rc==3)then
   ! Rejected increment ONLY. Abaqus still assembles before honoring cutback.
   ! Incoming DDSDDE is an output workspace (observed zero), not old tangent.
   ! Return finite elastic predictor stiffness; do not advance stress/history.
   ! PNEWDT<1 requires rejection; this C is NEVER an accepted material tangent.
   pnewdt=min(pnewdt,.5_real64)
   do i=1,6
    unit=0;unit(i)=1;ct(:,i)=elastic_stress(p,unit)
   enddo
   call volume_jacobian(s,ct,cj,i)
   if(i/=0)then
    rc=2;return
   endif
   c=cj
  endif
  return
 endif
 if(mu==0.and.all(de==0))then
  ! Legacy mu=0 Jacobian-only call must not initialize the material marker.
  st=s;nv=[vt(8:13),vt(37:38)]
 else
  call pack_geometry(new,vt,rc);if(rc/=0)return
 endif
 call volume_jacobian(st,ct,cj,rc);if(rc/=0)return
 out(1:39)=vt;out(40)=386;out(41:49)=reshape(qt,[9]);out(50:57)=nv
 out(58)=1049;out(59)=mu
 s=st;v=out;c=cj;rc=0
 if(present(diagnostics))diagnostics=diag
 end subroutine
end module

! Automatic input geometry adapter. Constitutive calculations are unchanged.
module fil_geometry_20261011
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 use geometry_provider_421,only:install_geometry,clear_geometry_table
 implicit none
 private
 public::load_fil_geometry
 type ascii_reader
  integer::unit,pos=81
  character(80)::line=''
 end type
contains
 subroutine next_char(r,c,rc)
 type(ascii_reader),intent(inout)::r
 character,intent(out)::c
 integer,intent(out)::rc
 rc=0
 if(r%pos>80)then
  read(r%unit,'(a80)',iostat=rc)r%line
  if(rc/=0)return
  r%pos=1
 endif
 c=r%line(r%pos:r%pos);r%pos=r%pos+1
 end subroutine
 subroutine chars(r,n,s,rc)
 type(ascii_reader),intent(inout)::r
 integer,intent(in)::n
 character(*),intent(out)::s
 integer,intent(out)::rc
 integer::i
 s=''
 do i=1,n
  call next_char(r,s(i:i),rc)
  if(rc/=0)return
 enddo
 end subroutine
 subroutine word(r,kind,iv,rv,sv,rc)
 type(ascii_reader),intent(inout)::r
 character,intent(out)::kind
 integer,intent(out)::iv,rc
 real(real64),intent(out)::rv
 character(8),intent(out)::sv
 character(32)::s
 integer::w
 iv=0;rv=0;sv=''
 call next_char(r,kind,rc);if(rc/=0)return
 select case(kind)
 case('I')
  call chars(r,2,s,rc);if(rc/=0)return
  read(s,*,iostat=rc)w;if(rc/=0)return
  if(w<1.or.w>10)then
   rc=2;return
  endif
  call chars(r,w,s,rc);if(rc/=0)return
  read(s,*,iostat=rc)iv
 case('D')
  call chars(r,22,s,rc);if(rc/=0)return
  read(s,*,iostat=rc)rv
 case('A')
  call chars(r,8,sv,rc)
 case default
  rc=2
 end select
 end subroutine
 subroutine record(r,iv,rv,sv,kinds,n,rc)
 type(ascii_reader),intent(inout)::r
 integer,intent(out)::iv(500),n,rc
 real(real64),intent(out)::rv(500)
 character(8),intent(out)::sv(500)
 character,intent(out)::kinds(500)
 character::c
 integer::i
 iv=0;rv=0;sv='';kinds='';n=0
 do
  call next_char(r,c,rc);if(rc/=0)return
  if(c=='*')exit
  if(c/=' ')then
   rc=2;return
  endif
 enddo
 call word(r,kinds(1),iv(1),rv(1),sv(1),rc);if(rc/=0)return
 n=iv(1)
 if(kinds(1)/='I'.or.n<2.or.n>500)then
  rc=2;return
 endif
 do i=2,n
  call word(r,kinds(i),iv(i),rv(i),sv(i),rc);if(rc/=0)return
 enddo
 if(kinds(2)/='I')rc=2
 end subroutine
 integer function node_count(kind) result(nn)
 character(*),intent(in)::kind
 nn=0
 select case(trim(kind))
 ! C3D8I has the same eight corner nodes; incompatible modes are internal to Abaqus.
 case('C3D8','C3D8R','C3D8I');nn=8
 case('C3D6');nn=6
 end select
 end function
 subroutine load_fil_geometry(path,rc,why)
 character(*),intent(in)::path
 integer,intent(out)::rc
 character(*),intent(out)::why
 type(ascii_reader)::r
 integer::iv(500),n,ios,pass,ne,nv,na,expected_e,expected_n,ie,inode,nn,i,j,k,stat
 real(real64)::rv(500),x(3,8)
 character(8)::sv(500)
 character::kinds(500)
 integer,allocatable::ids(:),conn(:,:),nid(:)
 real(real64),allocatable::xyz(:,:),data(:,:)
 call clear_geometry_table()
 rc=2;why='ASCII FIL unavailable or incomplete; request *FILE FORMAT, ASCII and mesh-wide *EL FILE / *NODE FILE.'
 ne=0;nv=0;na=0;expected_e=-1;expected_n=-1
 do pass=1,2
  open(newunit=r%unit,file=path,status='old',action='read',iostat=ios)
  if(ios/=0)return
  r%pos=81;ie=0;inode=0
  do
   call record(r,iv,rv,sv,kinds,n,ios)
   if(ios/=0)exit
   if(iv(2)==2000)exit ! Initial mesh only; never read evolving result coordinates.
   select case(iv(2))
   case(1921)
    if(n<8)then
     ios=2;exit
    endif
    if(pass==1)then
     expected_e=iv(7);expected_n=iv(8)
    endif
   case(1900)
    if(n<4.or.kinds(3)/='I'.or.kinds(4)/='A')then
     ios=2;exit
    endif
    nn=node_count(sv(4))
    if(pass==1)na=na+1
    if(nn==0)cycle ! Non-supported elements are not installed for CDPM2.
    if(n/=nn+4.or.iv(3)<=0)then
     ios=2;exit
    endif
    if(any(kinds(5:n)/='I').or.any(iv(5:n)<=0))then
     ios=2;exit
    endif
    if(pass==1)then
     ne=ne+1
    else
     ie=ie+1
     if(ie>ne)then
      ios=2;exit
     endif
     ids(ie)=iv(3);conn(:,ie)=0;conn(:nn,ie)=iv(5:n)
     data(1:3,ie)=[2._real64,real(nn,real64),0._real64]
    endif
   case(1901)
    if(n/=6.or.kinds(3)/='I'.or.any(kinds(4:6)/='D').or.iv(3)<=0)then
     ios=2;exit
    endif
    if(.not.all(ieee_is_finite(rv(4:6))))then
     ios=2;exit
    endif
    if(pass==1)then
     nv=nv+1
    else
     inode=inode+1
     if(inode>nv)then
      ios=2;exit
     endif
     nid(inode)=iv(3);xyz(:,inode)=rv(4:6)
    endif
   end select
  enddo
  close(r%unit)
  if(ios>0)return
  if(pass==1)then
   ! Reject a partially flushed header, rather than silently installing a subset.
   if(ne<1.or.nv<1.or.na/=expected_e.or.nv/=expected_n)return
   allocate(ids(ne),conn(8,ne),nid(nv),xyz(3,nv),data(27,ne),stat=stat)
   if(stat/=0)then
    why='FIL geometry allocation failed';return
   endif
   data=0
  else
   if(ie/=ne.or.inode/=nv)return
  endif
 enddo
 call sort_nodes(1,nv)
 do i=2,nv
  if(nid(i)<=nid(i-1))then
   why='Duplicate FIL node labels';return
  endif
 enddo
 do i=1,ne
  nn=int(data(2,i));x=0
  do j=1,nn
   k=find_node(conn(j,i))
   if(k==0)then
    why='Missing FIL connectivity node';return
   endif
   x(:,j)=xyz(:,k)
  enddo
  data(4:27,i)=reshape(x,[24])
 enddo
 call install_geometry(ids,data,rc)
 if(rc==0)then
  why='ASCII FIL initial geometry installed'
 else
  why='Invalid or duplicate FIL element geometry'
 endif
contains
 recursive subroutine sort_nodes(lo,hi)
 integer,intent(in)::lo,hi
 integer::a,b,p,tmp
 real(real64)::v(3)
 if(lo>=hi)return
 a=lo;b=hi;p=nid(lo+(hi-lo)/2)
 do
  do while(nid(a)<p);a=a+1;enddo
  do while(nid(b)>p);b=b-1;enddo
  if(a<=b)then
   tmp=nid(a);nid(a)=nid(b);nid(b)=tmp
   v=xyz(:,a);xyz(:,a)=xyz(:,b);xyz(:,b)=v
   a=a+1;b=b-1
  endif
  if(a>b)exit
 enddo
 if(lo<b)call sort_nodes(lo,b)
 if(a<hi)call sort_nodes(a,hi)
 end subroutine
 integer function find_node(label) result(pos)
 integer,intent(in)::label
 integer::a,b,m
 a=1;b=nv;pos=0
 do while(a<=b)
  m=a+(b-a)/2
  if(nid(m)==label)then
   pos=m;return
  endif
  if(nid(m)<label)then
   a=m+1
  else
   b=m-1
  endif
 enddo
 end function
 end subroutine
end module

subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
  use tangent_observations,only:reset_tangent_observations,report_tangent_observations
 use geometry_provider_421,only:load_geometry_table,clear_geometry_table,geometry_count
 use fil_geometry_20261011,only:load_fil_geometry
 use iso_fortran_env,only:real64
 implicit none
 integer,intent(in)::lop,lrestart,kstep,kinc
 real(real64),intent(in)::time(2),dtime
 character(2048),save::path=''
 logical,save::instance_table_pending=.false.
 logical,save::fil_pending=.false.
 character(256),save::output_base=''
 character(256)::jobname,why
 character(80)::magic
 integer::rc,n,stat,u,ios
 if(lop==0.or.lop==4)then
  !1107: initialize geometry and logger mutex also on restart LOP4.
  ! Material STATEV is restored by Abaqus, not reset here.
  ! Mutex1 is reserved for the whole logger configuration/counter/file resource.
  call MutexInit(1)
   call reset_tangent_observations()
  call clear_geometry_table()
  instance_table_pending=.false.
  fil_pending=.false.
  call get_environment_variable('CDPM2_GEOMETRY_TABLE',path,length=n,status=stat)
  if(stat/=0.or.n==0.or.n>len(path))then
   if(lop==4)then
    write(6,*)'CDPM2 restart requires original CDPM2_GEOMETRY_TABLE; no automatic deck guessing.'
    return
   endif
   if(n>len(path).or.(stat/=0.and.stat/=1))then
    write(6,*)'CDPM2 invalid geometry table environment; UMAT will refuse update.'
    return
   endif
   call GETOUTDIR(output_base,n)
   if(n<1.or.n>len(output_base))return
   call get_environment_variable('CDPM2_FIL_FILE',path,length=n,status=stat)
   if(stat==1.or.n==0)then
    call GETJOBNAME(jobname,n)
    if(n<1.or.n>len(jobname))return
    path=trim(output_base)//'/'//trim(jobname)//'.fil'
   else if(stat/=0.or.n>len(path))then
    return
   endif
   fil_pending=.true.
   write(6,*)'CDPM2 ASCII FIL geometry pending first increment: ',trim(path)
   return
  end if
  !1114: GETPARTINFO returned default flattened names at LOP0 on Abaqus2026.
  ! Its actual instance mapping was correct at LOP1. Preserve the strict
  ! round-trip check; defer V2 once to increment start, before material calls.
  magic=''
  open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
  if(ios==0)then
   read(u,'(a)',iostat=ios)magic
   close(u)
  endif
  if(ios==0.and.trim(magic)=='CDPM2_GEOMETRY_V2')then
   instance_table_pending=.true.
   write(6,*)'CDPM2 instance geometry pending first increment initialization.'
   return
  endif
  call load_geometry_table(trim(path),rc)
  write(6,*)'CDPM2 geometry initialization:',rc,geometry_count()
 else if(lop==1.and.fil_pending)then
  fil_pending=.false.
  call load_fil_geometry(trim(path),rc,why)
  write(6,*)'CDPM2 ASCII FIL geometry:',rc,geometry_count(),trim(why)
 else if(lop==1.and.instance_table_pending)then
  instance_table_pending=.false.
  fil_pending=.false.
  call load_geometry_table(trim(path),rc)
  write(6,*)'CDPM2 instance geometry initialization:',rc,geometry_count()
 else if(lop==3)then
   call MutexLock(1)
   call report_tangent_observations()
   call MutexUnlock(1)
  ! End-of-analysis cleanup only. Restart LOP4 reloads immutable geometry.
  call clear_geometry_table()
  instance_table_pending=.false.
 end if
end subroutine

subroutine cdpm2_map_element(internal_id,expected_label,part_name,local_id,back_id,rc)
 implicit none
 integer,intent(in)::internal_id,expected_label
 character(80),intent(out)::part_name
 integer,intent(out)::local_id,back_id,rc
 integer::jrc
 part_name='';local_id=-1;back_id=-1;rc=1
 if(internal_id<1.or.expected_label<1)return
 call GETPARTINFO(internal_id,1,part_name,local_id,jrc)
 if(jrc/=0)return
 rc=2
 if(local_id/=expected_label)return
 call GETINTERNAL(part_name,local_id,1,back_id,jrc)
 rc=3
 if(jrc/=0)return
 rc=4
 if(back_id/=internal_id)return
 rc=0
end subroutine

subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,drplde,drpldt, &
 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,ndi,nshr,ntens,nstatv, &
 props,nprops,coords,drot,pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
 use endpoint_payload_1079,only:viscous_payload=>endpoint_payload,real64,ieee_is_finite
 use trial_log_391,only:write_corot_trial
 use work_diagnostic_435,only:update_work
 use geometry_provider_421,only:get_element_geometry
 implicit none
 character(80),intent(in)::cmname
 integer,intent(in)::ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,jstep(4),kinc
 real(real64),intent(inout)::stress(ntens),statev(nstatv),ddsdde(ntens,ntens),sse,spd,scd,pnewdt
 real(real64),intent(inout)::rpl,ddsddt(ntens),drplde(ntens),drpldt
 real(real64),intent(in)::stran(ntens),dstran(ntens),time(2),dtime,temp,dtemp,predef(*),dpred(*)
 real(real64),intent(in)::props(nprops),coords(3),drot(3,3),celent,dfgrd0(3,3),dfgrd1(3,3)
 integer::rc,log_error,diagnostics(8)
 real(real64)::before(105),after(105),oldwork(3),newwork(3)
 integer::wrc
 integer::i
 real(real64)::initial_g,initial_lambda
 real(real64)::j0,j1,fullprops(50),vv(59),extraold(10)
 interface
  subroutine xit()
  end subroutine
 end interface
 if(trim(cmname)/='CDPM2_VISCOUS_1049'.or.ndi/=3.or.nshr/=3.or.ntens/=6.or. &
    nstatv/=62.or.nprops/=24.or.jstep(2)/=1.or.jstep(3)/=1.or.jstep(4)/=0)then
  call xit();return
 end if
 if(.not.all(ieee_is_finite(dfgrd0)).or..not.all(ieee_is_finite(dfgrd1)))then
  call xit();return
 end if
 j0=det3(dfgrd0);j1=det3(dfgrd1)
 if(.not.ieee_is_finite(j0).or..not.ieee_is_finite(j1))then
  call xit();return
 end if
 if(j0<=0.or.j1<=0)then
  call xit();return
 end if
 ! Abaqus can ask for initial stiffness with KINC=0 before UEXTERNALDB LOP1
 ! makes FIL geometry available. Only the pristine first-step query is elastic;
 ! history/restart/real increments must still use the unchanged material path.
 if(kinc==0.and.jstep(1)==1.and.all(time==0.0_real64).and. &
    all(stran==0.0_real64).and.all(dstran==0.0_real64).and. &
    all(stress==0.0_real64).and.all(statev==0.0_real64).and. &
    j0==1.0_real64.and.j1==1.0_real64)then
  if(.not.all(ieee_is_finite(props(1:2))))then
   call xit();return
  endif
  if(props(1)<=0.0_real64.or.props(2)<=-1.0_real64.or.props(2)>=0.5_real64)then
   call xit();return
  endif
  initial_g=props(1)/(2.0_real64*(1.0_real64+props(2)))
  initial_lambda=props(1)*props(2)/((1.0_real64+props(2))*(1.0_real64-2.0_real64*props(2)))
  ddsdde=0.0_real64
  ddsdde(1:3,1:3)=initial_lambda
  do i=1,3
   ddsdde(i,i)=initial_lambda+2.0_real64*initial_g
   ddsdde(i+3,i+3)=initial_g
  enddo
  rpl=0.0_real64;ddsddt=0.0_real64;drplde=0.0_real64;drpldt=0.0_real64
  ! Do not alter stress, STATEV, work/energy variables or PNEWDT.
  return
 endif
 fullprops(1:23)=props(1:23)
 call get_element_geometry(noel,fullprops(24:50),rc)
 if(rc/=0)then
  call xit();return
 end if
 oldwork=statev(50:52)
 if(.not.all(ieee_is_finite(oldwork)))then
  call xit();return
 end if
 before=[stress,statev,reshape(ddsdde,[36]),pnewdt]
 extraold=statev(53:62)
 vv=[statev(1:49),statev(53:62)]
 call viscous_payload(fullprops,props(24),dtime,stran,dstran,drot,j1,stress,vv,ddsdde,pnewdt,rc,diagnostics)
 if(rc==0)then
  statev(1:49)=vv(1:49);statev(53:62)=vv(50:59)
 endif
 if(rc==0)then
  newwork=oldwork
  call update_work(oldwork,before(1:6),stress,dstran,statev(2:7),statev(53:58),newwork,wrc)
  if(wrc/=0)then
   stress=before(1:6);statev=before(7:68);ddsdde=reshape(before(69:104),[6,6]);pnewdt=before(105);rc=2
  else
   statev(50:52)=newwork
  end if
 end if
 after=[stress,statev,reshape(ddsdde,[36]),pnewdt]
 call write_corot_trial(noel,npt,jstep,kinc,rc,time,dtime,props(24),stran,dstran,drot,dfgrd0,dfgrd1,before,after,log_error,diagnostics)
 if(log_error/=0)then
  call xit();return
 end if
 if(rc==2)then
  call xit();return
 end if
 if(rc/=0)return
 rpl=0;ddsddt=0;drplde=0;drpldt=0
contains
 real(real64) function det3(a)
 real(real64),intent(in)::a(3,3)
 det3=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
     -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
     +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
 end function
end subroutine

