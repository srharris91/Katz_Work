!> \brief
!! This subroutine computes the diffusive flux.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 04/25/2010
!!
!! Additional notes:\par
!!  none
!!
!! Source code:
!!   \include rhsSYSTEMDisFlux.f90


SUBROUTINE cusp(a,lx,ly,xv,ql,qr,f)
      
IMPLICIT NONE

INTEGER,PARAMETER :: nq=4
REAL,INTENT(IN   ) :: a,lx,ly,xv
REAL,INTENT(IN   ),DIMENSION(nq ) :: ql,qr
REAL,INTENT(  OUT),DIMENSION(nq ) :: f
REAL :: rhl,rhr,de,ua,va,ha,qq,cc,qs,cs,qm,rc0,rc,rp,l,fdp,fr,fl,pr,pl, &
        efix,fdp5,hl,hr,gm1


efix = .5
gm1  = .4

qq  = .5*(ql(2)*ql(2)+ql(3)*ql(3))/ql(1)
pl  = gm1*(ql(4)-qq)
qq  = .5*(qr(2)*qr(2)+qr(3)*qr(3))/qr(1)
pr  = gm1*(qr(4)-qq)
hl  = ql(4)+pl
hr  = qr(4)+pr

rhl = SQRT(ql(1))
rhr = SQRT(qr(1))
de  = 1./(rhl+rhr)
rhl = 1./rhl
rhr = 1./rhr
ua  =(ql(2)*rhl+qr(2)*rhr)*de
va  =(ql(3)*rhl+qr(3)*rhr)*de
ha  =(hl   *rhl+hr   *rhr)*de
qq  = .5*(ua*ua+va*va)
cc  = gm1*(ha-qq)
qs  = lx*ua+ly*va-xv
cs  = SQRT(cc*(lx*lx+ly*ly))
qm  = qs/cs

rc0 = efix*cs
rc  = ABS(qs)
IF (rc < rc0) rc = .5*(rc0+rc*rc/rc0)
rc  = rc*a
rp  = SIGN(a,qs)
IF (ABS(qm) < 1.) THEN
   l  = ABS(qs)-cs
   rp = rp*MAX(0.,(ABS(qs)+l)/(ABS(qs)-l))
END IF
rc   = rc-rp*qs
qs   =(lx*ql(2)+ly*ql(3))/ql(1)-xv
fdp5 =-qs*pl
fl   = rc+rp*qs
qs   =(lx*qr(2)+ly*qr(3))/qr(1)-xv
fdp5 = fdp5+qs*pr
fdp5 = fdp5*rp
fr   = rc+rp*qs
fdp  = rp*(pr-pl)
f(1) = fr*qr(1)-fl*ql(1)
f(2) = fr*qr(2)-fl*ql(2)+lx*fdp
f(3) = fr*qr(3)-fl*ql(3)+ly*fdp
f(4) = fr*qr(4)-fl*ql(4)+xv*fdp+fdp5


END SUBROUTINE cusp
