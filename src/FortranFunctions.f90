subroutine signalgenfortran(s,t,h,tf,tr,sigpts,gvar)
  implicit none
  integer(4),intent(in)  :: t,sigpts
  real(4),intent(in)     :: h,gvar,tf,tr
  real(4),intent(out)    :: s(sigpts)
  integer(4)             :: i
! f2py intent(in) t,sigpts,h,gvar,tf,tr
! f2py intent(out),depend(sigpts) s
  s = 0
  do i=t,min(sigpts,int(t+10*tf))
    s(i) = (exp((t-i)/tf)-exp((t-i)/tr))
  end do
  s = gvar*h*s
end subroutine

subroutine rollfortran(vecto,vect,t,gvar,h,npt)
  implicit none
  real(4),intent(in)        :: gvar,h
  integer(4),intent(in)     :: t,npt
  real(4),intent(in)        :: vect(npt)
  real(4),intent(out)       :: vecto(npt)
! f2py intent(in) npt,t,gvar,h
! f2py intent(in),depend(npt) vect
! f2py intent(out),depend(npt) vecto
  vecto = gvar*h*cshift(vect,-t,dim=1)
  vecto(1:t) = 0
end subroutine
