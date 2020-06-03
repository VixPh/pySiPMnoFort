pure subroutine signalgenfortran(s,t,h,tf,tr,sigpts,gvar)
  implicit none

  integer(4),intent(in)  :: t,sigpts
  real(4),intent(in)     :: h,gvar,tf,tr
  real(4),intent(out)    :: s(sigpts)
  integer(4)             :: i
!f2py intent(in) t,sigpts,h,gvar,tf,tr
!f2py intent(hide),depend(sigpts) s

  s = 0
  forall (i=t+1:min(sigpts,int(t+10*tf))) &
    s(i) = (exp((t-i+1)/tf)-exp((t-i+1)/tr))
  s = gvar*h*s

end subroutine signalgenfortran


pure subroutine rollfortran(vecto,vect,t,gvar,h,npt)
  implicit none

  integer(4),intent(in)     :: t,npt
  real(4),intent(in)        :: gvar,h,vect(npt)
  real(4),intent(out)       :: vecto(npt)

!f2py intent(in) npt,t,gvar,h
!f2py intent(in),depend(npt) vect
!f2py intent(hide),depend(npt) vecto

  vecto = cshift(vect,-t,dim=1)
  vecto = gvar*h*vecto
  vecto(1:t) = 0
end subroutine rollfortran


subroutine randintfortran(randout,ncells,nsig)
  implicit none

  integer(4),intent(in)     :: ncells,nsig
  integer(4),intent(out)    :: randout(nsig)
  real(4)                   :: temp(nsig)
  integer(4)                :: i
!f2py intent(in) ncells,nsig
!f2py intent(hide),depend(nsig) randout

  call random_number(temp)
  randout = nint(temp * ncells)

end subroutine randintfortran


subroutine randnfortran(vectout,mu,sigma,n)
  implicit none

  real(4),intent(in)    :: mu,sigma
  integer(4),intent(in) :: n
  real(4),intent(out)   :: vectout(n)
  real(4)               :: temp(n)
  real (4), parameter   :: pi = 3.141592653589793E+00
!f2py intent(in) mu,sigma,n
!f2py intent(hide),depend(n) vectout

  call random_number(vectout)
  call random_number(temp)

  where ( vectout < tiny(vectout) ) vectout = tiny(vectout)

  vectout = sqrt( - 2.0E+00 * log(vectout))*cos( 2.0E+00*pi*temp)*sigma+mu

end subroutine randnfortran


subroutine randpoissfortran(vectout,mu,n)
  implicit none

  integer(4),intent(in)     :: n
  real(4),intent(in)        :: mu
  integer(4),intent(out)    :: vectout(n)
  real(4)                   :: L,p(n),u
  integer(4)                :: i,j
!f2py intent(in) mu,n
!f2py intent(hide),depend(n) vectout

  L = exp(-mu)
  vectout = 0
  p = 1
  call random_number(u)

  do i = 1, n
    do while ( p(i) > L )
      vectout(i) = vectout(i) + 1
      call random_number(u)
      p(i) = p(i)*u
    end do
  end do

  vectout = vectout - 1

end subroutine


subroutine randexpfortran(vectout,mu,n)
  implicit none


  real(4),intent(in)      :: mu
  integer(4),intent(in)   :: n
  real(4),intent(out)     :: vectout(n)
!f2py intent(in) mu,nd
!f2py intent(hide),depend(n) vectout

  call random_number(vectout)
  vectout = log(1-vectout)*(-mu)

end subroutine randexpfortran


subroutine sortfortran(array)
  implicit none

  real(4), intent(inout) :: array(:)
  integer :: i,j,left,right,last
  real(4) :: temp,p,next

  last=size(array)

  p=0.5*(array(1)+array(last))
  if (array(1).gt.array(last)) then
     temp=array(last)
     array(last)=array(1)
     array(1)=temp
  endif

  left=1
  right=last
  temp=array(2)

  do i=2,last-1
     if (temp.lt.p) then
        do j=left,1,-1
           if (array(j).le.temp) exit
           array(j+1)=array(j)
        end do
        array(j+1)=temp
        temp=array(left+2)
        left=left+1
     else
        next=array(right-1)
        do j=right,last
           if (array(j).ge.temp) exit
           array(j-1)=array(j)
        end do
        array(j-1)=temp
        temp=next
        right=right-1
     endif
  end do

  end subroutine sortfortran
