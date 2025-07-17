program main

  implicit none

  integer, parameter :: grid_size = 100
  integer, parameter :: time_steps = 5000

  real, parameter :: dt = 0.02
  real, parameter :: dx = 1.0
  real, parameter :: g = 9.8

  real, parameter :: decay = 0.02
  real, parameter :: hmean = 10.0
  integer, parameter :: icenter = 25

  real, dimension(:), allocatable :: h
  real, dimension(:), allocatable :: u
  
  integer :: i , it, outunit

  !Allocate memory
  allocate(h(grid_size))
  allocate(u(grid_size))

  !Initialize h
  do i = 1, grid_size
     h(i) = exp( -decay *(i-icenter)**2 )
  end do

  !Intialize u
  u = 0.0

  ! write data
  open(newunit = outunit, file = 'heights.dat')
  call write_data(outunit,h)
  
  !Solve in time
  time_loop : do it = 1, time_steps
     u = u - ( u*diff(u) + g*diff(h) )*dt/dx
     h = h - diff(u*(hmean + h))*dt/dx
     call write_data(outunit,h)
  end do time_loop

contains

  function diff(x)
    real, intent(in) :: x(:)
    real :: diff(size(x))

    integer :: n

    n = size(x)

    diff(1) = x(2) - x(n)
    diff(n) = x(1) - x(n-1)

    diff(2:n-1) = x(3:n) - x(1:n-2)

    diff = 0.5*diff
    
  end function diff

  subroutine write_data(outunit,x)
    integer, intent(in) :: outunit
    real, intent(in) :: x(:)
    integer :: i

    do i = 1, size(x)
       write(outunit,*) i, x(i)
    end do
    write(outunit,'(2/)')
    
  end subroutine write_data
  
end program main
