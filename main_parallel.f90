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

#if defined(PARALLEL)
  real, dimension(:), allocatable :: h[:], h_global[:]
  real, codimension[:], dimension(:), allocatable :: u
  
  !Parallel variables
  integer :: tile_size, ix, ex, right ,left
  integer, dimension(:), allocatable :: ip_core, im_core
#elif defined(SERIAL)
  real, dimension(:), allocatable :: h
  real, dimension(:), allocatable :: u
#endif
  integer :: i , it, outunit
  character(:), allocatable :: filename

#ifdef PARALLEL
  tile_size = grid_size/num_images()
  ix = tile_size*(this_image()-1) + 1
  ex = tile_size*this_image()
  
  !Allocate memory
  allocate(h_global(grid_size)[*])
  allocate(h(0:tile_size+1)[*])
  allocate(u(0:tile_size+1)[*])

  !PBC
  allocate(ip_core(num_images()))
  allocate(im_core(num_images()))

  do i = 1, num_images()
     ip_core(i) = i + 1
     im_core(i) = i - 1
  end do
  ip_core(num_images()) = 1
  im_core(1) = num_images()

  right = ip_core(this_image())
  left  = im_core(this_image())

  
  !Initialize h
  h(1:tile_size) = [ (exp( -decay *(i-icenter)**2 ), i = ix, ex) ]
  h_global(ix:ex)[1] = h(1:tile_size)
  filename = "shallow_parallel.dat"
#elif defined(SERIAL)
  allocate(h(grid_size))
  allocate(u(grid_size))
  filename = "shallow_serial.dat"
  do i = 1, grid_size
     h(i) = exp( -decay *(i-icenter)**2 )
  end do
#endif
 
  !Intialize u
  u = 0.0
 
  ! write data
#if defined(PARALLEL)
  sync all
  if(this_image() == 1) then
     open(newunit = outunit, file = filename)
     call write_data(outunit,h_global)
  end if
#elif defined(SERIAL)
  open(newunit = outunit, file = filename)
  call write_data(outunit,h)
#endif
  
  !Solve in time
  time_loop : do it = 1, time_steps
#ifdef PARALLEL
     h(0)[right] = h(tile_size)
     h(tile_size+1)[left] = h(1)
     sync all
#endif
     
     u = u - ( u*diff(u) + g*diff(h) )*dt/dx

#ifdef PARALLEL
     sync all
     u(0)[right] = u(tile_size)
     u(tile_size+1)[left] = u(1)
     sync all
#endif
     
     h = h - diff(u*(hmean + h))*dt/dx
     
#ifdef PARALLEL
     h_global(ix:ex)[1] = h(1:tile_size)
     sync all
     if(this_image() == 1) then
        call write_data(outunit,h_global)
     end if
#elif defined(SERIAL)
     call write_data(outunit,h)
#endif
  end do time_loop

contains

  function diff(x)
    real, intent(in) :: x(:)
    real :: diff(size(x))

    integer :: n
    
    n = size(x)
#if defined(SERIAL)
    diff(1) = x(2) - x(n)
    diff(n) = x(1) - x(n-1)
#endif
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
