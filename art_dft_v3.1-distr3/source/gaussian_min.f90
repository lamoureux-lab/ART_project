! In this routine, we interface between ART and GAUSSIAN
subroutine min_converge_gau(success)
  use defs

  integer :: i, ret

  logical :: success
  integer, parameter :: FGAUSS = 21
  real*8, parameter :: ZERO = 0.0d0
  real*8, dimension(natoms) :: xx, yy, zz
  character(len=20) :: GAUSS   = 'art2gaussian.inp'
  character(len=20) :: GAUSSFORCE = 'gaussian2art'
  character(len=70) :: line
  character(len=10) :: string_natoms
  logical :: read_coordinates_done,read_energy_done
  real(kind=8),dimension(3) :: boxl

  boxl = box * scala

  ! We first write to a file the format requested by GAUSS
  call system('python create_header.py -k opt')

  open(unit=FGAUSS,file=GAUSS,status='unknown',action='write', position = 'append', iostat=ierror)

  ! Prepares gaussian input file coordinates
  do i = 1,NATOMS
  write(FGAUSS,"(a, f14.8, f14.8, f14.8)") typat(i),  x(i), y(i), z(i)
  enddo
  write(FGAUSS,*) !Blank line... Gaussian expects one blank line at the end of the input file.
  close(FGAUSS)

  ! We now call Gaussian do to the minimization

  call system('python execute_gaussian.py')
  
  ! We must now read the energy and positions from gaussian's output file
  open(unit=FGAUSS,file=GAUSSFORCE,status='old',action='read',iostat=ierror)
  do 
    !Gets the Gaussian output coordinates
    read(FGAUSS,'(A40)') line
    if ( line  == "outcoor:" ) then
      do i = 1, NATOMS
        read (FGAUSS,*) x(i), y(i), z(i)
        write (*,*) "pos_min", x(i),y(i),z(i)
      end do
    endif

    if (line == 'energy:') then
        read(FGAUSS,*) total_energy
    endif
    if (line == 'Stationary:') exit
  end do

  close(FGAUSS)

  success = .true.
  return

end subroutine
