! ============================================================================
! Program     : Convergence.f90
! Description : Extracts energy, force, stress, and computational time from OUTCAR
! Author      : Yida Yang
! Version     : tested
! ============================================================================

program convergence
  implicit none
  character(len=256) :: infile
 


  ! Define variables
  character(len=256) :: line  ! String to read lines from files
  character(len=10) :: ee,a,b,c,d,PP,PP1,PP2
  character (len=3)::conver(600)
  integer :: ios, step, NAT, j
  integer (kind=4) :: I
  real(8) :: energy, force_max, stress_max, ener(600), total_drift(600,3), cputime(600)
  real(8) :: de, force(600,3,3000), Position(3), force1(600,3000), forcemax(600), stress(600,9)
  real(8) :: stressmax(600), EDIFFG, VOLUME
  logical :: found
  if (command_argument_count() < 1) then
    print *, "Usage: compatition <OUTCAR file>"
    stop
  end if

  call get_command_argument(1, infile)
  open(unit=10, file=infile, status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open input file'
    stop
  end if
  ! Open OUTCAR file to read total energy
  open(unit=10, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  step = 0
  do
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit

    ! Extract nat
    if (index(line, 'NIONS') > 0) then
      read(line, *, iostat=ios) PP, PP, PP,PP,PP, nat, pp,pp,pp,pp,pp,nat
    endif
  enddo
  close(10)


  ! Open OUTCAR file to read total energy
  open(unit=10, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  step = 0
  do
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit

    ! Extract total energy
    if (index(line, 'energy  without entropy=') > 0) then
      step = step + 1
      read(line, *, iostat=ios) PP, PP, PP, ener(step)
    endif
  enddo
  close(10)

  ! Open OUTCAR file to read atomic forces
  open(unit=20, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  step = 0
  do
    read(20, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (index(line, 'POSITION') > 0) then
      step = step + 1
      read (20,*)  ! Skip one line

      DO I=1,NAT
        read(20,*) position(1:3), Force(STEP,1:3,i)  ! Read atomic positions and forces
      enddo

      read (20,*)  ! Skip one line
      read (20,*) a, b, total_drift(step,1:3)  ! Read total drift
    endif
  enddo
  close(20)

  ! Compute maximum forces
  DO I=1,STEP
    FORCE1 (I,1) = MAXVAL(Force(i,1,:))
    FORCE1 (I,2) = MAXVAL(Force(i,2,:))
    FORCE1 (I,3) = MAXVAL(Force(i,3,:))
  enddo

  DO I=1,STEP
    FORCEmax (I) = MAXVAL(Force1(i,:))
  enddo
  open(unit=20, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  do
    read(20, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (index(line, 'volume of cell :') > 0) then
      read(line, *, iostat=ios) PP, PP,pp,pp, VOLUME
    endif
  enddo
  close(20)
  ! Open OUTCAR file to read stress tensor
  open(unit=20, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  step = 0
  do
    read(20, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (index(line, 'stress matrix') > 0) then
      step = step + 1
      read (20,*) stress(step,1:3)
      read (20,*) stress(step,4:6)
      read (20,*) stress(step,7:9)
    endif
  enddo
  close(20)

  ! Compute maximum stress
  do i = 1, step
    stressmax(i) = maxval(stress(i,:))/VOLUME*160.21766208! GPA
  enddo

  ! Open OUTCAR file to read computational time
  open(unit=20, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  step = 0
  do
    read(20, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (index(line, 'LOOP+:') > 0) then
      step = step + 1
      read(line, *, iostat=ios) pp, PP1, PP2, PP, pp, pp, cputime(step)
    endif
  enddo
  close(20)

  ! Open OUTCAR file to read convergence criterion EDIFFG
  open(unit=20, file='OUTCAR', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Cannot open OUTCAR file'
    stop
  end if

  do
    read(20, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (index(line, 'EDIFFG') > 0) then
      read(line, *, iostat=ios) PP, PP, EDIFFG
    endif
  enddo
  close(20)
  do i=1,step
   if (Forcemax(i)<abs(EDIFFG))then
    conver(i) ="YES"
    endif
    if (Forcemax(i)>abs(EDIFFG))then
    conver(i)="NO"
   endif

  enddo

  ! Print results
  write (*,*) ' Step |   E (eV)  |  dE (eV)   | Fmax (eV/Å) | Smax (Gpa) | Time (s) | Average drift | Convergence |'
  write (*,*) '-------------------------------------------------------------------------------------------------------------'
  DO I=1,STEP
    if (i==1) then
      WRITE (*,"(2XI3,F13.6,5XA8,4F13.6,8XA4)") i, ener(i), "--------", FORCEmax(I), stressmax(i), cputime(i), &
          (total_drift(i,1) + total_drift(i,2) + total_drift(i,3)) / 3.0,conver(i)
    else
      WRITE (*,"(2XI3,6F13.6,8XA4)") i, ener(i), ener(i)-ener(i-1), FORCEmax(I), stressmax(i), cputime(i), &
          (total_drift(i,1) + total_drift(i,2) + total_drift(i,3)) / 3.0,conver (i)
    endif
  enddo

  WRITE (*,"(A7,F7.3)") "EDIFFG=", EDIFFG

end program convergence
