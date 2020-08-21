program amr2velcube_UsingInt

  implicit none

  real, allocatable, dimension(:,:,:,:) :: in_array
  real(kind=8), allocatable, dimension(:,:,:) :: out_array,trans_array

  character(len=128) :: intCube_file, out_file, ext_file
  character(len=1) :: look_sign
  integer :: nvar, nx, ny, nz, ix, iy, iz, iv, iiv, delta_iv
  integer :: front, back, delta, v_sign, use_v_cell=1
  real(kind=8) :: max_deltav=5.d0, deltav_inst=8.d0
  logical :: output_extinction=.false.

  ! Dust Values
  ! f_dust : fraction of gas metals in dust : ~30% Peeples et al. 2014
  ! dust_dens = 2.d0 - consistent more or less with carbonaceous chondrites, basalt etc.
  ! Qgeom = 1.5d0 - more or less true at small wavelengths. Don't currently have a non-constant implimentation of this.
  ! n_dust = -3.5 = slope of dust distribution. 3.5 = MRN distribution
  ! a_min = 0.005 microns = smallest grains
  ! a_max = 1 micron = largest grain
  real(kind=8) :: f_dust=0.3d0, screen_value=100.d0, dust_dens=3.d0, Qgeom=1.5d0, &
     n_dust = -3.5d0, a_min=5.d-7, a_max=1.d-4, oppac_factor

  real,parameter::pi=3.14159265359e0

  ! LOS velocity parameters:
  real(kind=8)::v_min_kms=-1d3, v_max_kms=1d3
  real(kind=8)::los_vel_kms, dvel, w_vel, val, wdop, Nmetl, vlos
  real(kind=4),dimension(:),allocatable::los_vel_kms_array
  logical::pix_clear
  real(kind=8)::tau
  real(kind=8),dimension(:,:),allocatable :: Av,Av2
  logical::do_doppler=.true.
  integer::nv=200 ! Number of equally spaced velocity bins (z-dimension)

  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,unit_l

  call read_params()

  dvel = (v_max_kms-v_min_kms)/nv
  if (use_v_cell .eq. 1) deltav_inst = dvel

  oppac_factor = 3.d0/4.d0 * f_dust * (4.d0+n_dust)/(3.d0+n_dust)* &
                      (a_max**(3.d0+n_dust)-a_min**(3.d0+n_dust)) / &
                      (a_max**(4.d0+n_dust)-a_min**(4.d0+n_dust))*Qgeom/dust_dens
  write(6,*) "Using an oppacity factor of ", oppac_factor

  open(42,file=intCube_file,status='old',form='unformatted')
  read(42) nvar,nz,ny,nx
  allocate(in_array(nvar,nz,ny,nx),out_array(nv,ny,nx))
  if (output_extinction) then
    allocate(Av(nx,ny),Av2(nx,ny))
    Av2 = 0.d0
  endif
  out_array = 0.d0


  read(42) in_array
  read(42)xxmin,xxmax
  read(42)yymin,yymax
  read(42)zzmin,zzmax
  read(42)unit_l
  close(42)

  write(6,*) "LOS Kinetic Doppler width is on (km/s)", minval(in_array(4,:,:,:)), maxval(in_array(4,:,:,:))
  write(6,*) "LOS Thermal Doppler width is on (km/s)", minval(in_array(5,:,:,:)), maxval(in_array(5,:,:,:))
  write(6,*) "Metal_dens on on ", minval(in_array(2,:,:,:)), maxval(in_array(2,:,:,:))

  if (look_sign .eq. '+') then
    front = 1
    back = nz
    delta = 1
    v_sign = 1.d0
  else
    front = nz
    back = 1
    delta = -1
    v_sign = -1.d0
  endif

  do ix=1,nx
    do iy=1,ny
      pix_clear=.true.
      tau = 0.d0
      do iz=front,back,delta
        if (pix_clear) then
          val  = 10.d0**in_array(1,iz,iy,ix)
          Nmetl= 10.d0**in_array(2,iz,iy,ix)
          vlos = v_sign*in_array(3,iz,iy,ix)
          wdop = in_array(4,iz,iy,ix)**2 + in_array(5,iz,iy,ix)**2
          if (.not. do_doppler) then
            wdop = 0.d0
          else
            wdop = sqrt(wdop + deltav_inst**2)
          endif
          iv = 1 + int( (vlos - v_min_kms) / dvel )
          !if(iz.lt.1 .or. iz.gt.nv) cycle ! Outside vel-range
          iv=max(iv,1)   ! Let outliers go into max and min
          iv=min(iv,nv)  ! velocities.
          ! Add to map
          delta_iv = int(max_deltav * wdop / dvel )
          w_vel = 0.d0
          do iiv=-delta_iv,delta_iv
            w_vel = w_vel + 1.d0/(1.d-100 + wdop) / sqrt(pi) * exp(-(iiv*dvel)**2/(1.d-100+wdop)**2)
          enddo
          do iiv=-delta_iv,delta_iv
            if ( (1 .le. iv+iiv) .and. (iv+iiv .le. nv) ) then
              out_array(iv+iiv,iy,ix) = out_array(iv+iiv,iy,ix) + val &
                      /(1.d-100 + wdop) / sqrt(pi) &
                      * exp(-(iiv*dvel)**2/(1.d-100+wdop)**2) / w_vel &
                      * exp(-tau)
              if (output_extinction) Av2(ix,iy) = Av2(ix,iy) + val &
                      /(1.d-100 + wdop) / sqrt(pi) &
                      * exp(-(iiv*dvel)**2/(1.d-100+wdop)**2) / w_vel
            endif
          enddo
          tau = tau + Nmetl*oppac_factor
          if (tau .gt. screen_value) pix_clear = .false.
          if (output_extinction) Av(ix,iy) = 1.086*tau 
        endif
      enddo
    enddo
  enddo
  deallocate(in_array)
  allocate(trans_array(nx,ny,nv))
  if (output_extinction) then
    do iy=1,ny
      do ix=1,nx
        if (Av2(ix,iy) .gt. 0.d0) then
          Av2(ix,iy) = -1.086d0*log(sum(out_array(:,iy,ix))/Av2(ix,iy))
        endif
      enddo
    enddo
  endif
  do ix=1,nx
    do iy=1,ny
      do iv=1,nv
        trans_array(ix,iy,iv) = out_array(iv,iy,ix)
      enddo
    enddo
  enddo
  deallocate(out_array)

  ! Generate velocities per bin:
  allocate(los_vel_kms_array(nv))
  do iv=1,nv
     los_vel_kms_array(iv) = &
          v_min_kms+(dble(iv-1))/dble(nv)*(v_max_kms-v_min_kms)
  end do
  ! Assign center velocities:
  los_vel_kms_array(:) = los_vel_kms_array(:) + 0.5*dvel

  open(42,file=out_file,status='unknown',form='unformatted')
  write(42) nx,ny,nv
  write(42) trans_array
  write(42)xxmin,xxmax
  write(42)yymin,yymax
  write(42)zzmin,zzmax
  ! Write conversion factor of x,y,z-ranges to Mpc physical
  write(42)unit_l
  write(42)los_vel_kms_array
  close(42)

  deallocate(los_vel_kms_array,trans_array)

  if (output_extinction) then
    open(42,file=ext_file,status='unknown')
    write(42,*) "x,y,Av,Av2"
    do ix=1,nx
      do iy=1,ny
        write(42,'(2(i3,","),es14.6,",",es14.6)') ix, iy, Av(ix,iy),Av2(ix,iy)
      enddo
    enddo
    close(42)
    deallocate(Av,Av2)
  endif
contains

  subroutine read_params

    implicit none

    integer       :: i,n,i_arg
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok

    n = iargc()
    if (n < 4) then
       print *, 'usage: amr2velcube  -inp  intermediate_cube'
       print *, '                    -out  output_file'
       print *, '                    [-sgn (+,-)]'
       print *, '                    [-nv n_velocity_bins] '
       print *, '                    [-vmi v_min_kms] '
       print *, '                    [-vma v_max_kms] '
       print *, '                    [-Nma screen_value]'
       print *, '                    [-fd f_dust]'
       print *, '                    [-dnd dust_dens]'
       print *, '                    [-Qgm Qgeom]'
       print *, '                    [-ami a_min]'
       print *, '                    [-ama a_max]'
       print *, '                    [-nd n_dust]'
       print *, '                    [-dvm max_deltav]'
       print *, '                    [-ext ext_file]'
       print *, 'ex: amr2velcube -inp IntCube_x_110.bin -out map.dat'// &
            &   ' -sgn + -nv 100 -vmi -2e3 -vma 2e3'
       print *, ' '
       stop
    end if

    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          intCube_file = trim(arg)
       case ('-out')
          out_file = trim(arg)
       case ('-sgn')
          look_sign = trim(arg)
       case ('-Nma')
          read (arg,*) screen_value
       case ('-fd')
          read (arg,*) f_dust
       case ('-dnd')
          read (arg,*) dust_dens
       case ('-ext')
          output_extinction = .true.
          ext_file = trim(arg)
       case ('-Qgm')
          read (arg,*) Qgeom
       case ('-nd')
          read (arg,*) n_dust
       case ('-ami')
          read (arg,*) a_min
       case ('-ama')
          read (arg,*) a_max
       case ('-nv') ! Number of velocity bins (3d dimension of cube)
          read (arg,*) nv
       case ('-vmi') ! Minimum LOS velocity included
          read (arg,*) v_min_kms
       case ('-vma') ! Maximum LOS velocity included
          read (arg,*) v_max_kms
       case ('-dvm') ! Delta-velocity LOS exponential width
           read (arg,*) max_deltav
       case ('-dvI') ! Instrumental velocity width
           read (arg,*) deltav_inst
       case ('-dop') ! Maximum LOS velocity included
           read (arg,*) i_arg
           if (i_arg .eq. 1) then
             do_doppler = .true.
           else
             do_doppler = .false.
           endif
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    return

  end subroutine read_params


end program amr2velcube_UsingInt
