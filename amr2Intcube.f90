program amr2map
  !-----------------------------------------------------------------------
  ! Ce programme calcule la carte de densite projetee pour les
  ! variables hydro d'une simulation RAMSES.
  ! Version F90 par R. Teyssier le 01/04/01.
  !-----------------------------------------------------------------------
  use ramses_info
  implicit none
  character(128)::type
  integer::domax=0, dotot=0, domw=0, read_rt=0
  integer::n, i, j, k, twotondim, nx, ny, nz, ilevel, iidim
  integer:: idim, jdim,kdim, ivar, lmax=0, ind, ipos, ngrida,delta_iv
  integer::icpu, ncpu_read
  real::boxlen
  real(KIND=8)::xcen=0.d0,ycen=0.d0,zcen=0.d0,rcut=0.d0

  integer::nx_sample=0, ny_sample=0, nz_sample=0,ngrid, nvarh, nvarRT !--joki
  integer::imin, imax, jmin, jmax, kmin, kmax
  integer::nx_full,ny_full,nz_full,lmin,nboundary,ngrid_current

  integer::ix,iy,iz,iix,iiy,iiz,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,ddx,dxline,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,xx,yy,dz_proj,tot_lum,all_lum
  real(KIND=8),dimension(:,:),allocatable::pos,xg
  real(KIND=8),dimension(:,:,:),allocatable::var
  real(KIND=4),dimension(:,:,:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::map,temp2
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=1)::response
  character(LEN=5)::ncharcpu
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_cell,do_tot,hydroOk=.true. !--joki
  logical::rtOk=.false. !--joki
  logical::is_nHmThreshold=.false.,is_nHlThreshold=.false.
  real(kind=8)::nH_max,nH_max_UU,nh_min,nh_min_UU
  character(LEN=1)::proj='z',proj_sign='+'
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  integer::tmp1,tmp2
  integer::deli2,delj2,delk2,ixc,iyc,izc,ix0,ix1,iy0,iy1,iz0,iz1,nix,niy,niz

  ! LOS velocity parameters:
  integer::i_los_vel ! Index in vars of line-of-sight velocity
  real(kind=8)::los_vel_kms,vel_x,vel_y,dl_x,dl_y,dl_z,du_x,du_y,du_z,dc,uu,ul

  logical :: do_doppler=.true., ave_pix=.true.
  real(kind=8) :: v_hubble(3), dist(3), x_cen(3)
  integer :: iiv


  type level
     integer::ilevel
     integer::ngrid
     real(KIND=8),dimension(:,:,:),pointer::map,rhoz,vlos,velx,vely,wdop
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  tot_lum = 0.d0
  all_lum = 0.d0

  call read_params
  call read_info(repository)                                  ! joki

  if (rcut .lt. 0.d0) then
    rcut = -rcut/unit_l
  endif
  if (rcut .gt. 0.d0) then
    xmin = xcen - rcut
    xmax = xcen + rcut

    ymin = ycen - rcut
    ymax = ycen + rcut

    zmin = zcen - rcut
    zmax = zcen + rcut
  endif  

  if (eps_star .gt. 0.) then
    t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star
  else if (t_star .lt. 0.d0) then
    stop "error. Either eps_star or t_star need to be set > 0"
  endif
  t0    = t_star*(1d9)
  nISM0  = n_star
  nCOM  = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*0.76/1.66d-24
  d0  = max(del_star,nISM0)
  write(6,*) "using a t0,d0 of ", t0/1.e9,d0

  call prepareEmissivity(type)

  x_cen(1) = (xmax+xmin)/2.d0
  x_cen(2) = (ymax+ymin)/2.d0
  x_cen(3) = (zmax+zmin)/2.d0

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     hydroOk=.false.
     !stop
  endif
  if(read_rt==1) rtOK=.true.
  nomfich=TRIM(repository)//'/rt_'//TRIM(nchar)//'.out00001'  !--joki
  inquire(file=nomfich, exist=ok) ! verify input file         !--joki
  if ( .not. ok ) then                                        !--joki
     print *,TRIM(nomfich)//' not found.'                     !--joki
     rtOk=.false.                                             !--joki
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
  endif

  ! ------------
  ! Cooling
  ! ------------
  call read_cool(cool_format,repository,nchar,cool_n,cool_T,cool_spec)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     i = int( log(0.5d0**0.333333d0/aexp)/log(2.d0) + 1.d0 )
     lmax=nlevelmax-i
     write(*,*) "Setting lmax = nlevelmax = ", lmax
  endif
  if ( 2.d0**lmax * max(xmax-xmin,ymax-ymin,zmax-zmin) .gt. 1000 ) then
    write(6,*) "I am concerned this is too big a cube: Getting a max number of ", &
       (2.d0**lmax * max(xmax-xmin,ymax-ymin,zmax-zmin))**3, "cells"
    write(6,*) "This is likely too much for memory? Proceed (enter 'Y')"
    read(5,*) response
    if (response .ne. 'Y') then
      stop "Ending"
    endif
  endif

  write(*,*)'time=',t
  write(*,*)'Working resolution =',2**nlevelmax
  do_tot=.true.
  zzmax=1.0
  zzmin=0.0
  if (proj=='x')then
     idim=2
     jdim=3
     kdim=1
     xxmin=ymin ; xxmax=ymax
     yymin=zmin ; yymax=zmax
     zzmin=xmin ; zzmax=xmax
     i_los_vel=2
  else if (proj=='y') then
     idim=3
     jdim=1
     kdim=2
     xxmin=zmin ; xxmax=zmax
     yymin=xmin ; yymax=xmax
     zzmin=ymin ; zzmax=ymax
     i_los_vel=3
  else
     idim=1
     jdim=2
     kdim=3
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
     zzmin=zmin ; zzmax=zmax
     i_los_vel=4
  end if
  dz_proj = 1d0  ! For a projection, e.g. column density

  if(TRIM(ordering).eq.'hilbert')then
     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,nlevelmax
        dx=0.5d0**ilevel
        if(dx.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax

     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do

     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do

     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1                                ! joki
     kmax=int(zzmax*dble(nz_full))+1                                ! joki
     allocate(grid(ilevel)%map(kmin:kmax,jmin:jmax,imin:imax))
     allocate(grid(ilevel)%rhoz(kmin:kmax,jmin:jmax,imin:imax))
     allocate(grid(ilevel)%vlos(kmin:kmax,jmin:jmax,imin:imax))
     allocate(grid(ilevel)%velx(kmin:kmax,jmin:jmax,imin:imax))
     allocate(grid(ilevel)%vely(kmin:kmax,jmin:jmax,imin:imax))
     allocate(grid(ilevel)%wdop(kmin:kmax,jmin:jmax,imin:imax))
     grid(ilevel)%map(:,:,:)=0.0
     grid(ilevel)%rhoz(:,:,:)=0.0
     grid(ilevel)%vlos(:,:,:)=0.0
     grid(ilevel)%velx(:,:,:)=0.0
     grid(ilevel)%vely(:,:,:)=0.0
     grid(ilevel)%wdop(:,:,:)=0.0
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin                                         ! joki
     grid(ilevel)%kmax=kmax                                         ! joki
  end do

  do ilevel=lmax+1,nlevelmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1                                ! joki
     kmax=int(zzmax*dble(nz_full))+1                                ! joki
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin                                         ! joki
     grid(ilevel)%kmax=kmax                                         ! joki
  enddo

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------
  ! Loop over processor files
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     do i=1,21
        read(10)
     end do
     ! Read grid numbers
     read(10)ngridlevel
     ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
     read(10)
     if(nboundary>0)then
        do i=1,2
           read(10)
        end do
        read(10)ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
     endif
     read(10)
     read(10)
     if(TRIM(ordering).eq.'bisection')then
        do i=1,5
           read(10)
        end do
     else
        read(10)
     endif
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nvarh=0
     if (hydroOk) then
        nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
        open(unit=11,file=nomfich,status='old',form='unformatted')
        read(11)
        read(11)nvarh
        read(11)
        read(11)
        read(11)
        read(11)
     endif

     nvarRT=0
     if (rtOk) then
        ! Open RT file and skip header
        nomfich=TRIM(repository)//'/rt_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
        open(unit=12,file=nomfich,status='old',form='unformatted')
        read(12)
        read(12)nvarRT
        read(12)
        read(12)
        read(12)
        read(12)
     endif

     ! Loop over levels
     do ilevel=1,nlevelmax

        ! Geometry
        dx=0.5**ilevel
        dxline=1.d0
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg(1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           if(hydroOk) allocate(var(1:ngrida,1:twotondim,1:nvarh+nvarRT))
           allocate(pos(1:ngrida,1:ndim))
           allocate(map(1:ngrida))
           allocate(temp2(1:ngrida))
           allocate(spec(1:ngrida,6))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

           ! Read AMR data
           if(ngridfile(j,ilevel)>0)then
              read(10) ! Skip grid index
              read(10) ! Skip next index
              read(10) ! Skip prev index
              ! Read grid center
              do iidim=1,ndim
                 if(j.eq.icpu)then
                    read(10)xg(:,iidim)
                 else
                    read(10)
                 endif
              end do
              read(10) ! Skip father index
              do ind=1,2*ndim
                 read(10) ! Skip nbor index
              end do
              ! Read son index
              do ind=1,twotondim
                 if(j.eq.icpu)then
                    read(10)son(:,ind)
                 else
                    read(10)
                 end if
              end do
              ! Skip cpu map
              do ind=1,twotondim
                 read(10)
              end do
              ! Skip refinement map
              do ind=1,twotondim
                 read(10)
              end do
           endif

           ! Read HYDRO data
           if (hydroOk) then
              read(11)
              read(11)tmp1
              if(rtOk)read(12) !--joki
              if(rtOk)read(12)tmp2 !--joki
              if(ngridfile(j,ilevel)>0)then
                 ! Read hydro variables
                 do ind=1,twotondim
                    do ivar=1,nvarh
                       if(j.eq.icpu)then
                          read(11)var(:,ind,ivar)
                       else
                          read(11)
                       end if
                    end do
                    do ivar=1,nvarRT !--begin joki
                       if(j.eq.icpu)then
                          read(12)var(:,ind,nvarh+ivar)
                       else
                          read(12)
                       end if
                    end do          !--end joki
                 end do
              end if
           endif
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

              ! Compute cell center
              do i=1,ngrida
                 pos(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                 pos(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                 pos(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              end do
              ! Check if cell is refined
              do i=1,ngrida
                 ref(i)=son(i,ind)>0.and.ilevel<nlevelmax
              end do

              ! Extract variable
              if(index(type, 'cpu') .eq. 1) then !------------------------
                 map = icpu
              else if(index(type, 'ref') .eq. 1) then !-------------------
                 map = ilevel
              else
                 call emissivity(var, map, ind, dx, ngrida, temp2)
              endif

              if(is_nHmThreshold) then ! Set to zero everything above nH-threshold
                 nH_max_UU=nH_max/unit_nH
                 where (var(:,ind,1) .gt. nH_max_UU) map=0.d0
              endif
              if(is_nHlThreshold) then ! Set to zero everything below nH-threshold
                 nH_min_UU=nH_min/unit_nH
                 where (var(:,ind,1) .lt. nH_min_UU) map=0.d0
              endif

              all_lum = all_lum + sum(map(:))*(dx*unit_l/3.08e21)**2

              ! Store data map
              do i=1,ngrida
                 ok_cell= .not.ref(i)
                 if(ok_cell)then
                    ix=int(pos(i,idim)*dble(nx_full))+1
                    iy=int(pos(i,jdim)*dble(ny_full))+1
                    iz=int(pos(i,kdim)*dble(nz_full))+1
                    ! Lower weight for edges:
                    weight=(min(pos(i,kdim)+dx/2.,zzmax)   &
                         -max(pos(i,kdim)-dx/2.,zzmin))/dx
                    weight=min(1.0d0,max(weight,0.0d0))
                    if(    ix>=grid(ilevel)%imin.and.&
                         & iy>=grid(ilevel)%jmin.and.&
                         & iz>=grid(ilevel)%kmin.and.&
                         & ix<=grid(ilevel)%imax.and.&
                         & iy<=grid(ilevel)%jmax.and.&
                         & iz<=grid(ilevel)%kmax)                     then
                       ! Find correct velocity bin:
                       dist(:) = (pos(i,:) - x_cen(:))
                       v_hubble(:) = dist(:) * hubble
                       los_vel_kms = (var(i,ind,i_los_vel) + v_hubble(kdim)) * unit_v / 1d5
                       vel_x = (var(i,ind,idim+1) + v_hubble(idim)) * unit_v / 1d5
                       vel_y = (var(i,ind,jdim+1) + v_hubble(jdim)) * unit_v / 1d5

                       if (ilevel .le. lmax) then
                         grid(ilevel)%map(iz,iy,ix) = &
                               grid(ilevel)%map(iz,iy,ix) &
                               + map(i) * weight

                         grid(ilevel)%rhoz(iz,iy,ix) = &
                               grid(ilevel)%rhoz(iz,iy,ix) &
                               + var(i,ind,1) * var(i,ind,iPT+1)*dx * weight &
                               * unit_d*scale*unit_l

                         grid(ilevel)%vlos(iz,iy,ix) = &
                               grid(ilevel)%vlos(iz,iy,ix) &
                               + los_vel_kms * map(i) * weight

                         grid(ilevel)%velx(iz,iy,ix) = &
                               grid(ilevel)%velx(iz,iy,ix) &
                               + vel_x * map(i) * weight

                         grid(ilevel)%vely(iz,iy,ix) = &
                               grid(ilevel)%vely(iz,iy,ix) &
                               + vel_y * map(i) * weight

                         grid(ilevel)%wdop(iz,iy,ix) = &
                               grid(ilevel)%wdop(iz,iy,ix) &
                               + sqrt(2.d0*kb*temp2(i)/mH)/1.d5 * map(i) * weight
                         tot_lum = tot_lum + map(i)*(dx*unit_l/3.08e21)**2
                       else
                         ndom = 2**lmax
                         iix=int(pos(i,idim)*dble(ndom))+1
                         iiy=int(pos(i,jdim)*dble(ndom))+1
                         iiz=int(pos(i,kdim)*dble(ndom))+1

                         grid(lmax)%map(iiz,iiy,iix) = &
                               grid(lmax)%map(iiz,iiy,iix) &
                               + map(i) * weight / (2.d0**(ilevel - lmax))**2

                         grid(lmax)%rhoz(iiz,iiy,iix) = &
                               grid(lmax)%rhoz(iiz,iiy,iix) &
                               + var(i,ind,1) * var(i,ind,iPT+1)*dx * weight  &
                               * unit_d*scale*unit_l / (2.d0**(ilevel - lmax))**2

                         grid(lmax)%vlos(iiz,iiy,iix) = &
                               grid(lmax)%vlos(iiz,iiy,iix) &
                               + los_vel_kms * map(i) * weight / (2.d0**(ilevel - lmax))**2

                         grid(lmax)%velx(iiz,iiy,iix) = &
                               grid(lmax)%velx(iiz,iiy,iix) &
                               + vel_x * map(i) * weight / (2.d0**(ilevel - lmax))**2

                         grid(lmax)%vely(iiz,iiy,iix) = &
                               grid(lmax)%vely(iiz,iiy,iix) &
                               + vel_y * map(i) * weight / (2.d0**(ilevel - lmax))**2

                         ! Note, 2kT/m is most probably speed. 
                         ! RMS speed is 3kT/m, which you should divide / 3 to get LOS. 
                         grid(lmax)%wdop(iiz,iiy,iix) = &
                               grid(lmax)%wdop(iiz,iiy,iix) &
                               + sqrt(2.d0*kb*temp2(i)/mH)/1.d5 * map(i) * weight / (2.d0**(ilevel - lmax))**2

                         tot_lum = tot_lum + map(i)*(dx*unit_l/3.08e21)**2
            
                       endif
                    endif
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,ref,map,pos)
           deallocate(temp2,spec)
           if(hydroOk) deallocate(var)
        end if

     end do
     ! End loop over levels

     close(10)
     if(hydroOk) close(11)
     if(rtOk) close(12)

  end do
  ! End loop over cpu

  write(6,*) "Found a total SFR of (Msun/yr)", tot_sfr
  write(6,*) "Used a total SFR of (Msun/yr)", tot_sfr2
  write(6,*) "Using a total luminosity of (erg/s)", tot_lum
  write(6,*) "Equivalent to a SFR of (Msun/yr)", tot_lum/1.26d41
  write(6,*) "Fed a total luminosity of (erg/s)", all_lum

  nx_full=2**lmax
  ny_full=2**lmax
  nz_full=2**lmax
  imin=int(xxmin*dble(nx_full))+1
  imax=int(xxmax*dble(nx_full))
  jmin=int(yymin*dble(ny_full))+1
  jmax=int(yymax*dble(ny_full))
  kmin=int(zzmin*dble(nz_full))+1
  kmax=int(zzmax*dble(nz_full))

  ! BEGIN JOKI: Readjust the range to match the cell sides so ------------
  !             that it is completely accurate:
  xxmin = dble(imin-1) / dble(nx_full)
  xxmax = dble(imax)   / dble(nx_full)
  yymin = dble(jmin-1) / dble(ny_full)
  yymax = dble(jmax)   / dble(ny_full)
  zzmin = dble(kmin-1) / dble(nz_full)
  zzmax = dble(kmax)   / dble(nz_full)
  ! END JOKI -------------------------------------------------------------

  do ix=imin,imax
     xmin=((ix-0.5)/2**lmax)
     do iy=jmin,jmax
        ymin=((iy-0.5)/2**lmax)
        do iz=kmin,kmax
          zmin=((iz-0.5)/2**lmax)
          do ilevel=1,lmax-1
            ndom=2**ilevel
            i=int(xmin*ndom)+1
            j=int(ymin*ndom)+1
            k=int(zmin*ndom)+1
            grid(lmax)%map(iz,iy,ix)= grid(lmax)%map(iz,iy,ix)               &
                + grid(ilevel)%map(k,j,i)/2.d0**(lmax - ilevel) ! They are column densities, so need to split them in half 
            grid(lmax)%rhoz(iz,iy,ix)= grid(lmax)%rhoz(iz,iy,ix)               &
                + grid(ilevel)%rhoz(k,j,i)/2.d0**(lmax - ilevel) ! They are column densities, so need to split them in half
            grid(lmax)%vlos(iz,iy,ix)= grid(lmax)%vlos(iz,iy,ix)               &
                + grid(ilevel)%vlos(k,j,i)
            grid(lmax)%velx(iz,iy,ix)= grid(lmax)%velx(iz,iy,ix)               &
                + grid(ilevel)%velx(k,j,i)
            grid(lmax)%vely(iz,iy,ix)= grid(lmax)%vely(iz,iy,ix)               &
                + grid(ilevel)%vely(k,j,i)
            grid(lmax)%wdop(iz,iy,ix)= grid(lmax)%wdop(iz,iy,ix)               &
                + grid(ilevel)%wdop(k,j,i)
          end do
        end do
     end do
  end do

  do ilevel=1,lmax-1
    deallocate(grid(ilevel)%map)
    deallocate(grid(ilevel)%rhoz)
    deallocate(grid(ilevel)%vlos)
    deallocate(grid(ilevel)%velx)
    deallocate(grid(ilevel)%vely)
    deallocate(grid(ilevel)%wdop)
  enddo

  ! Output to binary file
  nomfich=TRIM(outfich)
  write(*,*)'Writing to cube file '//TRIM(nomfich)
  write(*,*)'Resolution = ',imax-imin+1,jmax-jmin+1,kmax-kmin+1

  open(unit=20,file=nomfich,form='unformatted')
  if(nx_sample==0)then
     ! Size of map (nx,ny) is given by maxlevel
     write(20) 5, kmax-kmin+1, jmax-jmin+1, imax-imin+1
     allocate(toto(5,kmax-kmin+1,jmax-jmin+1,imax-imin+1))
     toto(1,:,:,:) = log10(1.d-100 + grid(lmax)%map(kmin:kmax,jmin:jmax,imin:imax))
     toto(2,:,:,:) = log10(1.d-100 + grid(lmax)%rhoz(kmin:kmax,jmin:jmax,imin:imax))
     grid(lmax)%vlos(kmin:kmax,jmin:jmax,imin:imax) = &
       grid(lmax)%vlos(kmin:kmax,jmin:jmax,imin:imax)/(1.d-100 + grid(lmax)%map(kmin:kmax,jmin:jmax,imin:imax))
     grid(lmax)%velx(kmin:kmax,jmin:jmax,imin:imax) = &
       grid(lmax)%velx(kmin:kmax,jmin:jmax,imin:imax)/(1.d-100 + grid(lmax)%map(kmin:kmax,jmin:jmax,imin:imax))
     grid(lmax)%vely(kmin:kmax,jmin:jmax,imin:imax) = &
       grid(lmax)%vely(kmin:kmax,jmin:jmax,imin:imax)/(1.d-100 + grid(lmax)%map(kmin:kmax,jmin:jmax,imin:imax))

     toto(3,:,:,:) = grid(lmax)%vlos(kmin:kmax,jmin:jmax,imin:imax)
     toto(4,:,:,:) = 0.d0
     do ix=imin,imax
       do iy=jmin,jmax
         do iz=kmin,kmax
           dc = grid(lmax)%map(iz,iy,ix) 
           dl_x = grid(lmax)%map(iz,iy,max(ix-1,imin)) ; du_x = grid(lmax)%map(iz,iy,min(ix+1,imax))           
           dl_y = grid(lmax)%map(iz,max(iy-1,jmin),ix) ; du_y = grid(lmax)%map(iz,min(iy+1,jmax),ix)           
           dl_z = grid(lmax)%map(max(iz-1,kmin),iy,ix) ; du_z = grid(lmax)%map(min(iz+1,kmax),iy,ix)           

           ! x differences:
           if (ix .gt. imin) then
             if (ix .lt. imax) then
               ul = (dl_x * grid(lmax)%velx(iz,iy,ix-1) + dc* grid(lmax)%velx(iz,iy,ix))/(dl_x + dc)
               uu = (du_x * grid(lmax)%velx(iz,iy,ix+1) + dc* grid(lmax)%velx(iz,iy,ix))/(du_x + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = (uu-ul)**2
               ul = (dl_x * grid(lmax)%vely(iz,iy,ix-1) + dc* grid(lmax)%vely(iz,iy,ix))/(dl_x + dc)
               uu = (du_x * grid(lmax)%vely(iz,iy,ix+1) + dc* grid(lmax)%vely(iz,iy,ix))/(du_x + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_x * grid(lmax)%vlos(iz,iy,ix-1) + dc* grid(lmax)%vlos(iz,iy,ix))/(dl_x + dc)
               uu = (du_x * grid(lmax)%vlos(iz,iy,ix+1) + dc* grid(lmax)%vlos(iz,iy,ix))/(du_x + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             else
               ul = (dl_x * grid(lmax)%velx(iz,iy,ix-1) + dc* grid(lmax)%velx(iz,iy,ix))/(dl_x + dc)
               uu = (-0.5d0*dl_x * grid(lmax)%velx(iz,iy,ix-1) + 1.5d0*dc* grid(lmax)%velx(iz,iy,ix))/(-0.5d0*dl_x + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = (uu-ul)**2
               ul = (dl_x * grid(lmax)%vely(iz,iy,ix-1) + dc* grid(lmax)%vely(iz,iy,ix))/(dl_x + dc)
               uu = (-0.5d0*dl_x * grid(lmax)%vely(iz,iy,ix-1) + 1.5d0*dc* grid(lmax)%vely(iz,iy,ix))/(-0.5d0*dl_x + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_x * grid(lmax)%vlos(iz,iy,ix-1) + dc* grid(lmax)%vlos(iz,iy,ix))/(dl_x + dc)
               uu = (-0.5d0*dl_x * grid(lmax)%vlos(iz,iy,ix-1) + 1.5d0*dc* grid(lmax)%vlos(iz,iy,ix))/(-0.5d0*dl_x + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             endif
           else
             ul = (-0.5d0*du_x * grid(lmax)%velx(iz,iy,ix+1) + 1.5d0*dc* grid(lmax)%velx(iz,iy,ix))/(-0.5d0*du_x + 1.5d0*dc)
             uu = (du_x * grid(lmax)%velx(iz,iy,ix+1) + dc* grid(lmax)%velx(iz,iy,ix))/(du_x + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = (uu-ul)**2
             ul = (-0.5d0*du_x * grid(lmax)%vely(iz,iy,ix+1) + 1.5d0*dc* grid(lmax)%vely(iz,iy,ix))/(-0.5d0*du_x + 1.5d0*dc)
             uu = (du_x * grid(lmax)%vely(iz,iy,ix+1) + dc* grid(lmax)%vely(iz,iy,ix))/(du_x + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             ul = (-0.5d0*du_x * grid(lmax)%vlos(iz,iy,ix+1) + 1.5d0*dc* grid(lmax)%vlos(iz,iy,ix))/(-0.5d0*du_x + 1.5d0*dc)
             uu = (du_x * grid(lmax)%vlos(iz,iy,ix+1) + dc* grid(lmax)%vlos(iz,iy,ix))/(du_x + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
           endif

           ! y differences:
           if (iy .gt. jmin) then
             if (iy .lt. jmax) then
               ul = (dl_y * grid(lmax)%velx(iz,iy-1,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(dl_y + dc)
               uu = (du_y * grid(lmax)%velx(iz,iy+1,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(du_y + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_y * grid(lmax)%vely(iz,iy-1,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(dl_y + dc)
               uu = (du_y * grid(lmax)%vely(iz,iy+1,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(du_y + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_y * grid(lmax)%vlos(iz,iy-1,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(dl_y + dc)
               uu = (du_y * grid(lmax)%vlos(iz,iy+1,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(du_y + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             else
               ul = (dl_y * grid(lmax)%velx(iz,iy-1,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(dl_y + dc)
               uu = (-0.5d0*dl_y * grid(lmax)%velx(iz,iy-1,ix) + 1.5d0*dc* grid(lmax)%velx(iz,iy,ix))/(-0.5d0*dl_y + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_y * grid(lmax)%vely(iz,iy-1,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(dl_y + dc)
               uu = (-0.5d0*dl_y * grid(lmax)%vely(iz,iy-1,ix) + 1.5d0*dc* grid(lmax)%vely(iz,iy,ix))/(-0.5d0*dl_y + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_y * grid(lmax)%vlos(iz,iy-1,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(dl_y + dc)
               uu = (-0.5d0*dl_y * grid(lmax)%vlos(iz,iy-1,ix) + 1.5d0*dc* grid(lmax)%vlos(iz,iy,ix))/(-0.5d0*dl_y + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             endif
           else
             ul = (-0.5d0*du_y * grid(lmax)%velx(iz,iy+1,ix) + 1.5d0*dc* grid(lmax)%velx(iz,iy,ix))/(-0.5d0*du_y + 1.5d0*dc)
             uu = (du_y * grid(lmax)%velx(iz,iy+1,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(du_y + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             ul = (-0.5d0*du_y * grid(lmax)%vely(iz,iy+1,ix) + 1.5d0*dc* grid(lmax)%vely(iz,iy,ix))/(-0.5d0*du_y + 1.5d0*dc)
             uu = (du_y * grid(lmax)%vely(iz,iy+1,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(du_y + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             ul = (-0.5d0*du_y * grid(lmax)%vlos(iz,iy+1,ix) + 1.5d0*dc* grid(lmax)%vlos(iz,iy,ix))/(-0.5d0*du_y + 1.5d0*dc)
             uu = (du_y * grid(lmax)%vlos(iz,iy+1,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(du_y + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
           endif

           ! z differences:
           if (iz .gt. kmin) then
             if (iz .lt. kmax) then
               ul = (dl_z * grid(lmax)%velx(iz-1,iy,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(dl_z + dc)
               uu = (du_z * grid(lmax)%velx(iz+1,iy,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(du_z + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_z * grid(lmax)%vely(iz-1,iy,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(dl_z + dc)
               uu = (du_z * grid(lmax)%vely(iz+1,iy,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(du_z + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_z * grid(lmax)%vlos(iz-1,iy,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(dl_z + dc)
               uu = (du_z * grid(lmax)%vlos(iz+1,iy,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(du_z + dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             else
               ul = (dl_z * grid(lmax)%velx(iz-1,iy,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(dl_z + dc)
               uu = (-0.5d0*dl_z * grid(lmax)%velx(iz-1,iy,ix) + 1.5d0*dc* grid(lmax)%velx(iz,iy,ix))/(-0.5d0*dl_z + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_z * grid(lmax)%vely(iz-1,iy,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(dl_z + dc)
               uu = (-0.5d0*dl_z * grid(lmax)%vely(iz-1,iy,ix) + 1.5d0*dc* grid(lmax)%vely(iz,iy,ix))/(-0.5d0*dl_z + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
               ul = (dl_z * grid(lmax)%vlos(iz-1,iy,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(dl_z + dc)
               uu = (-0.5d0*dl_z * grid(lmax)%vlos(iz-1,iy,ix) + 1.5d0*dc* grid(lmax)%vlos(iz,iy,ix))/(-0.5d0*dl_z + 1.5d0*dc)
               toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             endif
           else
             ul = (-0.5d0*du_z * grid(lmax)%velx(iz+1,iy,ix) + 1.5d0*dc* grid(lmax)%velx(iz,iy,ix))/(-0.5d0*du_z + 1.5d0*dc)
             uu = (du_z * grid(lmax)%velx(iz+1,iy,ix) + dc* grid(lmax)%velx(iz,iy,ix))/(du_z + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             ul = (-0.5d0*du_z * grid(lmax)%vely(iz+1,iy,ix) + 1.5d0*dc* grid(lmax)%vely(iz,iy,ix))/(-0.5d0*du_z + 1.5d0*dc)
             uu = (du_z * grid(lmax)%vely(iz+1,iy,ix) + dc* grid(lmax)%vely(iz,iy,ix))/(du_z + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
             ul = (-0.5d0*du_z * grid(lmax)%vlos(iz+1,iy,ix) + 1.5d0*dc* grid(lmax)%vlos(iz,iy,ix))/(-0.5d0*du_z + 1.5d0*dc)
             uu = (du_z * grid(lmax)%vlos(iz+1,iy,ix) + dc* grid(lmax)%vlos(iz,iy,ix))/(du_z + dc)
             toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) = toto(4,iz-kmin+1,iy-jmin+1,ix-imin+1) + (uu-ul)**2
           endif

         enddo
       enddo
     enddo
     toto(4,:,:,:) = sqrt(toto(4,:,:,:)/3.d0)
     toto(5,:,:,:) = grid(lmax)%wdop(kmin:kmax,jmin:jmax,imin:imax)/(1.d-100 + grid(lmax)%map(kmin:kmax,jmin:jmax,imin:imax))
     write(20)toto
     write(*,*)'Min and Max Emis = ',minval(toto(1,:,:,:),toto(1,:,:,:)>-50.),maxval(toto(1,:,:,:))
     write(*,*)'Min and Max Z-col= ',minval(toto(2,:,:,:)),maxval(toto(2,:,:,:))
     write(*,*)'Min and Max Vlos = ',minval(toto(3,:,:,:)),maxval(toto(3,:,:,:))
     write(*,*)'Min and Max Siglos = ',minval(toto(4,:,:,:)),maxval(toto(4,:,:,:))
     write(*,*)'Min and Max DopW = ',minval(toto(5,:,:,:)),maxval(toto(5,:,:,:))
     write(20)xxmin,xxmax
     write(20)yymin,yymax
     write(20)zzmin,zzmax
     ! Write conversion factor of x,y,z-ranges to Mpc physical
     write(20)unit_l/3.08d24
  else
     ! Size of map (nx,ny) is use-defined
     if(ny_sample==0)ny_sample=nx_sample
     if(nz_sample==0)nz_sample=nx_sample
     write(20) 5, nz_sample+1, ny_sample+1, nx_sample+1
     allocate(toto(5,0:nz_sample,0:ny_sample,0:nx_sample))

     if (ave_pix) then
       deli2 = max(int( dble(imax-imin+1)/dble(nx_sample+1)/2 ),1)
       delj2 = max(int( dble(jmax-jmin+1)/dble(ny_sample+1)/2 ),1)
       delk2 = max(int( dble(kmax-kmin+1)/dble(nz_sample+1)/2 ),1)
     endif

     do i=0,nx_sample
        if (ave_pix) then
           ixc = int(dble(i) * dble(imax-imin+1)/dble(nx_sample+1))+imin+deli2
           ix0=min(max(ixc-deli2,imin),ixc)
           ix1=int(dble(i+1) * dble(imax-imin+1)/dble(nx_sample+1))+imin+deli2
           ix1=min(max(ix1-deli2-1,ixc),imax)
           nix = ix1 - ix0 + 1
        else
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
        endif
        do j=0,ny_sample
           if (ave_pix) then
             iyc = int(dble(j) * dble(jmax-jmin+1)/dble(ny_sample+1))+jmin+delj2
             iy0=min(max(iyc-delj2,jmin),iyc)
             iy1=int(dble(j+1) * dble(jmax-jmin+1)/dble(ny_sample+1))+jmin+delj2
             iy1=min(max(iy1-delj2-1,iyc),jmax)
             niy = iy1 - iy0 + 1
           else
             iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
             iy=min(iy,jmax)
           endif
           do k=0,nz_sample
             if (ave_pix) then
               izc = int(dble(k) * dble(kmax-kmin+1)/dble(nz_sample+1))+kmin+delk2
               iz0=min(max(izc-delk2,kmin),izc)
               iz1=int(dble(k+1) * dble(kmax-kmin+1)/dble(nz_sample+1))+kmin+delk2
               iz1=min(max(iz1-delk2-1,izc),kmax)
               niz = iz1 - iz0 + 1

               toto(1,i,j,k)=log10(1.d-100 + sum(grid(lmax)%map(iz0:iz1,iy0:iy1,ix0:ix1))/nix/niy/niz)
               toto(2,i,j,k)=log10(1.d-100 + sum(grid(lmax)%rhoz(iz0:iz1,iy0:iy1,ix0:ix1))/nix/niy/niz)
               toto(3,i,j,k)=sum(grid(lmax)%vlos(iz0:iz1,iy0:iy1,ix0:ix1))/(1.d-100 + sum(grid(lmax)%map(iz0:iz1,iy0:iy1,ix0:ix1)))
               toto(4,i,j,k)=sum( grid(lmax)%vlos(iz0:iz1,iy0:iy1,ix0:ix1)**2 + grid(lmax)%velx(iz0:iz1,iy0:iy1,ix0:ix1)**2 + & 
                                   grid(lmax)%vely(iz0:iz1,iy0:iy1,ix0:ix1)**2 ) / & 
                                   (1.d-100 + sum(grid(lmax)%map(iz0:iz1,iy0:iy1,ix0:ix1)))
               toto(4,i,j,k) = sqrt( (toto(4,i,j,k) - toto(3,i,j,k)**2)/3.d0 ) 
               toto(5,i,j,k)=sum(grid(lmax)%wdop(iz0:iz1,iy0:iy1,ix0:ix1))/(1.d-100 + sum(grid(lmax)%map(iz0:iz1,iy0:iy1,ix0:ix1)))
             else
               iz=int(dble(k)/dble(nz_sample)*dble(kmax-kmin+1))+kmin
               iz=min(iz,kmax)

               toto(1,k,j,i)=log10(1.d-100 + grid(lmax)%map(iz,iy,ix))
               toto(2,k,j,i)=log10(1.d-100 + grid(lmax)%rhoz(iz,iy,ix))
               toto(3,k,j,i)=grid(lmax)%vlos(iz,iy,ix)/(1.d-100 + grid(lmax)%map(iz,iy,ix))
               toto(4,k,j,i)=0.d0
               toto(5,k,j,i)=grid(lmax)%wdop(iz,iy,ix)/(1.d-100 + grid(lmax)%map(iz,iy,ix))
             endif

           end do
        end do
     end do
     write(20)toto
     write(20)xxmin,xxmax
     write(20)yymin,yymax
     write(20)zzmin,zzmax
     ! Write conversion factor of x,y,z-ranges to Mpc physical
     write(20)unit_l/3.08d24
  endif
  close(20)
  deallocate(toto)
  deallocate(grid(lmax)%map)
  deallocate(grid(lmax)%rhoz)
  deallocate(grid(lmax)%vlos)
  deallocate(grid(lmax)%velx)
  deallocate(grid(lmax)%vely)
  deallocate(grid(lmax)%wdop)


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
       print *, 'usage: amr2map   -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-dir axis] '
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-typ type] '
       print *, '                 [-avp average_pix_sample] '
       print *, '                 [-fil filetype] '
       print *, '                 [-dop do_doppler] '
       print *, 'ex: amr2map -inp output_00001 -out map.dat'// &
            &   ' -dir z -xmi 0.1 -xma 0.7 -lma 12'
       print *, ' '
       print *, ' type :-1 = cpu number'
       print *, ' type : 0 = ref. level (default)'
       print *, ' type : 1 = gas density'
       print *, ' type : 2 = X velocity'
       print *, ' type : 3 = Y velocity'
       print *, ' type : 4 = Z velocity'
       print *, ' type : 5 = gas pressure'
       print *, ' type : 6 = gas metallicity'
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
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-dir')
          proj = trim(arg)
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-xc')
          read (arg,*) xcen
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-yc')
          read (arg,*) ycen
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-zc')
          read (arg,*) zcen
       case ('-rc')
          read (arg,*) rcut
       case ('-lma')
          read (arg,*) lmax
       case ('-S2H')
          read (arg,*) SFRtoHA
       case ('-nst')
          read (arg,*) n_star
       case ('-tst')
          read (arg,*) t_star
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-typ')
          type = trim(arg)
       case ('-fil')
          read (arg,*) filetype
       case ('-rt')
          read (arg,*) read_rt
       case ('-nHm')
          is_nHmThreshold=.true.
          read (arg,*) nH_max
       case ('-nHl')
          is_nHlThreshold=.true.
          read (arg,*) nH_min
       case ('-avp') ! Maximum LOS velocity included
           read (arg,*) i_arg
           if (i_arg .eq. 1) then
             ave_pix = .true.
           else
             ave_pix = .false.
           endif
       case ('-dop') ! Maximum LOS velocity included
           read (arg,*) i_arg
           if (i_arg .eq. 1) then
             do_doppler = .true.
           else
             do_doppler = .false.
           endif
       case ('-rdc')
           read(arg,*) cool_format
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params

end program amr2map

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title

!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
