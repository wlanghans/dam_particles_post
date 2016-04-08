module particle_data
  implicit none

  integer, parameter :: N_Scalar = 11!14


  character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qa","qv","qc","qi","qr","qs","qg","t","tabs","p"/)!,"ss","lh","b","conv_12","conv_13","conv_14","conv_x4","conv_x5","conv_x6","conv_y1"/)
  character(*), dimension(3), parameter :: Coord_var_name =(/"x","y","z"/)
  character(*), dimension(4), parameter ::   Vel_var_name =(/"u","v","w","vt"/)

  type particle_field
    integer :: g_id !global particle id
    integer :: Category
    integer :: natm
    logical :: Activity

    real(8) :: inactive_time
    real(8), dimension(3) :: Pos,Vel
    real(8) :: Vterm
    real(8), dimension(N_Scalar) :: scalar_var
  end type particle_field

  type(particle_field), allocatable, dimension(:), target :: part_data1, part_data2, part_data3
  type(particle_field), dimension(:), pointer :: part_new, part_new2, part_old, part_pt_tmp

  integer :: num_part

  real(8) :: dmass

  contains

  subroutine construct_particle
    implicit none
    integer :: i

    allocate( part_data1(num_part), part_data2(num_part) , part_data3(num_part))

    part_new => part_data1
    part_old => part_data2
    part_new2 => part_data3

    do i = 1, num_part
      part_new(i)%activity = .false.
      part_new2(i)%activity = .false.
      part_old(i)%activity = .false.
    enddo

  end subroutine construct_particle
    
  subroutine finalize_particle
    implicit none

    deallocate( part_data1, part_data2, part_data3)

  end subroutine finalize_particle


  subroutine compute_particle_mass(Nz,dz,dA)
    use mpi_info
    implicit none
    integer :: i,idx,N_t

    integer, dimension(:), allocatable :: N_j, local_mpi_int_buffer
    real(8), dimension(:), allocatable :: local_mass, local_mpi_real_buffer

    integer, intent(in) :: Nz
    real(8), intent(in) :: dz,dA

    allocate( N_j(Nz), local_mass(Nz) )
    allocate( local_mpi_int_buffer (Nz) )
    allocate( local_mpi_real_buffer(Nz) )

    N_j = 0; N_t = 0
    local_mass = 0.d0

    do i = 1, num_part
      idx = int( part_new(i)%Pos(3)/dz ) + 1
      if(idx < Nz+1) then
        N_j(idx) = N_j(idx) + 1
        local_mass(idx) = local_mass(idx) + part_new(i)%scalar_var(1)
      endif
    enddo

    call MPI_Allreduce(N_j       ,local_mpi_int_buffer ,Nz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    call MPI_Allreduce(local_mass,local_mpi_real_buffer,Nz,MPI_REAL8  ,MPI_SUM,MPI_COMM_WORLD,mpi_err)

    N_j        = local_mpi_int_buffer
    local_mass = local_mpi_real_buffer


    local_mass = dz*dA*local_mass/(dfloat(N_j)+1.d-13)
    N_t = sum(N_j)
    dmass = sum(local_mass)/dfloat(N_t)

    if(root) write(*,*) 'Mass of the particle:',dmass

    deallocate( N_j, local_mass )
    deallocate( local_mpi_int_buffer, local_mpi_real_buffer)

  end subroutine compute_particle_mass

  subroutine select_particles(step_in,N_sel,num_sel,rad_in,xc,yc,id_vec)
    use mpi_info
    implicit none

    integer,intent(in) :: num_sel,step_in
    integer,dimension(num_sel),intent(inout) :: id_vec
    integer,intent(inout) :: N_sel
    real,intent(in) :: rad_in,xc,yc

    integer,dimension(2) :: itmpo
    integer,parameter :: init_iseed=12433
    real(8),dimension(num_sel) :: rand
    integer,dimension(num_sel) :: rand_int
    integer,dimension(num_part) :: id_rad
    integer :: pcat
    real(8) :: px,py,pz,dist
    real(8),dimension(num_sel) :: x_vec,y_vec,z_vec,cat_vec, condrho_vec
    character(5) :: rank,step_char
    character(100) :: fname_out
    integer :: i, ierr, N_rad, id_in

    if (step_in.eq.0) then 

     N_rad=0
     do i = 1, num_part
      pz=part_new(i)%Pos(3)
      py=part_new(i)%Pos(2)
      px=part_new(i)%Pos(1)
      dist=sqrt((px-xc)**2.+(py-yc)**2.)

      !if (dist.gt.rad_in.or.(pz.gt.475..or.pz.lt.325.)) cycle
      if (dist.gt.rad_in.or.pz.gt.500.) cycle

      N_rad=N_rad+1
      id_rad(N_rad) = i
     end do


    ! initialize random number generator
    itmpo(1) = 1400; itmpo(2) = init_iseed+my_id
    call random_seed(put=itmpo)
    call random_number(rand)
    rand_int=INT(rand*float(N_rad))
    N_sel=0
    id_vec=0

     do i = 1, num_sel
    
      id_in = id_rad(rand_int(i))
    
      if (any(id_vec.eq.id_in)) cycle ! to avoid double counting 

      py=part_new(id_in)%Pos(2)
      px=part_new(id_in)%Pos(1)

      pcat = part_new(id_in)%Category
      pz=part_new(id_in)%Pos(3)

      N_sel = N_sel + 1
      x_vec(N_sel) = px
      y_vec(N_sel) = py
      z_vec(N_sel) = pz
      id_vec(N_sel) = id_in
      cat_vec(N_sel) = float(pcat)
      !condrho_vec(N_sel) = part_new(id_in)%scalar_var(1)*&
      !  (part_new(id_in)%scalar_var(4)+part_new(id_in)%scalar_var(5)+part_new(id_in)%scalar_var(6)+&
      !  part_new(id_in)%scalar_var(7)+part_new(id_in)%scalar_var(8))

     end do     

    else
    
     do i = 1, N_sel
       x_vec(i) = part_new(id_vec(i))%Pos(1)
       y_vec(i) = part_new(id_vec(i))%Pos(2)
       z_vec(i) = part_new(id_vec(i))%Pos(3)
       cat_vec(i) = float(part_new(id_vec(i))%Category)
      ! if (cat_vec(i).eq.8.) then
      !   condrho_vec(i) = -999.0
      ! else
      !   condrho_vec(i) = part_new(id_vec(i))%scalar_var(1)*&
      !    (part_new(id_vec(i))%scalar_var(4)+part_new(id_vec(i))%scalar_var(5)+part_new(id_vec(i))%scalar_var(6)+&
      !    part_new(id_vec(i))%scalar_var(7)+part_new(id_vec(i))%scalar_var(8))
      ! end if
        
     end do

    end if

    write(rank,'(i5)') my_id
    write(step_char,'(i5)') step_in
    fname_out ='./select_particles/animation_lt500/zlt500_radlt1000_data_'//trim(adjustl(rank))//'_'//trim(adjustl(step_char))//'.txt'
    !fname_out ='./select_particles/000_150/000_150_data_'//trim(adjustl(rank))//'_'//trim(adjustl(step_char))//'.txt'
    write(*,*) N_sel,' particles on rank ', my_id
    if (N_sel.gt.0) then
     open(unit=130,file=fname_out,status='replace',action=&
          'write', iostat=ierr)
     do i=1,N_sel
       write(130,'(f6.0,2x,f6.0,2x,f6.0,2x,f3.0)') x_vec(i),y_vec(i),z_vec(i),cat_vec(i)!,condrho_vec(i)
     end do
     close(unit=130)
    end if

  end subroutine select_particles 

   subroutine compute_particle_mass_subdomain_hori(Nz,Nx,Ny,dx,dy,dz,zl,zr,yl,yr,local_mass,N_j,icat)
    use mpi_info
    implicit none
    integer :: i,idx,idx_x, idx_y, pcat, N_t

    integer, dimension(:), allocatable :: local_mpi_int_buffer
    integer, intent(out),dimension(:), allocatable :: N_j
    real(8), dimension(:), allocatable :: local_mpi_real_buffer
    real(8), intent(out), dimension(:), allocatable :: local_mass
    real(8), intent(in) :: dx, dy

    integer, intent(in) :: Nz, Nx, Ny, icat
    real(8), intent(in) :: dz,zl,zr,yl,yr
    real(8) :: dA
    integer,dimension(2) :: inowater
    real(8),dimension(Nz*Ny,Nx) :: ind_xy
    real(8),dimension(Nx) :: sum_ind_xy

    allocate( N_j(Nx), local_mass(Nx) )
    allocate( local_mpi_int_buffer (Nx) )
    allocate( local_mpi_real_buffer(Nx) )

    inowater(1)=1
    inowater(2)=8

    N_j = 0; N_t = 0
    local_mass = 0.d0
    dA = (zr-zl) * (yr-yl)

    ind_xy=0.0

      do i = 1, num_part

      pcat = part_new(i)%Category

      if (icat.eq.pcat.or.(icat.eq.9.and..not.any(inowater.eq.pcat))) then

      idx   = int(part_new(i)%Pos(3) / dz)   + 1
      idx_x = int(part_new(i)%Pos(1) / dx)   + 1
      idx_y = int(part_new(i)%Pos(2) / dy)   + 1

      if(idx_x < Nx+1.and. (part_new(i)%Pos(3).gt.zl).and.(part_new(i)%Pos(3).le.zr) .and. &
            (part_new(i)%Pos(2).gt.yl).and.(part_new(i)%Pos(2).le.yr)) then
        N_j(idx_x) = N_j(idx_x) + 1
        if (icat.eq.9) then
        local_mass(idx_x) = local_mass(idx_x) + part_new(i)%scalar_var(1)* &
        (part_new(i)%scalar_var(3)+part_new(i)%scalar_var(4)+part_new(i)%scalar_var(5)+part_new(i)%scalar_var(6)+&
        part_new(i)%scalar_var(7)+part_new(i)%scalar_var(8))
        if (ind_xy(idx*idx_y,idx_x).ne.1.0) ind_xy(idx*idx_y,idx_x)=1.0
        else
        local_mass(idx_x) = local_mass(idx_x) + part_new(i)%scalar_var(1)*part_new(i)%scalar_var(pcat+1)
        if (ind_xy(idx*idx_y,idx_x).ne.1.0) ind_xy(idx*idx_y,idx_x)=1.0
        end if
      endif
      endif
    enddo

    do i=1,Nx
     sum_ind_xy(i)=sum(ind_xy(:,i))
    end do

    call MPI_Allreduce(N_j       ,local_mpi_int_buffer,Nx,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    call MPI_Allreduce(local_mass,local_mpi_real_buffer,Nx,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

    N_j        = local_mpi_int_buffer
    local_mass = local_mpi_real_buffer

    call MPI_Allreduce(sum_ind_xy,local_mpi_real_buffer,Nx,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    sum_ind_xy=local_mpi_real_buffer
    if (root) then
     do i=1,Nx
      write(*,*) i,'  ',sum_ind_xy(i)
     end do
    end if

       !mass per slice
    local_mass = dz*dx*dy*sum_ind_xy*local_mass
    if (root) then
     do i=1,Nx
      write(*,*) i,'  ',local_mass(i)
     end do
    end if
    do i=1,Nx
      local_mass(i) = local_mass(i)/(dfloat(N_j(i))+1.d-13)
    end do

    N_t = sum(N_j)
    dmass = sum(local_mass)/dfloat(N_t)

    if(root) write(*,*) 'Mass of the particle:',dmass


    deallocate( local_mpi_int_buffer, local_mpi_real_buffer)

  end subroutine compute_particle_mass_subdomain_hori



  subroutine compute_particle_mass_subdomain(Nz,Nx,Ny,dx,dy,dz,xl,xr,yl,yr,local_mass,N_j,icat)
    use mpi_info
    implicit none
    integer :: i,idx,idx_x, idx_y, pcat, N_t

    integer, dimension(:), allocatable :: local_mpi_int_buffer
    integer, intent(out),dimension(:), allocatable :: N_j
    real(8), dimension(:), allocatable :: local_mpi_real_buffer
    real(8), intent(out), dimension(:), allocatable :: local_mass
    real(8), intent(in) :: dx, dy

    integer, intent(in) :: Nz, Nx, Ny, icat
    real(8), intent(in) :: dz,xl,xr,yl,yr
    real(8) :: dA
    integer,dimension(2) :: inowater
    real(8),dimension(Nx*Ny,Nz) :: ind_xy
    real(8),dimension(Nz) :: sum_ind_xy
 
    allocate( N_j(Nz), local_mass(Nz) )
    allocate( local_mpi_int_buffer (Nz) )
    allocate( local_mpi_real_buffer(Nz) )

    inowater(1)=1
    inowater(2)=8

    N_j = 0; N_t = 0
    local_mass = 0.d0
    dA = (xr-xl) * (yr-yl)

    ind_xy=0.0

    do i = 1, num_part

      pcat = part_new(i)%Category

      if (icat.eq.pcat.or.(icat.eq.9.and..not.any(inowater.eq.pcat))) then

      idx   = int(part_new(i)%Pos(3) / dz)   + 1
      idx_x = int(part_new(i)%Pos(1) / dx)   + 1
      idx_y = int(part_new(i)%Pos(2) / dy)   + 1

      if(idx < Nz+1.and. (part_new(i)%Pos(1).gt.xl).and.(part_new(i)%Pos(1).le.xr) .and. &
            (part_new(i)%Pos(2).gt.yl).and.(part_new(i)%Pos(2).le.yr)) then
        N_j(idx) = N_j(idx) + 1
        if (icat.eq.9) then
        local_mass(idx) = local_mass(idx) + part_new(i)%scalar_var(1)* & 
        (part_new(i)%scalar_var(3)+part_new(i)%scalar_var(4)+part_new(i)%scalar_var(5)+part_new(i)%scalar_var(6)+&
        part_new(i)%scalar_var(7)+part_new(i)%scalar_var(8))
        if (ind_xy(idx_x*idx_y,idx).ne.1.0) ind_xy(idx_x*idx_y,idx)=1.0
        else
        local_mass(idx) = local_mass(idx) + part_new(i)%scalar_var(1)*part_new(i)%scalar_var(pcat+1)
        if (ind_xy(idx_x*idx_y,idx).ne.1.0) ind_xy(idx_x*idx_y,idx)=1.0
        end if
      endif
      endif
    enddo
  
    do i=1,Nz
     sum_ind_xy(i)=sum(ind_xy(:,i))
    end do
    
    call MPI_Allreduce(N_j       ,local_mpi_int_buffer ,Nz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    call MPI_Allreduce(local_mass,local_mpi_real_buffer,Nz,MPI_REAL8  ,MPI_SUM,MPI_COMM_WORLD,mpi_err)

    N_j        = local_mpi_int_buffer
    local_mass = local_mpi_real_buffer

    call MPI_Allreduce(sum_ind_xy,local_mpi_real_buffer,Nz,MPI_REAL8  ,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    sum_ind_xy=local_mpi_real_buffer
    if (root) then
     do i=1,Nz
      write(*,*) i,'  ',sum_ind_xy(i)
     end do
    end if

    !mass per layer
    local_mass = dz*dx*dy*sum_ind_xy*local_mass
    if (root) then
     do i=1,Nz
      write(*,*) i,'  ',local_mass(i)
     end do
    end if
    do i=1,Nz 
      local_mass(i) = local_mass(i)/(dfloat(N_j(i))+1.d-13)
    end do

    N_t = sum(N_j)
    dmass = sum(local_mass)/dfloat(N_t)

    if(root) write(*,*) 'Mass of the particle:',dmass


    deallocate( local_mpi_int_buffer, local_mpi_real_buffer)

  end subroutine compute_particle_mass_subdomain


  function get_activity(qc_cutoff, w_cutoff, w_in,q_in)
    implicit none

    real(8), intent(in) :: w_in, q_in,qc_cutoff, w_cutoff
    logical get_activity

    if(w_in > w_cutoff .and. q_in > qc_cutoff) then
      get_activity = .true.
    else
      get_activity = .false.
    endif
  end function get_activity

  function get_Thetae(q_in,T_in,P_in)
    implicit none

    !Define Thermodynamic Variables
    real(8), parameter :: rgasa =  287.04d0   ! Gas constant for dry air, J/kg/K
    real(8), parameter :: rgasv =  461.4d0    ! Gas constant for water vapor, J/kg/K
    real(8), parameter :: cva   =  719.d0     ! Heat capacity at constant volume for dry air, J/kg/K
    real(8), parameter :: cvv   = 1418.d0     ! Heat capacity at constant volume for water vapor, J/kg/K
    real(8), parameter :: cvl   = 4216.d0     ! Heat capacity at constant volume for liquid water, J/kg/K
    real(8), parameter :: cvs   = 2106.d0     ! Heat capacity at constant volume for solid water, J/kg/K
    real(8), parameter :: cpa   = cva + rgasa
    real(8), parameter :: cpv   = cvv + rgasv
    real(8), parameter :: E0v   = 2.374d6       ! Internal energy difference between vapor and liquid at the triple point, J/kg
    real(8), parameter :: ptrip    = 611.65d0   ! Triple point pressure, Pa
    real(8), parameter :: tabstrip = 273.16d0   ! Triple point temperature, K
    real(8), parameter :: pref     = 1.d5       ! Reference pressure, Pa
    real(8), parameter :: tabsref  = 3.d2       ! Reference temperature, K
    real(8), parameter :: s0v = E0v/tabstrip + rgasv
    real(8), parameter :: s0s = 3.337e5/tabstrip

    real(8), intent(in) :: q_in(3), T_in, P_in
    real(8) :: get_Thetae

    real(8) :: qa, qv, ql, qs, rv, rl, rs, pa, pv 

    qa = max(0.d0,q_in(1) )
    qv = max(0.d0,q_in(2) )
    ql = max(0.d0,q_in(3) )
    qs = 1.d0 - qv - ql - qa
 
    rv = qv / qa
    rl = ql / qa
    rs = qs / qa
    pa = rgasa/(rgasa+rv*rgasv)*P_in
    pv = rv*rgasv/(rgasa+rv*rgasv)*P_in

    get_Thetae = T_in * (pref/pa)**(rgasa/cpa) * &
                  (T_in/tabstrip)**((rv*cpv + rl*cvl + rs*cvs)/cpa) * &
                  (pv  /ptrip   )**(-rv*rgasv/cpa) * &
                  dexp( (rv*s0v - rs*s0s)/cpa)

  end function get_Thetae

end module particle_data
