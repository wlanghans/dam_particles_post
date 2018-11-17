module particle_data
  implicit none

 ! integer, parameter :: N_Scalar = 15
  integer, parameter :: N_Scalar = 6

  real(8), parameter :: qc_cutoff = 1.d-5, w_cutoff = -100.d0

!  character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qa","qv","qc","qr","qi","t","tabs","p","b","Fdyn","Fbuoy","Fbuoyd","Fbuoyv","Fbuoyc"/)

  character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qv","qc","qi","t","tabs"/)!,"Fdyn","Fbuoy","Fbuoyd","Fbuoyv","Fbuoyc"/)
  character(*), dimension(3), parameter :: Coord_var_name =(/"x","y","z"/)
  character(*), dimension(3), parameter ::   Vel_var_name =(/"u","v","w"/)

  type particle_field
    integer :: g_id !global particle id
    logical :: Activity

    real(8) :: inactive_time
    real(8), dimension(3) :: Pos,Vel
    real(8), dimension(N_Scalar) :: scalar_var
  end type particle_field

  type(particle_field), allocatable, dimension(:), target :: part_data1, part_data2
  type(particle_field), dimension(:), pointer :: part_new, part_old

  integer :: num_part

  real(8) :: dmass

  contains

  subroutine construct_particle
    implicit none
    integer :: i

    allocate( part_data1(num_part), part_data2(num_part) )

    part_new => part_data1
    part_old => part_data2

    do i = 1, num_part
      part_new(i)%activity = .false.
      part_old(i)%activity = .false.
    enddo

  end subroutine construct_particle
    
  subroutine finalize_particle
    implicit none

    deallocate( part_data1, part_data2)

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


    local_mass = local_mass/(dfloat(N_j)+1.d-13)
    N_t = sum(N_j)
    !WL2013
    dmass = dz*dA/dfloat(N_t)*sum(local_mass)
    !dmass = dz*dA*sum(local_mass)

    if(root) write(*,*) 'Mass of the particle:',dmass

    deallocate( N_j, local_mass )
    deallocate( local_mpi_int_buffer, local_mpi_real_buffer)

  end subroutine compute_particle_mass

  subroutine compute_particle_mass_subdomain(Nz,dx,dy,dz,xl,xr,yl,yr,local_mass,N_j)
    use mpi_info
    implicit none
    integer :: i,idx,idx_x, idx_y

    integer, dimension(:), allocatable :: local_mpi_int_buffer
    integer, intent(out),dimension(:), allocatable :: N_j
    real(8), dimension(:), allocatable :: local_mpi_real_buffer
    real(8), intent(out), dimension(:), allocatable :: local_mass
    real(8), intent(in) :: dx, dy

    integer, intent(in) :: Nz
    real(8), intent(in) :: dz,xl,xr,yl,yr
    real(8) :: dA
 
    allocate( N_j(Nz), local_mass(Nz) )
    allocate( local_mpi_int_buffer (Nz) )
    allocate( local_mpi_real_buffer(Nz) )

    N_j = 0
    local_mass = 0.d0
    dA = (xr-xl) * (yr-yl)

    do i = 1, num_part
      idx   = int(part_new(i)%Pos(3) / dz)   + 1
      idx_x = int(part_new(i)%Pos(1) / dx)   + 1
      idx_y = int(part_new(i)%Pos(2) / dy)   + 1

      if(idx < Nz+1.and. (part_new(i)%Pos(1).gt.xl).and.(part_new(i)%Pos(1).le.xr) .and. &
            (part_new(i)%Pos(2).gt.yl).and.(part_new(i)%Pos(2).le.yr)) then
        N_j(idx) = N_j(idx) + 1
        local_mass(idx) = local_mass(idx) + part_new(i)%scalar_var(1)*part_new(i)%scalar_var(2)
      endif
    enddo

    call MPI_Allreduce(N_j       ,local_mpi_int_buffer ,Nz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    call MPI_Allreduce(local_mass,local_mpi_real_buffer,Nz,MPI_REAL8  ,MPI_SUM,MPI_COMM_WORLD,mpi_err)

    N_j        = local_mpi_int_buffer
    local_mass = local_mpi_real_buffer

    !mass per layer
    local_mass = dz*dA*local_mass/(dfloat(N_j)+1.d-13)

    deallocate( local_mpi_int_buffer, local_mpi_real_buffer)

  end subroutine compute_particle_mass_subdomain


  function get_activity(w_in,q_in)
    implicit none

    real(8), intent(in) :: w_in, q_in
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
