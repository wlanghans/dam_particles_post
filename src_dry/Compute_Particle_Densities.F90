Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  real(8), parameter :: dt = 60.d0, dz=50.
  real(8), parameter :: dA = 4.d4**2.d0
  real(8), parameter :: max_height=5000.


  real(8) :: xl, xr, yl, yr


  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in
  integer :: sum_Np, Nt_tot
  integer :: N_Count, idx, idx_x, idx_y,ierr

  character(140):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_step,char_dA, param_in
  integer :: output_ncid
  integer :: dimids(2)

  integer :: t_varid, z_varid, zi_varid
  integer, dimension(3) :: var_out_id

  real(8), allocatable, dimension(:) :: sum_vertical_velocities,sum_vertical_velocities_cloud
  integer, allocatable, dimension(:) :: Particle_Number

  real :: time

  call initialize_mpi
  call initialize_core_files
  call construct_particle


  call getarg(3,param_in)
  step_int = char2int( trim(adjustl(param_in)) )

  Nz = int(max_height/dz)
  N_Count=0

  allocate( Particle_Number(Nz) )
  allocate( sum_vertical_velocities(Nz) )
  allocate( sum_vertical_velocities_cloud(Nz) )

  nstep=dt/dt_euler
  do step_out=1,N_Step
    nstep
    if (mod(nstep,step_int).eq.0) then


       N_Count=N_Count+1
    end if

  end do



  if(root) then
    write(*,*) "Writing Data"

    write(char_dz,'(i16)') int(dz)
    fname_out = trim(adjustl(output_dir))//'Particle_Number_Velocity_dz_'//trim(adjustl(char_dz))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"time",N_Count,dimids(2)) )
    call check_nc( nf90_def_dim(output_ncid,"zi",Nz+1,dimids(3)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(2),t_varid) )
    call check_nc( nf90_def_var(output_ncid,"zi",nf90_double,dimids(3),zi_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"             ,nf90_double,dimids(1:2),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"SUM_W"                  ,nf90_double,dimids(1:2),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"SUM_W_CLOUD"            ,nf90_double,dimids(1:2),var_out_id(3) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, zi_varid,(/(dfloat(i)*dz ,i=0,Nz   )/) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt ,i=1,N_Count )/) ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(Particle_Number)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Number_Density) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Specific_Number_Density) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

