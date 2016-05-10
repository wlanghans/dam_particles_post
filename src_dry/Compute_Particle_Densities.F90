Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  real(8), parameter :: dt = 60.d0, dz=500.
  real(8), parameter :: dA = 4.d4**2.d0
  real(8), parameter :: max_height=5000.


  real(8) :: xl, xr, yl, yr


  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in, Nz, step_start
  integer :: sum_Np, Nt_tot, Nz_ind
  integer :: N_Count, idx, idx_x, idx_y,ierr

  character(140):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_step,char_dA, param_in
  integer :: output_ncid
  integer :: dimids(2)

  integer :: t_varid, z_varid, zi_varid
  integer, dimension(3) :: var_out_id

  real(8), allocatable, dimension(:,:) :: sum_vertical_velocities,sum_vertical_velocities_cloud
  real(8), allocatable, dimension(:) :: local_mpi_real_buffer
  integer, allocatable, dimension(:,:) :: Particle_Number
  integer, allocatable, dimension(:) :: local_mpi_int_buffer

  integer :: lag_per_snap

  real :: time, time_firstsnap

  call initialize_mpi
  call initialize_core_files
  call construct_particle


  Nz = int(max_height/dz)
  if (root) write(*,*) 'Using nz = ',Nz,' vertical grid layers'



 
  time_firstsnap=150.*4.
  step_start = INT(time_firstsnap/dt)
  N_Count = 0 
  lag_per_snap = 10
  Nt_tot = (N_Step - step_start) / lag_per_snap + 1

  allocate( Particle_Number(Nz,Nt_tot) )
  allocate( sum_vertical_velocities(Nz,Nt_tot) )
  allocate( sum_vertical_velocities_cloud(Nz,Nt_tot) )
  allocate( local_mpi_int_buffer(Nz) )
  allocate( local_mpi_real_buffer(Nz) )

  Particle_Number = 0
  sum_vertical_velocities = 0.
  sum_vertical_velocities_cloud = 0.

  do step_out=step_start,N_Step,lag_per_snap
       N_Count = N_Count + 1

       call Read_Data(step_out, part_new)

       do i = 1, num_part
         
         Nz_ind = INT(part_new(i)%Pos(3)/dz) + 1
    
         if (Nz_ind.le.Nz) then
           Particle_Number (Nz_ind,N_Count) = Particle_Number (Nz_ind,N_Count) + 1
           sum_vertical_velocities (Nz_ind,N_Count) = sum_vertical_velocities (Nz_ind,N_Count) + part_new(i)%Vel(3)
           if ((part_new(i)%scalar_var(4)+part_new(i)%scalar_var(6)).gt.1.e-5) &
           sum_vertical_velocities_cloud (Nz_ind,N_Count) = sum_vertical_velocities_cloud (Nz_ind,N_Count) + part_new(i)%Vel(3)
         end if

       end do

        ! get sum
       call MPI_Allreduce(Particle_Number(:,N_Count),local_mpi_int_buffer,Nz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
       Particle_Number(:,N_Count)  = local_mpi_int_buffer
       call MPI_Allreduce(sum_vertical_velocities(:,N_Count),local_mpi_real_buffer,Nz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
       sum_vertical_velocities(:,N_Count) = local_mpi_real_buffer
       call MPI_Allreduce(sum_vertical_velocities_cloud(:,N_Count),local_mpi_real_buffer,Nz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
       sum_vertical_velocities_cloud(:,N_Count) = local_mpi_real_buffer


  end do





  if(root) then
    write(*,*) "Writing Data"

    write(char_dz,'(i16)') int(dz)
    fname_out = trim(adjustl(output_dir))//'Particle_Number_Velocity_dz_'//trim(adjustl(char_dz))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"time",Nt_tot,dimids(2)) )
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
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt,i=step_start,N_Step,lag_per_snap)/) ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(Particle_Number)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2),sum_vertical_velocities) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3),sum_vertical_velocities_cloud) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

