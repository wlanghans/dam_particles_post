program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  real(8), parameter :: dt = 60.d0
  real(8), parameter :: dA = 4.d4**2.d0, dz=0., dx=100., dy=100.
  real(8), parameter :: max_height=5000.
  integer, parameter :: nzm=115, Nx=400, Ny=400


  real(8) :: xl, xr, yl, yr, max_height_new


  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in, Nz, step_start
  integer :: sum_Np, Nt_tot, Nz_ind, Nx_ind, Ny_ind
  integer :: N_Count, idx, idx_x, idx_y,ierr

  character(140):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_step,char_dA, param_in
  integer :: output_ncid
  integer :: dimids(6)

  integer :: t_varid, z_varid, zi_varid, x_varid, y_varid, zs_varid
  integer, dimension(5) :: var_out_id

  real(8), allocatable, dimension(:,:) :: sum_vertical_velocities,sum_vertical_velocities_cloud
  real(8), allocatable, dimension(:) :: local_mpi_real_buffer, z, zi, dz_vector
  integer, allocatable, dimension(:,:) :: Particle_Number
  integer, allocatable, dimension(:,:,:,:) :: Nxy
  integer, allocatable, dimension(:) :: local_mpi_int_buffer

  integer :: lag_per_snap

  real :: time, time_firstsnap

  call initialize_mpi
  call initialize_core_files
  call construct_particle


  allocate(z(nzm))
  allocate(dz_vector(nzm))
  allocate(zi(nzm+1))
  
  call set_grid(z,zi,dz,dz_vector,nzm)

  Nz = minloc(abs(max_height-zi),1) - 1
  max_height_new = zi(Nz+1)
  if (root) write(*,*) 'Using Nz = ',Nz,' vertical grid layers with ztop= ',max_height_new

 
  time_firstsnap=0.*4.
  step_start = INT(time_firstsnap/dt) + 1
  N_Count = 0 
  lag_per_snap = 10
  Nt_tot = (N_Step - step_start) / lag_per_snap + 1

  allocate( Particle_Number(Nz,Nt_tot) )
  allocate( Nxy(Nx,Ny,10,Nt_tot) )
  allocate( sum_vertical_velocities(Nz,Nt_tot) )
  allocate( sum_vertical_velocities_cloud(Nz,Nt_tot) )
  allocate( local_mpi_int_buffer(Nz) )
  allocate( local_mpi_real_buffer(Nz) )

  Particle_Number = 0
  Nxy = 0 
  sum_vertical_velocities = 0.
  sum_vertical_velocities_cloud = 0.

  do step_out=step_start,N_Step,lag_per_snap
       N_Count = N_Count + 1

       call Read_Data(step_out, part_new)

       do i = 1, num_part
         
         Nz_ind = minloc(abs(part_new(i)%Pos(3) - zi),1 )
         Nx_ind = INT(part_new(i)%Pos(1)/dx) + 1
         Ny_ind = INT(part_new(i)%Pos(2)/dy) + 1
         if (part_new(i)%Pos(3).ge.zi(Nz_ind)) then
            Nz_ind = Nz_ind
         else
            Nz_ind = Nz_ind -1
         end if 
    
         if (Nz_ind.le.10) then
           Nxy(Nx_ind,Ny_ind,Nz_ind,N_Count) = Nxy(Nx_ind,Ny_ind,Nz_ind,N_Count) + 1
         end if    
 
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
    call check_nc( nf90_def_dim(output_ncid,"x",Nx,dimids(4)) )
    call check_nc( nf90_def_dim(output_ncid,"y",Ny,dimids(5)) )
    call check_nc( nf90_def_dim(output_ncid,"z_short",10,dimids(6)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(2),t_varid) )
    call check_nc( nf90_def_var(output_ncid,"zi",nf90_double,dimids(3),zi_varid) )
    call check_nc( nf90_def_var(output_ncid,"x",nf90_double,dimids(4),x_varid) )
    call check_nc( nf90_def_var(output_ncid,"y",nf90_double,dimids(5),y_varid) )
    call check_nc( nf90_def_var(output_ncid,"z_short",nf90_double,dimids(6),zs_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"             ,nf90_double,dimids(1:2),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"SUM_W"                  ,nf90_double,dimids(1:2),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"SUM_W_CLOUD"            ,nf90_double,dimids(1:2),var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"dz_vector"              ,nf90_double,dimids(1),var_out_id(4) ) )
    call check_nc( nf90_def_var( output_ncid,"N_Particle_3D"          ,nf90_double,dimids((/4,5,6,2/)),var_out_id(5) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, z(1:Nz) ) )
    call check_nc( nf90_put_var(output_ncid, zi_varid, zi(1:Nz+1) ) )
    call check_nc( nf90_put_var(output_ncid, zs_varid, z(1:10) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i-1)*dt,i=step_start,N_Step,lag_per_snap)/) ) )
    call check_nc( nf90_put_var(output_ncid, x_varid, (/(dfloat(i-1)*dx+0.5*dx,i=1,Nx)/) ) )
    call check_nc( nf90_put_var(output_ncid, y_varid, (/(dfloat(i-1)*dy+0.5*dy,i=1,Ny)/) ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(Particle_Number)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2),sum_vertical_velocities) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3),sum_vertical_velocities_cloud) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4),dz_vector(1:Nz)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(5),dfloat(Nxy)) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi


end Program 

  subroutine set_grid(z,zi,dz,dz_vector,nzm)

  use particle_data
  use mpi_info
  use core_info

  implicit none
  integer, intent(in) :: nzm
  real(8), dimension(nzm), intent(out)   :: z, dz_vector
  real(8), dimension(nzm+1), intent(out) :: zi
  real(8), intent(in) :: dz

  integer :: ios, k
 
 
      if (dz.le.0.) then

      open(8,file=trim(path)//'/grd',status='old',form='formatted',iostat=ios)
            if (ios .gt. 0) then
               if(root) write(*,*) 'Error in set_grid: Unable to open grd file'
               call finalize_mpi
            end if

            ios = 0
            k = 0
            do while (ios .le. 0 .and. k .lt. nzm)
               k = k + 1
               read(8,fmt=*,iostat=ios) z(k)
               if (ios .gt. 0 .and. k .le. 2) then
                  if(root) write(*,*) 'Error in set_grid: The grd file must contain at least two heights'
                  call finalize_mpi
               end if
            end do
            do while (k .le. nzm)
               z(k) = 2.*z(k-1) - z(k-2)
               k = k + 1
            end do
            close(8)

            ! Check that z is monotonically increasing
            if (any(z(2:nzm) .le. z(1:nzm-1))) then
               if(root) write(*,*) 'Error in set_grid: The heights in the grd file must be monotonically increasing'
               call finalize_mpi
            end if

      else
          z(1) = 0.5*dz
            do k = 2,nzm
               z(k) = z(k-1) + dz
            end do
      end if

         dz_vector(1) = 0.5*(z(1)+z(2))
         do k = 2,nzm-1
            dz_vector(k) = 0.5*( z(k+1) - z(k-1) )
         end do
         dz_vector(nzm) = z(nzm) - z(nzm-1)

         zi(1) = 0.
         do k = 2,nzm
            zi(k) = 0.5*( z(k-1) + z(k) )
         end do
         zi(nzm+1) = zi(nzm) + dz_vector(nzm)



  end subroutine set_grid

