Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: Nz = 164
  integer, parameter :: N_avg = 10

  real(8), parameter :: dz = 50.d0, dt = 120.d0
  real(8), parameter :: dA = 51.2d3**2.d0

  real(8), parameter :: Det_timer = 10.d0

  integer, parameter :: Nx = 256, Ny = 256
  real(8), parameter :: dx = 200.d0, dy = 200.d0
  real(8), parameter :: A_threshold = 0.5d0

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in, natm_jump
  integer :: sum_Np, Nt_tot
  integer :: idx, idx_x, idx_y,ierr

  character(128):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_tmp
  integer :: output_ncid
  integer :: dimids(4)

  integer :: t_varid, t2_varid, x_varid, y_varid
  integer, dimension(2) :: var_out_id

  real(8), allocatable, dimension(:,:,:) :: qvflux, aqvflux
  real(8), allocatable, dimension(:) :: gqvflux,local_mpi_real_buffer

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10)

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  Nt_tot = N_step/N_avg - 1

  allocate( qvflux(Nx,Ny,Nt_tot+1) )
  allocate( aqvflux(Nx,Ny,Nt_tot+1) )
  allocate( gqvflux(Nt_tot), local_mpi_real_buffer(Nt_tot))

  call Read_Data(1, part_new)
  !call compute_particle_mass(Nz,dz,dA)
  call getarg(3,param)
  read(param,*)dmass

  if(root) write(*,*) dmass

  qvflux = 0.0d0
  aqvflux = 0.0d0
  gqvflux = 0.0d0


  step_cnt = 1
  do step_out = 1, Nt_tot
    do step_in = 1, N_avg
      step_cnt = step_cnt + 1
      part_old = part_new

      call Read_Data(step_cnt, part_new)
      do i = 1, num_part

        if (part_new(i)%natm.gt.part_old(i)%natm) then

          natm_jump = part_new(i)%natm - part_old(i)%natm
          idx_x = int(part_new(i)%Pos(1) / dx) + 1
          idx_y = int(part_new(i)%Pos(2) / dy) + 1
          gqvflux(step_out) = gqvflux(step_out)  + dfloat(natm_jump)
          aqvflux(idx_x,idx_y,step_out+1) = aqvflux(idx_x,idx_y,step_out+1) + dfloat(natm_jump)
          if (step_in.eq.1 .and. step_out.eq.1) qvflux(idx_x,idx_y,step_out) = qvflux(idx_x,idx_y,step_out) + dfloat(natm_jump)
          if (step_in.eq.N_avg) qvflux(idx_x,idx_y,step_out+1) = qvflux(idx_x,idx_y,step_out+1) + dfloat(natm_jump)
 
        endif !natm increased
      enddo
    enddo

    if (step_out.eq.1) then
    call Reduce_Mtx(qvflux(1:Nx,1:Ny,step_out))
    call Reduce_Mtx(qvflux(1:Nx,1:Ny,step_out+1))
    else
    call Reduce_Mtx(qvflux(1:Nx,1:Ny,step_out+1))
    end if
    call Reduce_Mtx(aqvflux(1:Nx,1:Ny,step_out+1))

    if(root) write(*,*) step_out, step_cnt
  enddo  

  call MPI_Reduce(gqvflux,local_mpi_real_buffer,Nt_tot,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
  gqvflux= local_mpi_real_buffer
  gqvflux= dmass*gqvflux/dA/(dfloat(N_avg)*dt)
  qvflux = dmass*qvflux / (dx*dy*dt) 
  aqvflux = dmass*aqvflux / (dx*dy) / (dfloat(N_avg)*dt) 


  if(root) then
    write(*,*) "Writing Data"

    write(char_tmp,'(i16)') int(N_Avg* (dt + 1.d-10) )
    fname_out = trim(adjustl(output_dir))//'Surface_vaporflux_Tavg_'//trim(adjustl(char_tmp))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"x",Nx,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"y",Ny,dimids(2)) )
    call check_nc( nf90_def_dim(output_ncid,"time",Nt_tot+1,dimids(3)) )
    call check_nc( nf90_def_dim(output_ncid,"time_2",Nt_tot,dimids(4)) )

    call check_nc( nf90_def_var(output_ncid,"x",nf90_double,dimids(1),x_varid) )
    call check_nc( nf90_def_var(output_ncid,"y",nf90_double,dimids(2),y_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(3),t_varid) )
    call check_nc( nf90_def_var(output_ncid,"time_2",nf90_double,dimids(4),t2_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"dvflux"        ,nf90_double,dimids(1:3),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"global_dvflux",nf90_double,dimids(4),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"advflux",nf90_double,dimids(1:3),var_out_id(3) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, x_varid, (/(dfloat(i)*dx ,i=0,Nx    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, y_varid, (/(dfloat(i)*dy ,i=0,Ny    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt*N_avg/3600.d0/24.d0 ,i=0,Nt_tot)/) ) )
    call check_nc( nf90_put_var(output_ncid, t2_varid, (/((dfloat(i)*dt*N_avg +.5d0*dt*N_avg)/3600.d0/24.d0,i=0,Nt_tot-1)/) ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), qvflux) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), gqvflux) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), aqvflux) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

