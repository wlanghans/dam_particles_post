Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: Nz = 300
  integer, parameter :: N_avg = 15

  real(8), parameter :: dz = 50.d0, dt = 4.d0
  real(8), parameter :: dA = 2.d4**2.d0

  real(8), parameter :: Det_timer = 10.d0

  integer, parameter :: Nx = 200, Ny = 200
  real(8), parameter :: dx = 100.d0, dy = 100.d0
  real(8), parameter :: A_threshold = 0.5d0

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in
  integer :: sum_Np, Nt_tot
  integer :: N_Count, idx, idx_x, idx_y,ierr

  character(128):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_w,char_tmp
  integer :: output_ncid
  integer :: dimids(2)

  integer :: t_varid, z_varid
  integer, dimension(10) :: var_out_id

  real(8), allocatable, dimension(:,:) :: Cell_Particle, Entrain, Entrain_upp, Entrain_low, Lv_flux, inv_den, Perimeter

  integer, dimension(Nx,Ny,Nz) :: Act_3D, Cell_P_3D, Buffer_3D
  real(8), dimension(Nx,Ny,Nz) :: R_Act_3D

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), cell_volume

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  Nt_tot = N_step/N_avg - 1

  allocate( Cell_Particle(Nz,Nt_tot) )
  allocate( Entrain(Nz,Nt_tot) )
  allocate( Entrain_upp(Nz,Nt_tot) )
  allocate( Entrain_low(Nz,Nt_tot) )
  allocate( Lv_flux(Nz,Nt_tot) )
  allocate( inv_den(Nz,Nt_tot) )
  allocate( Perimeter(Nz,Nt_tot) )

  call Read_Data(1, part_new)
  call compute_particle_mass(Nz,dz,dA)

  if(root) write(*,*) 'w_cutoff',w_cutoff
  if(root) write(*,*) dmass

  Cell_Particle = 0.d0
  Entrain = 0.d0
  Entrain_low= 0.d0
  Entrain_upp= 0.d0
  Lv_flux = 0.d0
  inv_den = 0.d0
  Perimeter = 0.d0

  N_Count = 0
  do i = 1, num_part
    part_new(i)%inactive_time = 1.d6
  enddo

  step_cnt = 1
  do step_out = 1, Nt_tot
    Act_3D = 0
    Cell_P_3D = 0
    do step_in = 1, N_avg
      step_cnt = step_cnt + 1
      part_old = part_new

      call Read_Data(step_cnt, part_new)
      do i = 1, num_part
        part_new(i)%activity = get_activity( part_new(i)%Vel(3), part_new(i)%scalar_var(4) + part_new(i)%scalar_var(6) )
 
        part_new(i)%inactive_time = part_old(i)%inactive_time + dt
        r_tmp(1) = .5d0*( part_old(i)%Pos(3)+part_new(i)%Pos(3) )
        idx = int( r_tmp(1)/dz ) + 1

        if(idx < Nz+1) then

          idx_x = int(.5d0*(part_new(i)%Pos(1) + part_old(i)%Pos(1)) / dx) + 1
          idx_y = int(.5d0*(part_new(i)%Pos(2) + part_old(i)%Pos(2)) / dy) + 1
          Cell_P_3D(idx_x,idx_y,idx) = Cell_P_3D(idx_x,idx_y,idx) + 1
 
          if(part_new(i)%activity) then

            Act_3D(idx_x,idx_y,idx) = Act_3D(idx_x,idx_y,idx) + 1

            if(.not.part_old(i)%activity) &
              Entrain_upp(idx,step_out) = Entrain_upp(idx,step_out) + 1.d0
 
            if(.not.part_old(i)%activity .and. (part_new(i)%inactive_time > Det_timer) ) &
              Entrain(idx,step_out) = Entrain(idx,step_out) + 1.d0
 
            if(.not.part_old(i)%activity .and. (part_new(i)%inactive_time > 1.d4) ) &
              Entrain_low(idx,step_out) = Entrain_low(idx,step_out) + 1.d0
 
            part_new(i)%inactive_time = 0.d0
            Cell_Particle(idx,step_out) = Cell_Particle(idx,step_out) + 1.d0

            Lv_flux(idx,step_out) = Lv_flux(idx,step_out) + max(w_cutoff, .5d0*( part_new(i)%Vel(3) + part_old(i)%Vel(3) ) )
            inv_den(idx,step_out) = inv_den(idx,step_out) + 2.d0/( part_new(i)%Scalar_Var(1)*part_new(i)%Scalar_var(2) &
                                                                  +part_old(i)%Scalar_Var(1)*part_old(i)%Scalar_var(2) )

          endif
        endif
      enddo
    enddo

    call MPI_Allreduce(Act_3D,Buffer_3D,Nx*Ny*Nz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    Act_3D = Buffer_3D

    call MPI_Allreduce(Cell_P_3D,Buffer_3D,Nx*Ny*Nz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    Cell_P_3D = Buffer_3D

    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          R_Act_3D(i,j,k) = dfloat(Act_3D(i,j,k))/(dfloat(Cell_P_3D(i,j,k)) + 1.d-15)
        enddo
      enddo
    enddo

    do k = 1, Nz
      do j = my_id+2, Ny-1, Nproc
        do i = 2, Nx-1
          if(R_Act_3D(i,j,k) > A_threshold) then
            if(R_Act_3D(i-1,j  ,k) < A_threshold) Perimeter(k,step_out) = Perimeter(k,step_out) + dy
            if(R_Act_3D(i+1,j  ,k) < A_threshold) Perimeter(k,step_out) = Perimeter(k,step_out) + dy
            if(R_Act_3D(i  ,j-1,k) < A_threshold) Perimeter(k,step_out) = Perimeter(k,step_out) + dx
            if(R_Act_3D(i  ,j+1,k) < A_threshold) Perimeter(k,step_out) = Perimeter(k,step_out) + dx
          endif
        enddo
      enddo
    enddo

    if(root) write(*,*) step_out, step_cnt
  enddo  

  call Reduce_Mtx(Cell_Particle)
  call Reduce_Mtx(Lv_flux)
  call Reduce_Mtx(inv_den)
  call Reduce_Mtx(Entrain)
  call Reduce_Mtx(Entrain_low)
  call Reduce_Mtx(Entrain_upp)
  call Reduce_Mtx(Perimeter)

  if(root) then
    write(*,*) "Writing Data"

    write(char_tmp,'(i16)') int(N_Avg* (dt + 1.d-10) )
    write(char_dz,'(i16)') int(dz)
    if(w_cutoff >= 0.d0) then
      write(char_w,'(i16)') int(w_cutoff)
    else
      char_w = 'None'
    endif
    fname_out = trim(adjustl(output_dir))//'Lagrangian_Entrainment_dz_'//trim(adjustl(char_dz))//'_Tavg_'//trim(adjustl(char_tmp))//'_W_'//trim(adjustl(char_w))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"time",Nt_tot,dimids(2)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(2),t_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"        ,nf90_double,dimids,var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment_10_sec",nf90_double,dimids,var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment_upp_bd",nf90_double,dimids,var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment_low_bd",nf90_double,dimids,var_out_id(4) ) )
    call check_nc( nf90_def_var( output_ncid,"Vertical_flux"     ,nf90_double,dimids,var_out_id(5) ) )
    call check_nc( nf90_def_var( output_ncid,"inv_density"       ,nf90_double,dimids,var_out_id(6) ) )
    call check_nc( nf90_def_var( output_ncid,"Perimeter"         ,nf90_double,dimids,var_out_id(7) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt*N_avg +.5d0*dt*N_avg,i=0,Nt_tot-1)/) ) )

    cell_volume = dz*dA

    Cell_Particle = Cell_Particle/dfloat(N_avg)
    Lv_flux = Lv_flux/dfloat(N_avg)
    inv_den = inv_den/dfloat(N_avg)
    Entrain = Entrain/dfloat(N_avg)
    Entrain_low = Entrain_low/dfloat(N_avg)
    Entrain_upp = Entrain_upp/dfloat(N_avg)

    Lv_flux = dmass*Lv_flux/ cell_volume
    inv_den = dmass*inv_den/ cell_volume

    Entrain     = dmass*Entrain    /(cell_volume*dt)
    Entrain_low = dmass*Entrain_low/(cell_volume*dt)
    Entrain_upp = dmass*Entrain_upp/(cell_volume*dt)

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), Cell_Particle) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Entrain) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Entrain_upp) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4), Entrain_low) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(5), Lv_flux) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(6), inv_den) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(7), Perimeter) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

