Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: N_avg = 5

  real(8), parameter :: dz = 500.d0, dt = 4.d0
  real(8), parameter :: dA = 2.d4**2.d0

  real(8), parameter :: Det_timer = 10.d0

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

  real(8), allocatable, dimension(:) :: Cell_Particle, Entrain, Lv_flux, inv_den, Center_Mass

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), cell_volume, max_height

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  Nt_tot = N_step/N_avg - 1

  allocate( Cell_Particle(Nt_tot) )
  allocate( Entrain(Nt_tot) )
  allocate( Lv_flux(Nt_tot) )
  allocate( inv_den(Nt_tot) )
  allocate( Center_Mass(Nt_tot) )

  call Read_Data(1, part_new)
  call compute_particle_mass(300,50.d0,dA)

  if(root) write(*,*) 'w_cutoff',w_cutoff
  if(root) write(*,*) dmass

  Cell_Particle = 0.d0
  Entrain = 0.d0
  Lv_flux = 0.d0
  inv_den = 0.d0
  Center_Mass = 0.d0

  N_Count = 0
  do i = 1, num_part
    part_new(i)%inactive_time = 1.d6
  enddo

  step_cnt = 1
  do step_out = 1, Nt_tot

    max_height = 0.d0
    do step_in = 1, N_avg
      step_cnt = step_cnt + 1
      part_old = part_new

      call Read_Data(step_cnt, part_new)
      do i = 1, num_part
        part_new(i)%activity = get_activity( part_new(i)%Vel(3), part_new(i)%scalar_var(4) + part_new(i)%scalar_var(6) )
 
        part_new(i)%inactive_time = part_old(i)%inactive_time + dt

        if(part_new(i)%activity) &
          max_height = max(max_height, part_new(i)%Pos(3) )
      enddo

      call MPI_Allreduce(max_height,r_tmp(1),1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_height = r_tmp(1)
 
      do i = 1, num_part
        if(part_new(i)%activity) then

          if(part_new(i)%Pos(3) > max_height - dz) then

            Cell_Particle(step_out) = Cell_Particle(step_out) + 1.d0
            Center_Mass(step_out) = Center_Mass(step_out) + part_new(i)%Pos(3)
            if(.not.part_old(i)%activity .and. (part_new(i)%inactive_time > Det_timer) ) &
              Entrain(step_out) = Entrain(step_out) + 1.d0

            Lv_flux(step_out) = Lv_flux(step_out) + max(w_cutoff, .5d0*( part_new(i)%Vel(3) + part_old(i)%Vel(3) ) )
            inv_den(step_out) = inv_den(step_out) + 2.d0/( part_new(i)%Scalar_Var(1)*part_new(i)%Scalar_var(2) &
                                                          +part_old(i)%Scalar_Var(1)*part_old(i)%Scalar_var(2) )

          endif
          part_new(i)%inactive_time = 0.d0

        endif
      enddo

    enddo

    if(root) write(*,*) step_out, step_cnt
  enddo

  call Reduce_Mtx(Cell_Particle)
  call Reduce_Mtx(Center_Mass)
  call Reduce_Mtx(Lv_flux)
  call Reduce_Mtx(inv_den)
  call Reduce_Mtx(Entrain)

  Center_Mass = Center_Mass/(Cell_Particle + 1.d-15)

  if(root) then
    write(*,*) "Writing Data"

    write(char_tmp,'(i16)') int(N_Avg* (dt + 1.d-10) )
    write(char_dz,'(i16)') int(dz)
    write(char_w,'(i16)') int(w_cutoff)
    fname_out = trim(adjustl(output_dir))//'CT_Entrainment_dz_'//trim(adjustl(char_dz))//'_Tavg_'//trim(adjustl(char_tmp))//'_W_'//trim(adjustl(char_w))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"time",Nt_tot,dimids(1)) )

    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(1),t_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"        ,nf90_double,dimids(1),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment_10_sec",nf90_double,dimids(1),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Vertical_flux"     ,nf90_double,dimids(1),var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"inv_density"       ,nf90_double,dimids(1),var_out_id(4) ) )
    call check_nc( nf90_def_var( output_ncid,"Center_Mass"       ,nf90_double,dimids(1),var_out_id(5) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt*N_avg +.5d0*dt*N_avg,i=0,Nt_tot-1)/) ) )

    cell_volume = dz*dA

    Cell_Particle = Cell_Particle/dfloat(N_avg)
    Lv_flux = Lv_flux/dfloat(N_avg)
    inv_den = inv_den/dfloat(N_avg)
    Entrain = Entrain/dfloat(N_avg)

    Lv_flux = dmass*Lv_flux/ cell_volume
    inv_den = dmass*inv_den/ cell_volume

    Entrain     = dmass*Entrain    /(cell_volume*dt)

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), Cell_Particle) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Entrain) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Lv_flux) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4), inv_den) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(5), Center_Mass) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

