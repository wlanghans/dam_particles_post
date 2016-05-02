Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: Nz = 300
  integer, parameter :: N_avg = 5

  real(8), parameter :: dz = 50.d0, dt = 4.d0
  real(8), parameter :: dA = 2.d4**2.d0

  real(8), parameter :: Det_timer = dfloat(N_avg)*dt

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in
  integer :: sum_Np, Nt_tot
  integer :: N_Count, idx

  integer :: output_ncid
  integer :: dimids(2)

  integer :: t_varid, z_varid
  integer, dimension(10) :: var_out_id

  real(8), allocatable, dimension(:,:) :: Cell_Particle, Entrain, Lv_flux, inv_den

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), cell_volume

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  Nt_tot = N_step/N_avg - 1

  allocate( Cell_Particle(Nz,Nt_tot) )
  allocate( Entrain(Nz,Nt_tot) )
  allocate( Lv_flux(Nz,Nt_tot) )
  allocate( inv_den(Nz,Nt_tot) )

  call Read_Data(1, part_new)
  call compute_particle_mass(Nz,dz,dA)

  Cell_Particle = 0.d0
  Entrain = 0.d0
  Lv_flux = 0.d0
  inv_den = 0.d0

  N_Count = 0
  do i = 1, num_part
    part_new(i)%inactive_time = 10000.d0
  enddo

  step_cnt = 1
  do step_out = 1, Nt_tot
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
          if(part_new(i)%activity) then

            if(.not.part_old(i)%activity .and. (part_new(i)%inactive_time > Det_timer) ) &
              Entrain(idx,step_out) = Entrain(idx,step_out) + 1.d0
 
            part_new(i)%inactive_time = 0.d0
            Cell_Particle(idx,step_out) = Cell_Particle(idx,step_out) + 1.d0

            Lv_flux(idx,step_out) = Lv_flux(idx,step_out) + max(w_cutoff, .5d0*( part_new(i)%Vel(3) + part_old(i)%Vel(3) ) )
            inv_den(idx,step_out) = inv_den(idx,step_out) + 2.d0/( part_new(i)%Scalar_Var(1)*part_new(i)%Scalar_var(2) &
                                                                  +part_old(i)%Scalar_Var(1)*part_old(i)%Scalar_var(2) )

          endif
        endif
      enddo
    enddo
    if(my_id ==0) write(*,*) step_out, step_cnt
  enddo  

  call Reduce_Mtx(Cell_Particle)
  call Reduce_Mtx(Lv_flux)
  call Reduce_Mtx(inv_den)
  call Reduce_Mtx(Entrain)

  if(root) then
    write(*,*) "Writing Data"

    write(char_tmp,'(i5)') int(Det_timer)
    char_tmp = adjustl(char_tmp)
    fname_out = 'Entrainment_'//trim(char_tmp)//'_sec.nc'
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"time",Nt_tot,dimids(2)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(2),t_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"   ,nf90_double,dimids,var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment"  ,nf90_double,dimids,var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Vertical_flux",nf90_double,dimids,var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"inv_density"  ,nf90_double,dimids,var_out_id(4) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt*N_avg +.5d0*dt*N_avg,i=0,Nt_tot-1)/) ) )

    cell_volume = dz*dA

    Cell_Particle = Cell_Particle/dfloat(N_avg)
    Lv_flux = Lv_flux/dfloat(N_avg)
    inv_den = inv_den/dfloat(N_avg)
    Entrain = Entrain/dfloat(N_avg)

    Lv_flux = dmass*Lv_flux/ cell_volume
    inv_den = dmass*inv_den/ cell_volume
    Entrain = dmass*Entrain/(cell_volume*dt)

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), Cell_Particle) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Entrain) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Lv_flux) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4), inv_den) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

