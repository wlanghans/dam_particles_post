Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: N_avg = 5

  real(8), parameter :: dt = 60.d0, dz=50.
  real(8), parameter :: dA = 4.d4**2.d0

  real(8), parameter :: Det_timer = 10.d0
  real(8), parameter :: trigger_start=900.0, trigger_end=1800.,max_height=4000., min_height=150.

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in, int_tmp
  integer :: sum_Np, Nt_tot, N_active, Nz
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

  real(8) :: pi, dist, r_tmp(10), cell_volume, time

  call initialize_mpi
  call initialize_core_files
  call construct_particle


  call Read_Data(1, part_new)

  do i = 1, num_part
    part_new(i)%inactive_time = 0.
    part_new(i)%activity = .false.
  enddo

  N_active=0
  Nz = int(max_height/dz) 

  do step_out = 2, N_step

    time=step_out*dt
    part_old = part_new
    call Read_Data(step_out, part_new)
 
    
    do i = 1, num_part
      ! new initiation near sfc
      if ((.not.part_old(i)%activity)&
    .and.(part_old(i)%Vel(3).le.0.0.and.part_new(i)%Vel(3).gt.0.0).and. & ! w becomes positive
         (time.ge.trigger_start.and.time.le.trigger_end).and.           & ! within triggering time window   
         part_new(i)%Pos(3).le.min_height) then                           ! below 150 m
         part_new(i)%activity=.true. 
         part_new(i)%inactive_time=time
      end if
     
      ! initiated particle keeps rising actively, gets deactivated otherwise
      if (part_old(i)%activity.and.part_new(i)%Vel(3).ge.0.0) then 
          part_new(i)%activity=.true.
      elseif (part_old(i)%activity.and.part_new(i)%Vel(3).lt.0.0) then
          part_new(i)%activity=.false.
      end if
   
      !check if active particle is considered as triggered
      if (part_new(i)%activity.and.part_new(i)%Pos(3).ge.max_height) then
        ! swap sign to negative to indicate triggered particle
        if (part_new(i)%inactive_time.ge.0.0) part_new(i)%inactive_time=-part_new(i)%inactive_time
      end if
      
      ! count active but untriggered particles at end of simulation
      if (step_out.eq.N_step.and.part_new(i)%activity.and.part_new(i)%inactive_time.gt.0.0) N_active=N_active+1
    end do

  end do

  call MPI_Allreduce(N_active,int_tmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      N_active = int_tmp
  if (masterproc) write(*,*) int_tmp ' particles are active and untriggered at end of simulation'

  ! now get sum up all vertical profiles for all triggered particles
  do step_out=2,N_step
      do i = 1, num_part
        if(part_new(i)%inactive_time.lt.0.0) then
           
        end if
     
      end do
  end do


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

