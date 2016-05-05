Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  real(8), parameter :: dt = 60.d0, dz=50.
  real(8), parameter :: dA = 4.d4**2.d0

  real(8), parameter :: Det_timer = 10.d0
  real(8), parameter :: trigger_start=900.0, trigger_end=1800.,max_height=5000., min_height=150.

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in, int_tmp
  integer :: sum_Np, Nt_tot, N_active, Nz
  integer :: N_Count, idx, idx_x, idx_y,ierr, Nz_ind, Nz_ind_old
 

  character(128):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_w,char_tmp
  integer :: output_ncid
  integer :: dimids(2)

  integer :: z_varid, zi_varid
  integer, dimension(10) :: var_out_id

  integer, dimension(:), allocatable :: N_triggered, local_mpi_int_buffer
  real(8), dimension(:), allocatable :: Fdyn, Fbuoy, w2, local_mpi_real_buffer, local_mpi_real_buffer2

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), time, dz_part, ztop, fdyn_tmp, fbuoy_tmp, w2_tmp

  call initialize_mpi
  call initialize_core_files
  call construct_particle


  call Read_Data(1, part_new)

  do i = 1, num_part
    part_new(i)%inactive_time = 0.
    part_new(i)%activity = .false.
    part_new(i)%Vel(1) = 0.
  enddo

  N_active=0
  Nz = int(max_height/dz) 
 
  if (root) write(*,*) 'Using nz = ',Nz,' vertical grid layers'

  allocate(N_triggered(Nz+1))
  allocate(w2(Nz+1))
  allocate(Fdyn(Nz))
  allocate(Fbuoy(Nz))
  allocate(local_mpi_int_buffer(Nz+1))
  allocate(local_mpi_real_buffer(Nz))
  allocate(local_mpi_real_buffer2(Nz+1))
  Fdyn  = 0.
  Fbuoy = 0.
  N_triggered=0
  w2=0.
  N_Count=0

  do step_out = 2, N_step

    time=step_out*dt
    part_old = part_new
    call Read_Data(step_out, part_new)
 
    
    do i = 1, num_part

      ! retain triggering time
      part_new(i)%Vel(1) = part_old(i)%Vel(1) 

      ! new initiation below min_height
      if ((.not.part_old(i)%activity)&
         .and.(part_old(i)%Vel(3).le.0.0.and.part_new(i)%Vel(3).gt.0.0).and. & ! w becomes positive
         (time.ge.trigger_start.and.time.le.trigger_end).and.           & ! within triggering time window   
         part_new(i)%Pos(3).lt.min_height) then                           ! below 150 m
         ! make sure it did not get triggered before
         if (part_new(i)%inactive_time.ge.0.) then
           part_new(i)%activity=.true. 
           part_new(i)%inactive_time=time
         end if
      end if
     
      ! initiated particle keeps rising actively, gets deactivated otherwise
      if (part_old(i)%activity.and.part_new(i)%Vel(3).ge.0.0) then 
          part_new(i)%activity=.true.
      elseif (part_old(i)%activity.and.part_new(i)%Vel(3).lt.0.0) then
          part_new(i)%activity=.false.
      end if
   
      !check if active particle is considered as triggered
      if (part_new(i)%activity.and.part_new(i)%Pos(3).gt.max_height) then
        if (part_new(i)%inactive_time.gt.0.0) then 
           ! swap sign to negative to indicate triggered particle
           part_new(i)%inactive_time=-part_new(i)%inactive_time
           ! save triggering time
           part_new(i)%Vel(1) = time
           ! increase total number of triggered particles
           N_Count=N_Count+1
        end if
      end if
      
      ! count active but untriggered particles at end of simulation
      if (step_out.eq.N_step.and.part_new(i)%activity.and.part_new(i)%inactive_time.gt.0.0) N_active=N_active+1
    end do

  end do

  call MPI_Allreduce(N_Count,int_tmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      N_Count = int_tmp
  if (root) write(*,*) N_Count ,' particles are triggered'
  call MPI_Allreduce(N_active,int_tmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      N_active = int_tmp
  if (root) write(*,*) N_active ,' particles are active AND untriggered at end of simulation'


  ! now sum up properties of all trajectories of triggered particles
  do step_out=2,N_step
    time=step_out*dt
    part_old = part_new
    call Read_Data(step_out, part_new)
      do i = 1, num_part

      ! retain triggering time (ie, time first above max_height)
      part_new(i)%Vel(1) = part_old(i)%Vel(1) 

        if(part_new(i)%inactive_time.lt.0.0.and.                                 &    ! got triggered
         time.ge.abs(part_new(i)%inactive_time).and.time.le.part_new(i)%Vel(1)) then  ! data point between
                                                                                      ! initiation and triggering

           Nz_ind     = INT(part_new(i)%Pos(3)/dz) + 1
           Nz_ind_old = INT(part_old(i)%Pos(3)/dz) + 1

           ! particle is initiated in layer
           if (time.eq.abs(part_new(i)%inactive_time)) then
              N_triggered( Nz_ind   ) = N_triggered(Nz_ind) + 1
              dz_part=part_new(i)%pos(3)-(Nz_ind-1)*dz 
              ! lin interpolate forcings to centerpoint (between bottom of dz layer and current Pos(3))
              fbuoy_tmp = part_new(i)%scalar_var(12) - 0.5 * dz_part * & 
              (part_new(i)%scalar_var(12)-part_old(i)%scalar_var(12))/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
              fdyn_tmp = part_new(i)%scalar_var(11) - 0.5 * dz_part * & 
              (part_new(i)%scalar_var(11)-part_old(i)%scalar_var(11))/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
              ! part_old(i)%Vel(3)<=0, thus a plus sign below
              w2_tmp = (part_new(i)%Vel(3))**2 - dz_part * &
              ((part_new(i)%Vel(3))**2+(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
              Fbuoy(Nz_ind)       = Fbuoy(Nz_ind) + dz_part * fbuoy_tmp
              Fdyn(Nz_ind)        = Fdyn(Nz_ind)  + dz_part * fdyn_tmp
              w2(Nz_ind)          = w2(Nz_ind)    + w2_tmp
           end if

           ! particle is still in same layer, then just add work done in layer
           if (Nz_ind_old.eq.Nz_ind) then
             ! again lin interpolate forcings to centerpoint
              dz_part=max(0.,part_new(i)%pos(3)-part_old(i)%pos(3))                
              Fbuoy(Nz_ind)       = Fbuoy(Nz_ind) + &
              0.5*(part_new(i)%scalar_var(12)+part_old(i)%scalar_var(12)) * dz_part
              Fdyn(Nz_ind)        = Fdyn(Nz_ind)  + &
              0.5*(part_new(i)%scalar_var(11)+part_old(i)%scalar_var(11)) * dz_part
           end if

           !particle rises into layer from below, but does not get initiated in
           !this layer
           if (Nz_ind_old.lt.Nz_ind.and.time.ne.abs(part_new(i)%inactive_time)) then
              if (Nz_ind.gt.Nz) then
                 !particle above or at max_height
                 ztop=max_height
                 if (Nz_ind_old.lt.Nz) then
                   N_triggered( Nz_ind_old+1:Nz  ) = N_triggered(Nz_ind_old+1:Nz) + 1
                 end if
                 N_triggered( Nz+1 ) = N_triggered( Nz+1 ) + 1
                 w2_tmp = (part_old(i)%Vel(3))**2 + ( Nz*dz-part_old(i)%Pos(3)) * &
                 ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                 w2(Nz+1) = w2(Nz+1) + w2_tmp
              else
                 !particle below max_height
                 ztop= part_new(i)%Pos(3)
                 N_triggered( Nz_ind_old+1:Nz_ind   ) = N_triggered(Nz_ind_old+1:Nz_ind) + 1
              end if

              !linear interpolate forcings and w2 and incrementally add them to height bins
              do k=Nz_ind_old,min(Nz_ind,Nz)
                 if (k.eq.Nz_ind_old) then
                   dz_part= Nz_ind_old*dz-part_old(i)%Pos(3)
                   fbuoy_tmp = part_old(i)%scalar_var(12) + 0.5 * dz_part *&
                   (part_new(i)%scalar_var(12) - part_old(i)%scalar_var(12))/(part_new(i)%Pos(3) - part_old(i)%Pos(3))
                   fdyn_tmp = part_old(i)%scalar_var(11) + 0.5 * dz_part *&
                   (part_new(i)%scalar_var(11) - part_old(i)%scalar_var(11))/(part_new(i)%Pos(3) - part_old(i)%Pos(3))
                 elseif (k.gt.Nz_ind_old.and.k.lt.Nz_ind) then
                   dz_part=dz
                   fbuoy_tmp = part_old(i)%scalar_var(12) + ( Nz_ind_old*dz-part_old(i)%Pos(3) + (k-Nz_ind_old-1)*dz + 0.5 * dz) *&
                   (part_new(i)%scalar_var(12) - part_old(i)%scalar_var(12))/(part_new(i)%Pos(3) - part_old(i)%Pos(3))
                   fdyn_tmp = part_old(i)%scalar_var(11) + ( Nz_ind_old*dz-part_old(i)%Pos(3) + (k-Nz_ind_old-1)*dz + 0.5 * dz) *&
                   (part_new(i)%scalar_var(11) - part_old(i)%scalar_var(11))/(part_new(i)%Pos(3) - part_old(i)%Pos(3))
                   w2_tmp = (part_old(i)%Vel(3))**2 + ( (k-1)*dz-part_old(i)%Pos(3)) * &
                   ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                   w2(k) = w2(k) + w2_tmp
                 else
                   dz_part= ztop - (k-1) * dz
                   fbuoy_tmp = part_old(i)%scalar_var(12) + (part_new(i)%Pos(3) - part_old(i)%Pos(3) - 0.5 * dz_part) *&
                   (part_new(i)%scalar_var(12) - part_old(i)%scalar_var(12))/(part_new(i)%Pos(3) - part_old(i)%Pos(3))
                   fdyn_tmp = part_old(i)%scalar_var(11) + (part_new(i)%Pos(3) - part_old(i)%Pos(3) - 0.5 * dz_part) *&
                   (part_new(i)%scalar_var(11) - part_old(i)%scalar_var(11))/(part_new(i)%Pos(3) - part_old(i)%Pos(3))
                   w2_tmp = (part_old(i)%Vel(3))**2 + ( (k-1)*dz-part_old(i)%Pos(3)) * &
                   ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                   w2(k) = w2(k) + w2_tmp
                 end if
                 Fbuoy(k) = Fbuoy(k) +  fbuoy_tmp * dz_part
                 Fdyn(k)  = Fdyn(k)  +  fdyn_tmp  * dz_part
              end do
           end if

        end if


      end do
  end do

  ! get average profiles
  call MPI_Allreduce(N_triggered,local_mpi_int_buffer,Nz+1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  N_triggered  = local_mpi_int_buffer
  call MPI_Allreduce(Fbuoy,local_mpi_real_buffer,Nz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  Fbuoy = local_mpi_real_buffer
  call MPI_Allreduce(Fdyn,local_mpi_real_buffer,Nz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  Fdyn = local_mpi_real_buffer
  call MPI_Allreduce(w2,local_mpi_real_buffer2,Nz+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  w2 = local_mpi_real_buffer2

  Fbuoy = Fbuoy/dfloat(N_triggered(1:Nz))/dz
  Fdyn  = Fdyn /dfloat(N_triggered(1:Nz))/dz
  w2    = w2   /dfloat(N_triggered(1:Nz+1))

  if(root) then
    write(*,*) "Writing Data"

    write(char_dz,'(i16)') int(dz)
    fname_out = trim(adjustl(output_dir))//'Triggered_profile_dz_'//trim(adjustl(char_dz))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"zi",Nz+1,dimids(2)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"zi",nf90_double,dimids(2),zi_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_triggered"       ,nf90_double,dimids(2),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Fdyn"              ,nf90_double,dimids(1),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Fbuoy"             ,nf90_double,dimids(1),var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"w2"                ,nf90_double,dimids(2),var_out_id(4) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, (/(dfloat(i)*dz +.5d0*dz, i=0,Nz-1)/) ) )
    call check_nc( nf90_put_var(output_ncid, zi_varid, (/(dfloat(i)*dz, i=0,Nz)/) ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(N_triggered)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Fdyn ) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Fbuoy) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4), w2   ) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

