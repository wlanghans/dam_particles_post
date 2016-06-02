Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  real(8), parameter :: dt = 60.d0, dz=0.
  real(8), parameter :: dA = 4.d4**2.d0 ! z_w2 has to be <=max_height
                                                      ! and > min_height
  integer, parameter :: nzm = 48, Nz_ind_w2=7          ! how many grid level are used zi(nzm+1) is then the max_height,
							! and indicate index of layer above analysis interface
  integer, parameter :: s_dyn =7, s_buoy=8, s_bw=1, s_dw=2, s_w2=3

  real(8), parameter :: Det_timer = 10.d0, dw = 0.01
  real(8), parameter :: trigger_start=600.0, trigger_end=1800., min_height=300.

  logical, parameter :: dojpdf=.false.

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in, int_tmp
  integer :: sum_Np, Nt_tot, N_active, Nz, Nb, Nd  ! choose Nb, Nd such that dbf and dfd are around 0.01
  integer :: N_Count, idx, idx_x, idx_y,ierr, Nz_ind, Nz_ind_old, Nz_ind_top, Nb_ind, Nd_ind
 

  character(128):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_w,char_dw
  integer :: output_ncid
  integer :: dimids(4)

  integer :: z_varid, zi_varid, bw_varid, dw_varid
  integer, dimension(10) :: var_out_id

  integer, dimension(:), allocatable :: N_triggered, local_mpi_int_buffer
  real(8), dimension(:,:), allocatable :: N_bd, w2bd, Buffer_2D
  real(8), dimension(:), allocatable :: Fdyn, Fbuoy, w2, local_mpi_real_buffer, local_mpi_real_buffer2, z, zi, dz_vector

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), time, dz_part, ztop, fdyn_tmp, fbuoy_tmp, w2_tmp, &
             max_bw, max_dw, min_bw, min_dw, max_height, z_w2
  real(8) :: fb_int, fd_int
 

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  allocate(z(nzm))
  allocate(dz_vector(nzm))
  allocate(zi(nzm+1))

  call set_grid(z,zi,dz,dz_vector,nzm)

  max_height = zi(nzm+1)

  if (Nz_ind_w2 .gt.nzm+1) then
    if (root) write(*,*) 'ERROR: Nz_ind_w2 outside of grid'
    call finalize_mpi
  else
    z_w2 = zi(Nz_ind_w2)
  end if

  call Read_Data(1, part_new)

  do i = 1, num_part
    part_new(i)%inactive_time = 0.
    part_new(i)%activity = .false.
    part_new(i)%Vel(1) = 0.
  enddo

  N_active=0
  Nz = nzm 

 
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
  Nz_ind_top=0
  max_bw=0.0
  max_dw=0.0
  min_bw=1.d6
  min_dw=1.d6

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
      if (step_out.eq.N_step.and.part_new(i)%activity.and.part_new(i)%inactive_time.gt.0.0) then 
         N_active=N_active+1
      end if
      if (step_out.eq.N_step) then 
         part_new(i)%Vel(2)=0.0
         part_new(i)%scalar_var(s_bw)=0.0
         part_new(i)%scalar_var(s_dw)=0.0
         part_new(i)%scalar_var(s_w2)=0.0
      end if
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
        ! retain max height reached so far in trajectory
        part_new(i)%Vel(2) = part_old(i)%Vel(2) 
        ! retain work done so far from buoyancy and mechanical forcing
        part_new(i)%Scalar_var(s_bw) = part_old(i)%scalar_var(s_bw) 
        part_new(i)%Scalar_var(s_dw) = part_old(i)%scalar_var(s_dw) 
        ! retain w at height z_w2
        part_new(i)%Scalar_var(s_w2) = part_old(i)%scalar_var(s_w2) 

        
        Nz_ind = minloc( abs(part_new(i)%Pos(3)-zi),1 ) 
        if (zi(Nz_ind).gt.part_new(i)%Pos(3)) Nz_ind=Nz_ind-1
        Nz_ind_old = minloc( abs(part_old(i)%Pos(3)-zi),1 ) 
        if (zi(Nz_ind_old).gt.part_old(i)%Pos(3)) Nz_ind_old=Nz_ind_old-1

        if(part_new(i)%inactive_time.lt.0.0.and.                                 &    ! got triggered
        time.ge.abs(part_new(i)%inactive_time).and.time.le.part_new(i)%Vel(1)) then  ! data point between

        Nz_ind_top = minloc( abs(part_new(i)%Vel(2)-zi),1 ) 
        if (zi(Nz_ind_top).gt.part_new(i)%Vel(2)) Nz_ind_top=Nz_ind_top-1

           ! particle is initiated in layer
           if (time.eq.abs(part_new(i)%inactive_time)) then
              N_triggered( Nz_ind   ) = N_triggered(Nz_ind) + 1
              ! use average for buoy and dyn forcing
              fbuoy_tmp = 0.5 *   (part_new(i)%scalar_var(s_buoy)+part_old(i)%scalar_var(s_buoy))
              fdyn_tmp =  0.5 *   (part_new(i)%scalar_var(s_dyn)-part_old(i)%scalar_var(s_dyn))
              ! part_old(i)%Vel(3)<=0, thus a plus sign below in gradient
               dz_part = -((part_new(i)%Vel(3))**2) * &
              (part_new(i)%Pos(3)-part_old(i)%Pos(3))/ ((part_new(i)%Vel(3))**2+(part_old(i)%Vel(3))**2)
              Fbuoy(Nz_ind)       = Fbuoy(Nz_ind) + dz_part * fbuoy_tmp
              Fdyn(Nz_ind)        = Fdyn(Nz_ind)  + dz_part * fdyn_tmp
              w2(Nz_ind)          = w2(Nz_ind)    + 0.0
              part_new(i)%Scalar_var(s_dw) = part_new(i)%Scalar_var(s_dw) + dz_part * fdyn_tmp
              part_new(i)%Scalar_var(s_bw) = part_old(i)%scalar_var(s_bw) + dz_part * fbuoy_tmp
              part_new(i)%Scalar_var(s_w2) = 0.0
              max_bw = max(max_bw,part_new(i)%Scalar_var(s_bw))
              max_dw = max(max_dw,part_new(i)%Scalar_var(s_dw))
           end if

           if (Nz_ind.ge.Nz_ind_top .and. time.ne.abs(part_new(i)%inactive_time)) then !particle either in top layer or higher

           ! particle is still in same layer, then just add work done in layer
           if (Nz_ind.eq.Nz_ind_top) then
              if (part_new(i)%Pos(3).gt.part_new(i)%Vel(2)) then
               ! again lin interpolate forcings to centerpoint
               ! the clipping here means that im just ignoring any downward segment of the trajectories
                dz_part=part_new(i)%Pos(3)-part_new(i)%Vel(2)                
                Fbuoy(Nz_ind)       = Fbuoy(Nz_ind) + &
                0.5*(part_new(i)%scalar_var(s_buoy)+part_old(i)%scalar_var(s_buoy)) * dz_part
                Fdyn(Nz_ind)        = Fdyn(Nz_ind)  + &
                0.5*(part_new(i)%scalar_var(s_dyn)+part_old(i)%scalar_var(s_dyn)) * dz_part
                part_new(i)%Vel(2) = max(part_new(i)%Pos(3),part_new(i)%Vel(2))
                if (Nz_ind.lt.Nz_ind_w2) then
                  part_new(i)%Scalar_var(s_bw) = part_new(i)%Scalar_var(s_bw) + dz_part * &
                     0.5*(part_new(i)%scalar_var(s_buoy)+part_old(i)%scalar_var(s_buoy))
                  part_new(i)%Scalar_var(s_dw) = part_old(i)%scalar_var(s_dw) + dz_part * &
                     0.5*(part_new(i)%scalar_var(s_dyn)+part_old(i)%scalar_var(s_dyn))
                end if
              end if

           else  !particle rises into layer from below

              if (Nz_ind.gt.Nz) then
                 !particle above or at max_height
                 ztop=max_height
                 if (Nz_ind_top.lt.Nz) then
                   N_triggered( Nz_ind_top+1:Nz  ) = N_triggered(Nz_ind_top+1:Nz) + 1
                 end if
                 N_triggered( Nz+1 ) = N_triggered( Nz+1 ) + 1
                 w2_tmp = (part_old(i)%Vel(3))**2 + ( zi(nzm+1)-part_old(i)%Pos(3)) * &
                 ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                 w2(Nz+1) = w2(Nz+1) + w2_tmp
                 if (Nz_ind_w2.eq.(Nz+1)) then
                    part_new(i)%Scalar_var(s_w2) = w2_tmp
                 end if
              else
                 !particle below max_height
                 ztop= part_new(i)%Pos(3)
                 N_triggered( Nz_ind_top+1:Nz_ind   ) = N_triggered(Nz_ind_top+1:Nz_ind) + 1
              end if

              fbuoy_tmp = 0.5 *(part_new(i)%scalar_var(s_buoy) + part_old(i)%scalar_var(s_buoy))
              fdyn_tmp = 0.5 *(part_new(i)%scalar_var(s_dyn) + part_old(i)%scalar_var(s_dyn))
              !linear interpolate w2 
              do k=Nz_ind_top,min(Nz_ind,Nz)
                 if (k.eq.Nz_ind_top) then
                   dz_part= zi(Nz_ind_top+1)-part_new(i)%Vel(2)
                 elseif (k.gt.Nz_ind_top.and.k.lt.Nz_ind) then
                   dz_part=zi(k+1) - zi(k)
                   w2_tmp = (part_old(i)%Vel(3))**2 + ( zi(k)-part_old(i)%Pos(3)) * &
                   ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                   w2(k) = w2(k) + w2_tmp
                 else
                   dz_part= ztop - zi(k)
                   w2_tmp = (part_old(i)%Vel(3))**2 + ( zi(k)-part_old(i)%Pos(3)) * &
                   ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                   w2(k) = w2(k) + w2_tmp
                 end if
                 Fbuoy(k) = Fbuoy(k) +  fbuoy_tmp * dz_part
                 Fdyn(k)  = Fdyn(k)  +  fdyn_tmp  * dz_part
                 if (k.lt.Nz_ind_w2) then
                  part_new(i)%Scalar_var(s_bw) = part_new(i)%Scalar_var(s_bw) + dz_part * fbuoy_tmp
                  part_new(i)%Scalar_var(s_dw) = part_old(i)%scalar_var(s_dw) + dz_part * fdyn_tmp
                  max_bw = max(max_bw,part_new(i)%Scalar_var(s_bw))
                  max_dw = max(max_dw,part_new(i)%Scalar_var(s_dw))
                 elseif (k.gt.Nz_ind_top.and.k.eq.Nz_ind_w2) then
                  part_new(i)%Scalar_var(s_w2) = w2_tmp
                 end if
              end do
           end if
           end if

        ! store maximum height and index
        part_new(i)%Vel(2) = max(part_new(i)%Pos(3),part_new(i)%Vel(2))

        end if


      end do

      if (step_out.eq.N_step) then
         min_bw = min(min_bw,part_new(i)%Scalar_var(s_bw))
         min_dw = min(min_dw,part_new(i)%Scalar_var(s_dw))
      end if
  end do

   
  if (root) write(*,*) 'Summed vertical profiles obtained'

  ! get average profiles
  call MPI_Allreduce(N_triggered,local_mpi_int_buffer,Nz+1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  N_triggered  = local_mpi_int_buffer
  call MPI_Allreduce(Fbuoy,local_mpi_real_buffer,Nz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  Fbuoy = local_mpi_real_buffer
  call MPI_Allreduce(Fdyn,local_mpi_real_buffer,Nz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  Fdyn = local_mpi_real_buffer
  call MPI_Allreduce(w2,local_mpi_real_buffer2,Nz+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  w2 = local_mpi_real_buffer2

  Fbuoy = Fbuoy/dfloat(N_triggered(1:Nz))/dz_vector(1:Nz)
  Fdyn  = Fdyn /dfloat(N_triggered(1:Nz))/dz_vector(1:Nz)
  w2    = w2   /dfloat(N_triggered(1:Nz+1))

  if (root) write(*,*) 'Average vertical profiles obtained'

  if (dojpdf) then
  if (root) write(*,*) 'Now Working on jpdf'

  ! find max/min buoy and dyn work
  call MPI_Allreduce(max_bw,fbuoy_tmp,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpi_err)
  max_bw = fbuoy_tmp
  call MPI_Allreduce(max_dw,fdyn_tmp,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpi_err)
  max_dw = fdyn_tmp
  call MPI_Allreduce(min_bw,fbuoy_tmp,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpi_err)
  min_bw = fbuoy_tmp
  call MPI_Allreduce(min_dw,fdyn_tmp,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpi_err)
  min_dw = fdyn_tmp

  ! span grid for jpdf
  Nb = min(2,INT((max_bw-min_bw) / dw ) + 1)
  Nd = min(2,INT((max_dw-min_dw) / dw ) + 1)
  max_bw = min_bw + real(Nb) * dw
  max_dw = min_dw + real(Nd) * dw
 
  if (root) then 
    write(*,*) 'Min buoy work:', min_bw
    write(*,*) 'Min dyn  work:', min_dw
    write(*,*) 'Nb:', Nb
    write(*,*) 'Max buoy work:', max_bw
    write(*,*) 'Max dyn  work:', max_dw
    write(*,*) 'Nd:', Nd
  end if

  deallocate(local_mpi_int_buffer)
  deallocate(local_mpi_real_buffer)
  deallocate(local_mpi_real_buffer2)
  allocate(N_bd(Nb,Nd))
  allocate(w2bd(Nb,Nd))
  if (root) write(*,*) 'Allocated N_bd, w2bd'

  N_bd=0.d0
  w2bd=0.d0

  ! now sum up properties of all trajectories of triggered particles
  do i = 1, num_part
  
     if (part_new(i)%inactive_time.lt.0.0) then ! particle gets triggered

        Nb_ind = min(Nb,INT((part_new(i)%Scalar_var(s_bw)-min_bw)/dw) + 1)
        Nd_ind = min(Nd,INT((part_new(i)%Scalar_var(s_dw)-min_dw)/dw) + 1 )
   
        N_bd(Nb_ind,Nd_ind) = N_bd(Nb_ind,Nd_ind) + 1.
        w2bd(Nb_ind,Nd_ind) = w2bd(Nb_ind,Nd_ind) + part_new(i)%Scalar_var(s_w2)
 
     end if
  end do

  end if

  call finalize_particle
  call finalize_core_files

  if (allocated(local_mpi_int_buffer)) deallocate(local_mpi_int_buffer)
  if (allocated(local_mpi_real_buffer)) deallocate(local_mpi_real_buffer)
  if (allocated(local_mpi_real_buffer2)) deallocate(local_mpi_real_buffer2)

  if (dojpdf) then

  ! sum up over tasks
  if (root) write(*,*) 'Reduce_Mtx'
  allocate(Buffer_2D(Nb,Nd))
  call MPI_Allreduce(N_bd,Buffer_2D,Nb*Nd,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  N_bd = Buffer_2D
  call MPI_Allreduce(w2bd,Buffer_2D,Nb*Nd,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  w2bd = Buffer_2D
  deallocate(Buffer_2D)

  ! divide by number in each bin
  if (root) write(*,*) 'Divide by N_bd'
  do i=1,Nb
  do j=1,Nd
    w2bd(i,j) = w2bd(i,j)/(N_bd(i,j)+1.d-6)
  end do
  end do

  if (root) write(*,*) 'jpdf obtained'
  end if !dojpdf



  if(root) then
    write(*,*) "Writing Data"

    write(char_w,'(i16)') int(z_w2)
    write(char_dw,'(i16)') int(1000.*dw)
    if (dojpdf) then
      fname_out = trim(adjustl(output_dir))//'Triggered_profile_budgetz_'//trim(adjustl(char_w))//'_dw_'//trim(adjustl(char_dw))//'.nc'
    else
      fname_out = trim(adjustl(output_dir))//'Triggered_profile_budgetz_'//trim(adjustl(char_w))//'.nc'
    end if
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"zi",Nz+1,dimids(2)) )
    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"zi",nf90_double,dimids(2),zi_varid) )
    if (dojpdf) then
     call check_nc( nf90_def_dim(output_ncid,"bw",Nb,dimids(3)) )
     call check_nc( nf90_def_dim(output_ncid,"dw",Nd,dimids(4)) )
     call check_nc( nf90_def_var(output_ncid,"bw",nf90_double,dimids(3),bw_varid) )
     call check_nc( nf90_def_var(output_ncid,"dw",nf90_double,dimids(4),dw_varid) )
    end if


!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_triggered"       ,nf90_double,dimids(2),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Fdyn"              ,nf90_double,dimids(1),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Fbuoy"             ,nf90_double,dimids(1),var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"w2"                ,nf90_double,dimids(2),var_out_id(4) ) )
    if (dojpdf) then
      call check_nc( nf90_def_var( output_ncid,"w2_bd"             ,nf90_double,dimids(3:4),var_out_id(5) ) )
      call check_nc( nf90_def_var( output_ncid,"N_bd"              ,nf90_double,dimids(3:4),var_out_id(6) ) )
    end if

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, z(1:Nz) ) )
    call check_nc( nf90_put_var(output_ncid, zi_varid, zi(1:Nz+1) ) )
    if (dojpdf) then
      call check_nc( nf90_put_var(output_ncid, bw_varid, (/(min_bw+0.5d0*dw+dfloat(i)*dw, i=1,Nb)/) ) )
      call check_nc( nf90_put_var(output_ncid, dw_varid, (/(min_dw+0.5d0*dw+dfloat(i)*dw, i=1,Nd)/) ) )
    end if

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(N_triggered)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Fdyn ) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Fbuoy) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4), w2   ) )
    if (dojpdf) then
      call check_nc( nf90_put_var(output_ncid, var_out_id(5), w2bd   ) )
      call check_nc( nf90_put_var(output_ncid, var_out_id(6), N_bd   ) )
    end if

    call check_nc( nf90_close(output_ncid) )

  endif

  if (allocated(N_bd)) deallocate(N_bd)
  if (allocated(w2bd)) deallocate(w2bd)
  if (allocated(N_triggered)) deallocate(N_triggered)
  if (allocated(w2)) deallocate(w2)
  if (allocated(Fdyn)) deallocate(Fdyn)
  if (allocated(Fbuoy)) deallocate(Fbuoy)
  if (allocated(z)) deallocate(z)
  if (allocated(zi)) deallocate(zi)
  if (allocated(dz_vector)) deallocate(dz_vector)
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




