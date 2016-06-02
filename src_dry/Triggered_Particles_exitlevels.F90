Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  real(8), parameter :: dt = 60.d0, dz=0., Lv=2.5d6, cp=1005., grav= 9.81
  real(8), parameter :: dA = 4.d4**2.d0 ! z_w2 has to be <=max_height
                                                      ! and > min_height
  integer , parameter :: nlayer = 1, Nz_ind_w2=7, s_tabs=6, s_qv=2, s_dyn=7,s_buoy=8, &
                    s_qc=3, s_qi=4, s_dw= 1, s_bw=5, s_w2 = 9
!  integer, dimension(nlayer), parameter :: nzm = (/28, 33, 42, 48/)
  integer, dimension(nlayer), parameter :: nzm = (/48/)
  logical, parameter :: dojpdf=.false.

  real(8), parameter :: Det_timer = 10.d0, dw = 0.01
  real(8), parameter :: trigger_start=600.0, trigger_end=1800., min_height=300.

  integer :: i,j,k, k_ind
  integer :: itmp(10), step, step_cnt, step_out, step_in, int_tmp
  integer :: sum_Np, Nt_tot, N_active, Nz, Nb, Nd  ! choose Nb, Nd such that dbf and dfd are around 0.01
  integer :: idx, idx_x, idx_y,ierr, Nz_ind, Nz_ind_old, Nz_ind_top, Nb_ind, Nd_ind
 

  character(128):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_w,char_dw
  integer :: output_ncid
  integer :: dimids(4)

  integer :: z_varid, zi_varid, bw_varid, dw_varid, ze_varid
  integer, dimension(10) :: var_out_id

  integer, dimension(:,:), allocatable :: N_triggered
  integer, dimension(:), allocatable :: local_mpi_int_buffer
  real(8), dimension(:,:,:), allocatable :: N_bd
  real(8), dimension(:,:), allocatable :: bd_tmp
  real(8), dimension(:), allocatable ::  local_mpi_real_buffer, local_mpi_real_buffer2, z, zi, dz_vector
  real(8), dimension(:,:), allocatable :: Fdyn, Fbuoy, w2, mse
  real(8), dimension(:,:,:), allocatable :: w2bd

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), time, dz_part, ztop, fdyn_tmp, fbuoy_tmp, w2_tmp, &
             max_bw, max_dw, min_bw, min_dw, max_height(nlayer)
  real(8) :: fb_int, fd_int, z_w2
 

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  allocate(z(nzm(nlayer)))
  allocate(dz_vector(nzm(nlayer)))
  allocate(zi(nzm(nlayer)+1))

  call set_grid(z,zi,dz,dz_vector,nzm(nlayer))

  do i=1,nlayer
    max_height(i) = zi(nzm(i)+1)
    if (root) write(*,*) 'z= ',max_height(i)
  end do

  z_w2 = zi(Nz_ind_w2)
 

  call Read_Data(1, part_new)

  do i = 1, num_part
    part_new(i)%inactive_time = 0.
    part_new(i)%activity = .false.
    part_new(i)%Vel(1) = 0.
  enddo

  N_active=0
  Nz = nzm(nlayer) 

 
  if (root) write(*,*) 'Using nz = ',Nz,' vertical grid layers'

  allocate(N_triggered(nlayer,nzm(nlayer)+1))
  allocate(w2(nlayer,nzm(nlayer)+1))
  allocate(Fdyn(nlayer,nzm(nlayer)))
  allocate(Fbuoy(nlayer,nzm(nlayer)))
  allocate(mse(nlayer,nzm(nlayer)))
  allocate(local_mpi_int_buffer(nzm(nlayer)+1))
  allocate(local_mpi_real_buffer(nzm(nlayer)))
  allocate(local_mpi_real_buffer2(nzm(nlayer)+1))
  Fdyn  = 0.
  Fbuoy = 0.
  mse = 0.
  N_triggered=0
  w2=0.
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
      if (part_new(i)%inactive_time.eq.0.0 &                              ! make sure its in "not initated" mode
         .and.(part_old(i)%Vel(3).le.0.0.and.part_new(i)%Vel(3).gt.0.0).and. & ! w becomes positive
         (time.ge.trigger_start.and.time.le.trigger_end).and.           & ! within triggering time window   
         part_new(i)%Pos(3).lt.min_height) then                           ! below min_height

         part_new(i)%inactive_time=time

      end if

      do k=1,nlayer
      !check if formerly active particle is considered as triggered
      if (part_old(i)%activity.and.part_new(i)%Pos(3).ge.max_height(k)) then
           ! swap sign to negative to indicate triggered particle
           if (part_new(i)%inactive_time.gt.0.0) part_new(i)%inactive_time=-part_new(i)%inactive_time
           ! save triggering height
           part_new(i)%Vel(1) = part_new(i)%Pos(3)
      end if
      end do
     
      if (part_new(i)%inactive_time.lt.0.0) then ! triggered particle
          if (part_old(i)%activity.and.part_new(i)%Vel(3).gt.0.0.and.&
            (part_new(i)%scalar_var(s_qc)+part_new(i)%scalar_var(s_qi)  ).gt.1.d-5) then
            part_new(i)%activity=.true. ! stays entrained
          else
            part_new(i)%activity=.false.! stays detrained or detrains
          end if
      elseif (part_new(i)%inactive_time.gt.0.0) then ! initiated but not triggered yet
          if (part_new(i)%Vel(3).gt.0.0.and.&
            (part_new(i)%scalar_var(s_qc)+part_new(i)%scalar_var(s_qi)  ).gt.1.d-5) then
            part_new(i)%activity=.true. ! entrained
          elseif (part_new(i)%Vel(3).le.0.0) then
            part_new(i)%inactive_time=0.  ! set particle back to "not initiated" particle
            part_new(i)%activity=.false.  ! and unactive
          end if
      end if

      
      ! count triggered and not yet detrained particles at end of simulation
      if (step_out.eq.N_step.and.part_new(i)%activity.and.part_new(i)%inactive_time.lt.0.0) then 
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

  call MPI_Allreduce(N_active,int_tmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      N_active = int_tmp
  if (root) write(*,*) N_active ,' particles are triggered and not yet detrained at end of simulation'


  ! now sum up properties of all trajectories of triggered particles
  do step_out=2,N_step
    time=step_out*dt
    part_old = part_new
    call Read_Data(step_out, part_new)
      do i = 1, num_part

        ! retain triggering height (ie, maximum height)
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
        time.ge.abs(part_new(i)%inactive_time).and.part_new(i)%Vel(2).lt.part_new(i)%Vel(1)) then  ! data point below max trigger height

        Nz_ind_top = minloc( abs(part_new(i)%Vel(2)-zi),1 ) 
        if (zi(Nz_ind_top).gt.part_new(i)%Vel(2)) Nz_ind_top=Nz_ind_top-1

           ! particle is initiated in layer
           if (time.eq.abs(part_new(i)%inactive_time)) then
              do k=nlayer,1,-1
                if (part_new(i)%Vel(1).ge.max_height(k)) then 
                   N_triggered(k, Nz_ind   ) = N_triggered(k,Nz_ind) + 1
                   ! use average for buoy and dyn forcing
                   fbuoy_tmp = 0.5 *   (part_new(i)%scalar_var(s_buoy)+part_old(i)%scalar_var(s_buoy))
                   fdyn_tmp =  0.5 *   (part_new(i)%scalar_var(s_dyn)-part_old(i)%scalar_var(s_dyn))
                   ! part_old(i)%Vel(3)<=0, thus a plus sign below in gradient
                   dz_part = -((part_new(i)%Vel(3))**2) * &
                   (part_new(i)%Pos(3)-part_old(i)%Pos(3))/ ((part_new(i)%Vel(3))**2+(part_old(i)%Vel(3))**2)
                   Fbuoy(k,Nz_ind)       = Fbuoy(k,Nz_ind) + dz_part * fbuoy_tmp
                   Fdyn(k,Nz_ind)        = Fdyn(k,Nz_ind)  + dz_part * fdyn_tmp
                   mse(k,Nz_ind)         = mse(k,Nz_ind) + (part_new(i)%Pos(3)-zi(Nz_ind)) * 0.5 *&
                   ((part_new(i)%scalar_var(s_tabs)+part_new(i)%scalar_var(s_qv)*Lv/cp +grav/cp*part_new(i)%Pos(3) ) +&
                   (part_old(i)%scalar_var(s_tabs)+part_old(i)%scalar_var(s_qv)*Lv/cp +grav/cp*part_old(i)%Pos(3) ))
                   w2(k,Nz_ind)          = w2(k,Nz_ind)    + 0.0
                   part_new(i)%Scalar_var(s_bw) = part_new(i)%Scalar_var(s_bw) + dz_part * fbuoy_tmp
                   part_new(i)%Scalar_var(s_dw) = part_old(i)%scalar_var(s_dw) + dz_part * fdyn_tmp
                   part_new(i)%Scalar_var(s_w2) = 0.0
                   max_bw = max(max_bw,part_new(i)%Scalar_var(s_bw))
                   max_dw = max(max_dw,part_new(i)%Scalar_var(s_dw))
                   exit
                end if
              end do
           end if

           if (Nz_ind.ge.Nz_ind_top .and. time.ne.abs(part_new(i)%inactive_time)) then !particle either in current top layer or higher

           ! particle is still in same layer, then just add work done in layer
           if (Nz_ind.eq.Nz_ind_top) then
              if (part_new(i)%Pos(3).gt.part_new(i)%Vel(2)) then
                do k=nlayer,1,-1
                  if (part_new(i)%Vel(1).ge.max_height(k)) then 
                    ! again lin interpolate forcings to centerpoint
                    ! the clipping here means that im just ignoring any downward segment of the trajectories
                     dz_part=part_new(i)%Pos(3)-part_new(i)%Vel(2)                
                     Fbuoy(k,Nz_ind)       = Fbuoy(k,Nz_ind) + &
                     0.5*(part_new(i)%scalar_var(s_buoy)+part_old(i)%scalar_var(s_buoy)) * dz_part
                     Fdyn(k,Nz_ind)        = Fdyn(k,Nz_ind)  + &
                     0.5*(part_new(i)%scalar_var(s_dyn)+part_old(i)%scalar_var(s_dyn)) * dz_part
                     mse(k,Nz_ind)         = mse(k,Nz_ind) + dz_part * 0.5 *&
                     ((part_new(i)%scalar_var(s_tabs)+part_new(i)%scalar_var(s_qv)*Lv/cp +grav/cp*part_new(i)%Pos(3) ) +&
                     (part_old(i)%scalar_var(s_tabs)+part_old(i)%scalar_var(s_qv)*Lv/cp +grav/cp*part_old(i)%Pos(3) ))
                     part_new(i)%Vel(2) = max(part_new(i)%Pos(3),part_new(i)%Vel(2))
                     if (Nz_ind.lt.Nz_ind_w2) then
                       part_new(i)%Scalar_var(s_bw) = part_new(i)%Scalar_var(s_bw) + dz_part * &
                          0.5*(part_new(i)%scalar_var(s_buoy)+part_old(i)%scalar_var(s_buoy))
                       part_new(i)%Scalar_var(s_dw) = part_old(i)%scalar_var(s_dw) + dz_part * &
                          0.5*(part_new(i)%scalar_var(s_dyn)+part_old(i)%scalar_var(s_dyn))
                       max_bw = max(max_bw,part_new(i)%Scalar_var(s_bw))
                       max_dw = max(max_dw,part_new(i)%Scalar_var(s_dw))
                     end if
                     exit
                  end if
                end do
              end if

           else  !particle rises into layer from below

              do k=nlayer,1,-1
                if (part_new(i)%Vel(1).ge.max_height(k)) then 
                   k_ind=k
                   exit
                end if
              end do

              if (zi(Nz_ind).ge.max_height(k_ind)) then
                 !particle above or at max_height
                 ztop=max_height(k_ind)
                 if (Nz_ind_top.lt.nzm(k_ind)) then
                   N_triggered( k_ind,Nz_ind_top+1:nzm(k_ind)  ) = N_triggered(k_ind,Nz_ind_top+1:nzm(k_ind)) + 1
                 end if
                 N_triggered( k_ind,nzm(k_ind)+1 ) = N_triggered( k_ind,nzm(k_ind)+1 ) + 1
                 w2_tmp = (part_old(i)%Vel(3))**2 + ( zi(nzm(k_ind)+1)-part_old(i)%Pos(3)) * &
                 ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                 w2(k_ind,nzm(k_ind)+1) = w2(k_ind,nzm(k_ind)+1) + w2_tmp
                 if (Nz_ind_w2.eq.(nzm(k_ind)+1)) then
                    part_new(i)%Scalar_var(7) = w2_tmp
                 end if
              else
                 !particle below max_height
                 ztop= part_new(i)%Pos(3)
                 N_triggered( k_ind,Nz_ind_top+1:Nz_ind   ) = N_triggered(k_ind,Nz_ind_top+1:Nz_ind) + 1
              end if

              fbuoy_tmp = 0.5 *(part_new(i)%scalar_var(s_buoy) + part_old(i)%scalar_var(s_buoy))
              fdyn_tmp = 0.5 *(part_new(i)%scalar_var(s_dyn) + part_old(i)%scalar_var(s_dyn))
              !linear interpolate w2 
              do k=Nz_ind_top,min(Nz_ind,nzm(k_ind))
                 if (k.eq.Nz_ind_top) then
                   dz_part= zi(Nz_ind_top+1)-part_new(i)%Vel(2)
                 elseif (k.gt.Nz_ind_top.and.k.lt.Nz_ind) then
                   dz_part=zi(k+1) - zi(k)
                   w2_tmp = (part_old(i)%Vel(3))**2 + ( zi(k)-part_old(i)%Pos(3)) * &
                   ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                   w2(k_ind,k) = w2(k_ind,k) + w2_tmp
                 else
                   dz_part= ztop - zi(k)
                   w2_tmp = (part_old(i)%Vel(3))**2 + ( zi(k)-part_old(i)%Pos(3)) * &
                   ((part_new(i)%Vel(3))**2-(part_old(i)%Vel(3))**2)/(part_new(i)%Pos(3)-part_old(i)%Pos(3))
                   w2(k_ind,k) = w2(k_ind,k) + w2_tmp
                 end if
                 Fbuoy(k_ind,k) = Fbuoy(k_ind,k) +  fbuoy_tmp * dz_part
                 Fdyn(k_ind,k)  = Fdyn(k_ind,k)  +  fdyn_tmp  * dz_part
                 mse(k_ind,k)         = mse(k_ind,k) + dz_part * 0.5 *&
                 ((part_new(i)%scalar_var(s_tabs)+part_new(i)%scalar_var(s_qv)*Lv/cp +grav/cp*part_new(i)%Pos(3) ) +&
                 (part_old(i)%scalar_var(s_tabs)+part_old(i)%scalar_var(s_qv)*Lv/cp +grav/cp*part_old(i)%Pos(3) ))
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

        ! store maximum height
        part_new(i)%Vel(2) = max(part_new(i)%Pos(3),part_new(i)%Vel(2))

        end if


      end do

      if (step_out.eq.N_step) then
         min_bw = min(min_bw,part_new(i)%Scalar_var(s_bw))
         min_dw = min(min_dw,part_new(i)%Scalar_var(s_dw))
      end if
  end do

  do k=1,nlayer
  ! get average profiles
  call MPI_Allreduce(N_triggered(k,:),local_mpi_int_buffer,nzm(nlayer)+1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  N_triggered(k,:)  = local_mpi_int_buffer
  call MPI_Allreduce(Fbuoy(k,:),local_mpi_real_buffer,nzm(nlayer),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  Fbuoy(k,:) = local_mpi_real_buffer
  call MPI_Allreduce(Fdyn(k,:),local_mpi_real_buffer,nzm(nlayer),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  Fdyn(k,:) = local_mpi_real_buffer
  call MPI_Allreduce(mse(k,:),local_mpi_real_buffer,nzm(nlayer),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  mse(k,:) = local_mpi_real_buffer
  call MPI_Allreduce(w2(k,:),local_mpi_real_buffer2,nzm(nlayer)+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  w2(k,:) = local_mpi_real_buffer2

  Fbuoy(k,:) =Fbuoy(k,:)/(dfloat(N_triggered(k,1:nzm(nlayer)))+1.d-6)/dz_vector(1:Nz)
  Fdyn(k,:)  = Fdyn(k,:)/(dfloat(N_triggered(k,1:nzm(nlayer)))+1.d-6)/dz_vector(1:Nz)
  mse(k,:)   = mse(k,:)/(dfloat(N_triggered(k,1:nzm(nlayer)))+1.d-6)/dz_vector(1:Nz)
  w2(k,:)    = w2(k,:)   /(dfloat(N_triggered(k,1:nzm(nlayer)+1))+1.d-6)
  end do

  if (dojpdf) then
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
  Nb = INT((max_bw-min_bw) / dw ) + 1
  Nd = INT((max_dw-min_dw) / dw ) + 1
  max_bw = min_bw + real(Nb) * dw
  max_dw = min_dw + real(Nd) * dw

  allocate(N_bd(nlayer,Nb,Nd))
  allocate(bd_tmp(Nb,Nd))
  allocate(w2bd(nlayer,Nb,Nd))

  N_bd=0.
  w2bd=0.0
 

  ! now sum up properties of all trajectories of triggered particles for each exit level
  do i = 1, num_part
  
     if (part_new(i)%inactive_time.lt.0.0) then ! particle gets triggered

        do k=nlayer,1,-1
          if (part_new(i)%Vel(1).ge.max_height(k)) then 
             k_ind=k
             exit
          end if
        end do

        Nb_ind = INT((part_new(i)%Scalar_var(s_bw)-min_bw)/dw) + 1
        Nd_ind = INT((part_new(i)%Scalar_var(s_dw)-min_dw)/dw) + 1
   
        N_bd(k_ind,Nb_ind,Nd_ind) = N_bd(k_ind,Nb_ind,Nd_ind) + 1.
        w2bd(k_ind,Nb_ind,Nd_ind) = w2bd(k_ind,Nb_ind,Nd_ind) + part_new(i)%Scalar_var(s_w2)
 
     end if
  end do

  ! sum up over tasks
  do k=1,nlayer
  call MPI_Allreduce(N_bd(k,:,:),bd_tmp,Nb*Nd,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  N_bd(k,:,:) = bd_tmp
  call MPI_Allreduce(w2bd(k,:,:),bd_tmp,Nb*Nd,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  w2bd(k,:,:) = bd_tmp
  end do

  ! divide by number in each bin
  if (root) write(*,*) 'Divide by N_bd'
  do i=1,Nb
  do j=1,Nd
    w2bd(:,i,j) = w2bd(:,i,j)/(N_bd(:,i,j)+1.d-6)
  end do
  end do

  end if

  if(root) then
    write(*,*) "Writing Data"

    write(char_w,'(i16)') int(z_w2)
    write(char_dw,'(i16)') int(1000.*dw)
    if (dojpdf) then
    fname_out = trim(adjustl(output_dir))//'Triggered_multiprofile_budgetz_'//trim(adjustl(char_w))//'_dw_'//trim(adjustl(char_dw))//'.nc'
    else
    fname_out = trim(adjustl(output_dir))//'Triggered_multiprofile_budgetz_'//trim(adjustl(char_w))//'.nc'
    end if
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"zi",Nz+1,dimids(2)) )
    call check_nc( nf90_def_dim(output_ncid,"ze",nlayer,dimids(5)) )
    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"zi",nf90_double,dimids(2),zi_varid) )
    call check_nc( nf90_def_var(output_ncid,"ze",nf90_double,dimids(5),ze_varid) )
    if (dojpdf) then
    call check_nc( nf90_def_dim(output_ncid,"bw",Nb,dimids(3)) )
    call check_nc( nf90_def_dim(output_ncid,"dw",Nd,dimids(4)) )
    call check_nc( nf90_def_var(output_ncid,"bw",nf90_double,dimids(3),bw_varid) )
    call check_nc( nf90_def_var(output_ncid,"dw",nf90_double,dimids(4),dw_varid) )
    end if

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_triggered"       ,nf90_double,dimids((/5,2/)),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Fdyn"              ,nf90_double,dimids((/5,1/)),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Fbuoy"             ,nf90_double,dimids((/5,1/)),var_out_id(3) ) )
    call check_nc( nf90_def_var( output_ncid,"w2"                ,nf90_double,dimids((/5,2/)),var_out_id(4) ) )
    call check_nc( nf90_def_var( output_ncid,"mse"               ,nf90_double,dimids((/5,1/)),var_out_id(7) ) )
    if (dojpdf) then
    call check_nc( nf90_def_var( output_ncid,"w2_bd"             ,nf90_double,dimids((/5,3,4/)),var_out_id(5) ) )
    call check_nc( nf90_def_var( output_ncid,"N_bd"              ,nf90_double,dimids((/5,3,4/)),var_out_id(6) ) )
    end if

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, z(1:Nz) ) )
    call check_nc( nf90_put_var(output_ncid, zi_varid, zi(1:Nz+1) ) )
    call check_nc( nf90_put_var(output_ncid, ze_varid, max_height ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(N_triggered)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Fdyn ) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Fbuoy) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(4), w2   ) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(7), mse) )
    if (dojpdf) then
    call check_nc( nf90_put_var(output_ncid, var_out_id(5), w2bd   ) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(6), N_bd   ) )
    call check_nc( nf90_put_var(output_ncid, bw_varid, (/(min_bw+0.5*dw+dfloat(i)*dw, i=0,Nb-1)/) ) )
    call check_nc( nf90_put_var(output_ncid, dw_varid, (/(min_dw+0.5*dw+dfloat(i)*dw, i=0,Nd-1)/) ) )
    end if

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




