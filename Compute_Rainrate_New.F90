Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: N_avg = 30
  integer, parameter :: N_fil = 120


  real(8), parameter :: dt = 120.d0
  real(8), parameter :: dA = 51.2d3**2.d0

  real(8), dimension(N_fil) :: Det_timer 

  real(8), parameter :: w_cutoff=-100.d0, tau_f=50.d0*60.d0

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in,step_in2,step_end
  integer :: sum_Np, Nt_tot
  integer :: N_Count, ierr

  character(128):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_tmp
  integer :: output_ncid
  integer :: dimids(2)

  real(8) :: ent_time

  integer, dimension(5),parameter :: icat_cond=(/3,4,5,6,7/)

  integer :: t_varid, fil_varid
  integer, dimension(2) :: var_out_id

  real(8), allocatable, dimension(:) :: rain_upp, local_mpi_real_buffer
  real(8), allocatable, dimension(:,:) :: rain

  logical :: entrainedonce,entrainfound,incloud

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  Nt_tot = N_step/N_avg - 1

  do i=1,N_fil
     Det_timer(i) = dfloat(i) * 4.0d0 * dt
  end do

  allocate( rain(N_fil,Nt_tot) )
  allocate( rain_upp(Nt_tot) )

  call Read_Data(N_step+1, part_old)
  !call compute_particle_mass(Nz,dz,dA)
  call getarg(3,param)
  read(param,*)dmass

  do i = 1, num_part
    part_old(i)%Vterm = 0.0d0
  end do

  if(root) write(*,*) 'w_cutoff',w_cutoff
  if(root) write(*,*) dmass

  rain = 0.d0
  rain_upp= 0.d0

  N_Count = 0

  step_cnt = N_step
  do step_out = Nt_tot,1,-1


    do step_in = 1, N_avg
      step_cnt = step_cnt - 1
      part_new = part_old

      call Read_Data(step_cnt, part_old)
      do i = 1, num_part


        
        !check all raining particles if they entrain for first time
        if (part_new(i)%Category.eq.1 .and. part_new(i)%Vterm.gt.0.0d0) then
          part_old(i)%activity = get_activity(1.0d-05,w_cutoff,part_old(i)%Vel(3), part_old(i)%scalar_var(4) + part_old(i)%scalar_var(5) )
          if (.not.part_old(i)%activity.and.part_new(i)%activity) then
             part_old(i)%inactive_time=0.0d0
          end if
          if (.not.part_old(i)%activity.and..not.part_new(i)%activity.and..not.any(part_old(i)%Category.eq.icat_cond)) then
            part_old(i)%inactive_time=part_new(i)%inactive_time+dt
          end if
          if (part_old(i)%activity.or.any(part_old(i)%Category.eq.icat_cond).and.part_new(i)%inactive_time.lt.tau_f) then
            part_new(i)%inactive_time=0.0d0 
          end if
        end if
        

        ! determine all new raining particles
        if (part_new(i)%Category.eq.8 .and. part_old(i)%Category.ne.8) then ! raining particle
           rain_upp(step_out) = rain_upp(step_out) + 1.0d0
           part_new(i)%Category=1             
           part_new(i)%Vterm=dfloat(step_cnt + 1) ! store rainout time
           part_new(i)%inactive_time=0.0d0
        end if

            

           

           step_end=step_cnt-INT((maxval(Det_timer)+1.0d-30)/dt)
           if (step_end.ge.2) then

           !search for time of first entrainment
           entrainedonce=.False.
           entrainfound=.False.
           incloud=.False.
           do step_in2=step_cnt-1,step_end,-1
              call Read_Data(step_in2, part_new2)
              part_new2(i)%activity = get_activity(1.0d-05,w_cutoff,part_new2(i)%Vel(3), part_new2(i)%scalar_var(4) + part_new2(i)%scalar_var(5) )
               
              if (.not.part_new2(i)%activity.and.incloud) then
                 if (.not. entrainedonce) then
                   part_new2(i)%inactive_time= dt* dfloat(step_in2)
                   entrainedonce=.True.
                   ent_time= part_new2(i)%inactive_time
                 else
                   part_new2(i)%inactive_time= ent_time - dt* dfloat(step_in2)
                   if (part_new2(i)%inactive_time.ge.tau_f) then
                     part_new2(i)%inactive_time=dt* dfloat(step_cnt)-ent_time
                     entrainfound=.True.
                     exit
                   else
                     ent_time=dt* dfloat(step_in2)
                   end if
                 end if
              end if


              if (.not.part_new2(i)%activity) then
                incloud=.False.
              elseif (part_new2(i)%activity) then
                incloud=.True.
              end if
           end do


           if (entrainfound) then
            do j=1,N_fil 
               if (part_new2(i)%inactive_time.lt.Det_timer(j)) rain(j,step_out) = rain(j,step_out) + 1.0d0
            end do
           end if

           else ! time outside maximum filter time
             rain(1:N_fil,step_out)=0.0d0
           end if
        end if 

        if (part_new(i)%Category.eq.1) part_old(i)%Category=1
        if (part_new(i)%Category.eq.1) part_old(i)%Vterm=part_new(i)%Vterm
      enddo

    enddo

    if(root) write(*,*) step_out, step_cnt
  enddo  

  call Reduce_Mtx(rain)
  call Reduce_Mtx(rain_upp)
  call MPI_Reduce(rain_upp,local_mpi_real_buffer,Nt_tot,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
  rain_upp= local_mpi_real_buffer  

  if(root) then
    write(*,*) "Writing Data"

    write(char_tmp,'(i16)') int(N_Avg* (dt + 1.d-10) )
    fname_out = trim(adjustl(output_dir))//'Lagrangian_Rainrate_'//'Tavg_'//trim(adjustl(char_tmp))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"time",Nt_tot,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"filtertime",N_fil,dimids(2)) )

    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(1),t_varid) )
    call check_nc( nf90_def_var(output_ncid,"filtertime",nf90_double,dimids(2),fil_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"rain_fil",  nf90_double,dimids((/2,1/)),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"rain_upp_bd",nf90_double,dimids(1),var_out_id(2) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt*N_avg +.5d0*dt*N_avg,i=0,Nt_tot-1)/) ) )
    call check_nc( nf90_put_var(output_ncid, fil_varid, Det_timer  ) )

    rain = rain/dfloat(N_avg)/dt / dA * dmass
    rain_upp = rain_upp/dfloat(N_avg)/dt / dA * dmass

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), rain) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), rain_upp) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

