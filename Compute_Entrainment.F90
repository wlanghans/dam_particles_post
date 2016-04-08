Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: Nz = 50
  integer, parameter :: N_avg = 30
  integer, parameter :: N_fil = 120


  real(8), parameter :: dz = 400.d0, dt = 120.d0
  real(8), parameter :: dA = 51.2d3**2.d0

  real(8), dimension(N_fil) :: Det_timer 

  integer, parameter :: Nx = 200, Ny = 200
  real(8), parameter :: dx = 100.d0, dy = 100.d0, w_cutoff=-100.d0
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

  integer, dimension(5),parameter :: icat_cond=(/3,4,5,6,7/)

  integer :: t_varid, z_varid, fil_varid
  integer, dimension(7) :: var_out_id

  real(8), allocatable, dimension(:,:) :: Cell_Particle, Entrain_upp, Lv_flux, inv_den
  real(8), allocatable, dimension(:,:,:) :: Entrain

  real(8), dimension(3) :: Pos

  real(8) :: pi, dist, r_tmp(10), cell_volume

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  Nt_tot = N_step/N_avg - 1

  do i=1,N_fil
     Det_timer(i) = dfloat(i) * 4.0d0 * dt
  end do

  allocate( Cell_Particle(Nz,Nt_tot) )
  allocate( Entrain(N_fil,Nz,Nt_tot) )
  allocate( Entrain_upp(Nz,Nt_tot) )
  allocate( Lv_flux(Nz,Nt_tot) )
  allocate( inv_den(Nz,Nt_tot) )

  call Read_Data(1, part_new)
  !call compute_particle_mass(Nz,dz,dA)
  call getarg(3,param)
  read(param,*)dmass

  if(root) write(*,*) 'w_cutoff',w_cutoff
  if(root) write(*,*) dmass

  Cell_Particle = 0.d0
  Entrain = 0.d0
  Entrain_upp= 0.d0
  Lv_flux = 0.d0
  inv_den = 0.d0

  N_Count = 0
  do i = 1, num_part
    part_new(i)%inactive_time = maxval(Det_timer) + 10000.d0 
  enddo

  step_cnt = 1
  do step_out = 1, Nt_tot


    do step_in = 1, N_avg
      step_cnt = step_cnt + 1
      part_old = part_new

      call Read_Data(step_cnt, part_new)
      do i = 1, num_part


        if (part_new(i)%natm.ne.part_old(i)%natm) then
           part_new(i)%inactive_time = maxval(Det_timer) + 10000.d0
           r_tmp(1) = part_new(i)%Pos(3) 
           part_old(i)%activity=.False.
           part_old(i)%Vel(3)=0.0d0
           part_old(i)%Vterm=0.0d0
           part_old(i)%Scalar_var(2)=part_new(i)%Scalar_var(2)
           part_old(i)%Scalar_var(1)=part_new(i)%Scalar_var(1)
        else
           part_new(i)%inactive_time = part_old(i)%inactive_time + dt
           r_tmp(1) = .5d0*( part_old(i)%Pos(3)+part_new(i)%Pos(3) )
        end if
        
        if (part_new(i)%Category.ne.8) then

        part_new(i)%activity = get_activity(1.0d-05,w_cutoff, part_new(i)%Vel(3), part_new(i)%scalar_var(4) + part_new(i)%scalar_var(5) )
 
        idx = int( r_tmp(1)/dz ) + 1

        if(idx < Nz+1) then

          if(part_new(i)%activity.and.(part_old(i)%Category.eq.2.or.part_old(i)%Category.eq.8)) then

            if(.not.part_old(i)%activity) &
              Entrain_upp(idx,step_out) = Entrain_upp(idx,step_out) + 1.d0
 
            if (dfloat(step_cnt) * dt .ge. maxval(Det_timer)) then

            do j=1,N_fil
              if(.not.part_old(i)%activity .and. (part_new(i)%inactive_time > Det_timer(j) ) ) then
                  Entrain(j,idx,step_out) = Entrain(j,idx,step_out) + 1.d0
              else
                exit
              end if
            end do

            end if
 
            part_new(i)%inactive_time = 0.d0
            Cell_Particle(idx,step_out) = Cell_Particle(idx,step_out) + 1.d0

            Lv_flux(idx,step_out) = Lv_flux(idx,step_out) + .5d0*( part_new(i)%Vel(3) + part_old(i)%Vel(3) ) &
                                                          + .5d0*( part_new(i)%Vterm + part_old(i)%Vterm ) 
            inv_den(idx,step_out) = inv_den(idx,step_out) + 2.d0/( part_new(i)%Scalar_Var(1)*part_new(i)%Scalar_var(2) &
                                                                  +part_old(i)%Scalar_Var(1)*part_old(i)%Scalar_var(2) )

          endif
        endif
        endif !Cat==8
        if (any(part_new(i)%Category.eq.icat_cond)) part_new(i)%inactive_time = 0.d0
      enddo
    enddo

    if(root) write(*,*) step_out, step_cnt
  enddo  

  call Reduce_Mtx(Cell_Particle)
  call Reduce_Mtx(Lv_flux)
  call Reduce_Mtx(inv_den)
  do j=1,N_fil
    call Reduce_Mtx(Entrain(j,1:Nz,1:Nt_tot))
  end do
  call Reduce_Mtx(Entrain_upp)

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
    call check_nc( nf90_def_dim(output_ncid,"filtertime",N_fil,dimids(3)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(2),t_varid) )
    call check_nc( nf90_def_var(output_ncid,"filtertime",nf90_double,dimids(3),fil_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"        ,nf90_double,dimids(1:2),var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment_fil",  nf90_double,dimids((/3,1,2/)),var_out_id(2) ) )
    call check_nc( nf90_def_var( output_ncid,"Entrainment_upp_bd",nf90_double,dimids(1:2),var_out_id(3) ) )

    call check_nc( nf90_def_var( output_ncid,"Vertical_flux"     ,nf90_double,dimids(1:2),var_out_id(4) ) )
    call check_nc( nf90_def_var( output_ncid,"inv_density"       ,nf90_double,dimids(1:2),var_out_id(5) ) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, (/(dfloat(i)*dt*N_avg +.5d0*dt*N_avg,i=0,Nt_tot-1)/) ) )
    call check_nc( nf90_put_var(output_ncid, fil_varid, Det_timer  ) )

    cell_volume = dz*dA

    Cell_Particle = Cell_Particle/dfloat(N_avg)
    Lv_flux = Lv_flux/dfloat(N_avg)
    inv_den = inv_den/dfloat(N_avg)
    Entrain = Entrain/dfloat(N_avg)
    Entrain_upp = Entrain_upp/dfloat(N_avg)

    Lv_flux = dmass*Lv_flux/ cell_volume
    inv_den = dmass*inv_den/ cell_volume

    Entrain     = dmass*Entrain    /(cell_volume*dt)
    Entrain_upp = dmass*Entrain_upp/(cell_volume*dt)

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), Cell_Particle) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Entrain) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Entrain_upp) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(4), Lv_flux) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(5), inv_den) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

