Program Parallel_Statistics
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: Nz = 200
  real(8), parameter :: dz = 100.d0, dt = 4.d0, dA = 2.d4**2.d0
  real(8), parameter :: t_avg = 3300.d0

  integer, parameter :: Step_Avg = int( (t_avg+1.d-10)/dt )

  integer :: i,j,k
  integer :: itmp(10), step_old, step_new
  integer :: sum_Np, Nt_tot
  integer :: Avg_Count, idx_i, idx_j

  real(8), dimension(Nz,Nz) :: T_Mtx

  real(8) :: pi, dist, r_tmp(10), cell_volume

! Output Netcdf parameters
  integer :: output_ncid
  integer :: dimids(2)

  integer :: dim_varid(2)
  integer, dimension(10) :: var_out_id


  call initialize_mpi
  call initialize_core_files
  call construct_particle

  call Read_Data(1, part_new)
  call compute_particle_mass(Nz,dz,dA)

  if(root) write(*,*) 'Mass:', dmass

  T_Mtx = 0.d0; Avg_Count = 0
  step_old = 1; step_new = step_old + Step_Avg

  do while(step_new <= N_step) 

    call Read_Data(step_old, part_old)
    call Read_Data(step_new, part_new)

    do i = 1, num_part
      idx_i = int(part_old(i)%Pos(3) / dz) + 1
      idx_j = int(part_new(i)%Pos(3) / dz) + 1

      if(idx_i < Nz+1 .and. idx_j < Nz+1 ) &
        T_Mtx(idx_i,idx_j) = T_Mtx(idx_i,idx_j) + 1.d0
    enddo

    Avg_Count = Avg_Count + 1
    step_old = step_old+1; step_new = step_new+1
    if( mod(Avg_Count,10) == 0 .and. my_id ==0) write(*,*) step_old, step_new
  enddo  

  call Reduce_Mtx(T_Mtx)
  T_Mtx = T_Mtx / dfloat(Avg_Count) * dmass

  if(root) then
    write(*,*) "Writing Data"

    write(char_tmp,'(i5)') int(t_avg)
    char_tmp = adjustl(char_tmp)
    fname_out = 'T_Mtx_'//trim(char_tmp)//'_sec.nc'
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z_org" ,Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"z_dest",Nz,dimids(2)) )

    call check_nc( nf90_def_var(output_ncid,"z_org" ,nf90_double,dimids(1),dim_varid(1)) )
    call check_nc( nf90_def_var(output_ncid,"z_dest",nf90_double,dimids(2),dim_varid(2)) )

!  Define Variables
    call check_nc( nf90_def_var(output_ncid,"T_Matrix",nf90_double,dimids,var_out_id(1)) )

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, dim_varid(1), (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz-1)/) ) )
    call check_nc( nf90_put_var(output_ncid, dim_varid(2), (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz-1)/) ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), T_Mtx) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

