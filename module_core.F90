module core_info
  use netcdf
  use mpi_info
  implicit none

  character(120) :: param, path
  character(5)  :: rank

  character(120), allocatable, dimension(:) :: fname_in

  integer, parameter :: data_type = nf90_real, data_bytes = 4
!  integer, parameter :: data_type = nf90_double, data_bytes = 8

  integer :: N_infile, N_local_file

  integer, dimension(:), allocatable :: data_ncid, local_count
  integer, dimension(:), allocatable ::  gid_dimid,  gid_varid
  integer, dimension(:), allocatable :: time_dimid, time_varid, Cat_id,Vterm_id, Natm_id
  integer, dimension(:,:), allocatable ::  Pos_id, Vel_id, Scalar_id

  real(8), dimension(:,:), allocatable :: tmp_read

  integer :: N_step, Np_local

  contains

  subroutine initialize_core_files
    use particle_data
    implicit none
    integer :: i, j, id_start, id_end, Np_local

    call getarg(1,path)
    path =trim(adjustl(path))

    call getarg(2,param)
    N_infile = char2int( trim(adjustl(param)) )

    if(N_infile < Nproc) then
      if(root) write(*,*) 'Error: Num. of files is larger than N. Procs.'
      call finalize_mpi
    endif

    N_local_file = N_infile/Nproc
    id_start = N_local_file*my_id
    id_end   = id_start + N_local_file - 1
    if(my_id == Nproc - 1) then
      id_end = N_infile - 1
      N_local_file = id_end - id_start + 1
    endif

    allocate(  fname_in (N_local_file) )
    allocate(  data_ncid(N_local_file) )
    allocate(  gid_dimid(N_local_file), gid_varid(N_local_file) )
    allocate( time_dimid(N_local_file),time_varid(N_local_file) )
    allocate( Pos_id(3,N_local_file), Vel_id(4,N_local_file), Cat_id(N_local_file), Vterm_id(N_local_file), Natm_id(N_local_file) )
    allocate( Scalar_id(N_Scalar,N_local_file) )
    allocate( local_count(N_local_file) )
    
    do i = 1, N_local_file
      write(rank,'(i5)') id_start + (i-1)
      fname_in(i) = trim(path)//'/Lagrangian_data_'//trim(adjustl(rank))//'.nc'
    enddo

    do i = 1, N_local_file
      if (root) write(*,*) 'Opening file ', fname_in(i)
      call check_nc( nf90_open(fname_in(i), nf90_nowrite, data_ncid(i)) )
      call check_nc( nf90_inq_dimid(data_ncid(i),"g_id",  gid_dimid(i)) )
      call check_nc( nf90_inq_dimid(data_ncid(i),"time", time_dimid(i)) )
    enddo

    ! Identify the number of the local particles and total time step
    Np_local = 0; local_count = 0
    do i = 1, N_local_file
      call check_nc( nf90_inquire_dimension(data_ncid(i), gid_dimid(i),len=local_count(i)) )
      call check_nc( nf90_inquire_dimension(data_ncid(i),time_dimid(i),len=N_step   ) )
      Np_local = Np_local+local_count(i)
      N_step=N_step-1 ! subtract initial condition
    enddo

    num_part = Np_local

    allocate ( tmp_read(Np_local,3) )

    ! Get Netcdf_Varid
    do j = 1, N_local_file
      do i = 1, 3
        call check_nc( nf90_inq_varid(data_ncid(j), Coord_var_name(i), Pos_id(i,j)) )
        call check_nc( nf90_inq_varid(data_ncid(j),   Vel_var_name(i), Vel_id(i,j)) )
      enddo
        call check_nc( nf90_inq_varid(data_ncid(j),   Vel_var_name(4), Vel_id(4,j)) )
        call check_nc( nf90_inq_varid(data_ncid(j), "Category" , Cat_id(j)) )
        call check_nc( nf90_inq_varid(data_ncid(j), "N_Atm" , Natm_id(j)) )
      do i = 1, N_Scalar
        call check_nc( nf90_inq_varid(data_ncid(j), var_name(i), Scalar_id(i,j)) )
      enddo
    enddo

  end subroutine initialize_core_files


  subroutine Read_Data(t_loc, p_list)
    use particle_data
    implicit none
    integer, intent(in) :: t_loc
    integer :: i,j,k,idx

    type(particle_field), dimension(:), pointer, intent(inout) :: p_list

    idx = 0
    do k = 1, N_local_file

      do i = 1, 3
        call check_nc( nf90_get_var(data_ncid(k), Pos_id(i,k), tmp_read(:,i), start=(/1,t_loc/), count=(/local_count(k),1/)) )
        do j = 1, local_count(k)
          p_list(idx+j)%Pos(i) = tmp_read(j,i)
        enddo
      enddo

      do i = 1, 3
        call check_nc( nf90_get_var(data_ncid(k), Vel_id(i,k), tmp_read(:,i), start=(/1,t_loc/), count=(/local_count(k),1/)) )
        do j = 1, local_count(k)
          p_list(idx+j)%Vel(i) = tmp_read(j,i)
        enddo
      enddo 

      call check_nc( nf90_get_var(data_ncid(k), Vel_id(4,k), tmp_read(:,1), start=(/1,t_loc/), count=(/local_count(k),1/)) )
      do j = 1, local_count(k)
        p_list(idx+j)%Vterm = tmp_read(j,1)
      enddo

      call check_nc( nf90_get_var(data_ncid(k), Cat_id(k), tmp_read(:,1), start=(/1,t_loc/), count=(/local_count(k),1/)) )
      do j = 1, local_count(k)
        p_list(idx+j)%Category = tmp_read(j,1)
      enddo

      call check_nc( nf90_get_var(data_ncid(k), Natm_id(k), tmp_read(:,1),start=(/1,t_loc/), count=(/local_count(k),1/)) )
      do j = 1, local_count(k)
        p_list(idx+j)%natm = tmp_read(j,1)
      enddo

      do i = 1,N_Scalar
        call check_nc( nf90_get_var(data_ncid(k), Scalar_id(i,k), tmp_read(:,1), start=(/1,t_loc/), count=(/local_count(k),1/)) )
        do j = 1, local_count(k)
          p_list(idx+j)%Scalar_var(i) = tmp_read(j,1)
        enddo
      enddo

      idx = idx+local_count(k)
    enddo
 
  end subroutine Read_Data

  subroutine finalize_core_files
    implicit none
    integer :: i

    do i = 1, N_local_file
      call check_nc( nf90_close(data_ncid(i)) )
    enddo

    deallocate(  fname_in, data_ncid )
    deallocate(  gid_dimid, gid_varid )
    deallocate( time_dimid,time_varid )
    deallocate( Pos_id, Vel_id, Scalar_id, Cat_id, Natm_id )
    deallocate( local_count )

    deallocate( tmp_read )

  end subroutine finalize_core_files

  subroutine check_nc(status)
    integer, intent(in) :: status
    
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      call exit()
    endif
  end subroutine check_nc


  function char2int(c)
    implicit none

    character(*), intent(in) :: c

    integer :: char2int
    integer :: i
    character(100) :: c2

    c2 = trim(c)   ! Eliminate prefixed blanks
    char2int = 0
    do i = 1,len(trim(c2))
      char2int = 10*char2int + index('0123456789',c2(i:i)) - 1
    end do

  end function char2int

end module core_info
