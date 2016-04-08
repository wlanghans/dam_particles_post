Program Parallel_Merge_Lagrangian_Data
  use mpi
  use netcdf
  implicit none

  type particle_field
    integer :: g_id !global particle id

    real(8), dimension(3) :: Pos,Vel
    real(8), allocatable, dimension(:) :: scalar_var
    real(8) :: Vterm
    integer :: Category 
    integer :: Natm 
  end type particle_field

  type(particle_field) :: p_info

  character(100) :: param, path_in, fname_in, fname_out, fname_restart, path_out
  character(5) :: rank, rst_str

  integer :: my_id, Nproc, ierr

  integer :: i,j,k, iskip, idesti, idx, rst
  integer :: id_start,id_end, itmp(10), step
  integer :: Np_tot,Np_local,N_infile,N_Local_file, ios,sum_Np, N_restart
  integer :: file_start, file_end
  integer, parameter :: fread_unit = 100, fread_param = 1000

  integer, parameter :: N_Scalar = 11
  character(*), parameter :: restart_name = "Lagrangian_Restart.nc"
  character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qa","qv","qc","qi","qr","qs","qg","t","tabs","p"/)!,"ss","lh","b","conv_12","conv_13","conv_14","conv_x4","conv_x5","conv_x6","conv_y1"/)
  character(*), dimension(N_Scalar), parameter :: var_unit =(/"kg/m^3","kg/kg","kg/kg","kg/kg","kg/kg","kg/kg","kg/kg","kg/kg","K","K","Pa"/)!,"1","J/m^3/s","m/s^2","1/s","1/s","1/s","1/s","1/s","1/s","1/s"/)
  character(*), dimension(N_Scalar), parameter :: var_lname =(/"Density","Dry air","Water vapor","Liquid water","ice","Rain","snow","graupel","Moist Potential Temperature","Temperature","Pressure"/)!,"SuperSaturation","Latent Heat","Buoyancy","Cond to Liquid","Cond to Ice", "Cond to Snow", "Liquid/Ice to Rain","Liquid/Ice to Snow","Liquid/Ice to Graupel","Rain/Snow/Graupel to Vapor"/)

  character(*), dimension(3), parameter :: Coord_var_name =(/"x","y","z"/)
  character(*), dimension(3), parameter :: Coord_var_unit =(/"m","m","m"/)

  character(*), dimension(3), parameter ::   Vel_var_name =(/"u","v","w"/)
  character(*), dimension(4), parameter ::   Vel_var_unit =(/"m/s","m/s","m/s","m/s"/)

  integer :: output_ncid, restart_ncid
  integer ::  gid_dimid, gid_varid,cat_varid, natm_varid
  integer :: time_dimid,time_varid
  integer :: dimids(2), Pos_id(3), Vel_id(4)
  integer :: Restart_Dim_id(2),Restart_Dvar_id(2)
  integer :: Restart_Var_id(N_Scalar+6)
  integer, dimension(N_Scalar) :: Scalar_id 
  integer, allocatable, dimension(:) :: p_ios, Np_frac

  integer, allocatable, dimension(:) :: N_Send, N_Recv, N_Sum_Send, N_Sum_Recv, mpi_id, Send_Counter, mpi_desti

  integer, parameter :: data_type = nf90_real, data_bytes = 4
!  integer, parameter :: data_type = nf90_double, data_bytes = 8

  integer :: Local_Read_in, buffer_size
  real(8), dimension(:,:), allocatable :: Pos,Vel,Scalar_Var
  real(8), dimension(:), allocatable ::  R_Read_Buffer, R_Send_Buffer, R_Recv_Buffer
  integer, dimension(:), allocatable :: Int_Read_Buffer,Int_Read_Buffer2, Int_Read_Buffer3, Int_Send_Buffer,Int_Send_Buffer2, Int_Send_Buffer3, Int_Recv_Buffer, Int_Recv_Buffer2,Int_Recv_Buffer3, Category, Natm
  integer, dimension(:), allocatable :: Send_Request, Recv_Request
  integer, dimension(:,:), allocatable :: mpi_status

  real(8) :: ttime, t_tmp

#if defined(_CRAYFTN) || defined(ibm)
external :: iargc
integer :: iargc
#endif

  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,ierr)

  allocate(N_Send(Nproc), N_Recv(Nproc))
  allocate(N_Sum_Send(Nproc),N_Sum_Recv(Nproc),Send_Counter(Nproc))
  allocate(mpi_id(Nproc-1))
  allocate(Send_Request(Nproc-1),Recv_Request(Nproc-1))
  allocate( mpi_status(MPI_STATUS_SIZE,Nproc-1) )

  itmp(1) = 0
  do i = 0, Nproc-1
    if(i /= my_id) then
      itmp(1) = itmp(1) + 1
      mpi_id( itmp(1) ) = i
    endif
  enddo

  call getarg(1,param)
  path_in = trim(adjustl(param))

  !WL2013
  call getarg(2,param)
  path_out = trim(adjustl(param))
  
  call getarg(3,param)
  N_infile = char2int( trim(adjustl(param)) )

  call getarg(4,param)
  N_restart = char2int( trim(adjustl(param)) )
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Open Data files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  itmp(1) = N_infile/Nproc

  file_start = my_id*itmp(1)
  file_end   = file_start+itmp(1)-1
  if(my_id == Nproc-1) file_end = N_infile-1
  N_local_file = file_end - file_start+1

  rst=1
  write(rst_str,'(i5)') rst
  do i = 1, file_end - file_start + 1
    write(rank,'(i5)') file_start + i - 1
    fname_in = trim(path_in)//'part_'//trim(adjustl(rst_str))//'/LT_data_tmp_'//trim(adjustl(rank))//'.dat'
    open(unit = fread_unit+i, file = fname_in,form='unformatted',status='old', iostat=ios)
    if(ios .gt. 0) then
      write(*,*) 'Error in reading the data file :', fname_in
      call exit()
    endif
    read(fread_unit+i) Np_tot,itmp(1)         !total particle number, number of scalars
  enddo

  Np_local = Np_tot/Nproc
  id_start = Np_local*my_id + 1
  id_end   = id_start -1 + Np_local
  if(my_id == Nproc-1) then
    id_end = Np_tot
    Np_local = id_end - id_start + 1
  endif

  allocate(Pos(Np_local,3),Vel(Np_local,4),Scalar_Var(Np_local,N_scalar),Category(Np_local),Natm(Np_local))

  call MPI_Allreduce(Np_local,sum_Np,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(my_id == 0) write(*,*) Np_tot, sum_Np 
  if(sum_Np /= Np_tot) then
    write(*,*) 'Something is wrong!',Np_tot,sum_Np
    call MPI_FINALIZE(ierr)
  endif


  allocate(p_info%scalar_var(N_scalar))

  allocate(p_ios(N_local_file))
  allocate(Np_frac(N_local_file))

  Local_Read_in = 0
  do i = 1,file_end - file_start + 1
    read(fread_unit+i) itmp(10), ttime        !nstep_in, dfloat(nstep_in)*dt
    read(fread_unit+i) Np_frac(i)             !P_list%Ne (number of particles/processor)
    Local_Read_in = Local_Read_in + Np_frac(i)
  enddo

  call MPI_Allreduce(Local_Read_in,buffer_size,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr) !maximum number of total particles/processor

  iskip = N_Scalar+7
  buffer_size = int(real(buffer_size)*1.25)
  allocate(R_Read_Buffer(buffer_size*iskip),R_Send_Buffer(buffer_size*4),R_Recv_Buffer(buffer_size*4))
  allocate(Int_Read_Buffer(buffer_size),Int_Read_Buffer2(buffer_size),Int_Send_Buffer(buffer_size),Int_Send_Buffer2(buffer_size),Int_Recv_Buffer(buffer_size),Int_Recv_Buffer2(buffer_size),Int_Read_Buffer3(buffer_size),Int_Send_Buffer3(buffer_size),Int_Recv_Buffer3(buffer_size))
  allocate(mpi_desti(buffer_size))


  if(my_id == 0) write(*,*) "Total number of particles:", Np_tot


  step = 1
  do while(.true.)

    call check_nc( nf90_put_var(output_ncid, time_varid, real(ttime,data_bytes), start=(/step/)) )

    call MPI_Allreduce(Local_Read_in,buffer_size,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 
    if (my_id == 0) write(*,*) 'Particle number is ', buffer_size
    if (my_id == 0.and.buffer_size.ne.Np_tot) write(*,*) 'Error in Particle number'

    itmp(1) = 0
    do i = 1, file_end - file_start + 1
      do j = 1,Np_frac(i)
        itmp(1) = itmp(1) + 1
        idx = (itmp(1)-1)*iskip
        read(fread_unit+i) Int_Read_Buffer( itmp(1) ), Int_Read_Buffer2( itmp(1) ), Int_Read_Buffer3( itmp(1) ), R_Read_Buffer(idx+1:idx+iskip)
      enddo !end do over Number of particle
    enddo !end do over N_local_file

    Local_Read_in = 0; p_ios = 0
    do i = 1, file_end - file_start + 1
      read(fread_unit+i,iostat = p_ios(i)) itmp(10), t_tmp
      read(fread_unit+i,iostat = p_ios(i)) Np_frac(i)
      Local_Read_in = Local_Read_in + Np_frac(i)
    enddo

    if(sum(p_ios) /= 0.and.rst.eq.N_restart) then
      exit
    elseif (sum(p_ios) /= 0.and.rst.ne.N_restart) then
      rst=rst+1
      write(rst_str,'(i5)') rst
     
      Local_Read_in = 0
      do i = 1, file_end - file_start + 1
        close(fread_unit+i)
        write(rank,'(i5)') file_start + i - 1
        fname_in =trim(path_in)//'part_'//trim(adjustl(rst_str))//'/LT_data_tmp_'//trim(adjustl(rank))//'.dat'
        open(unit = fread_unit+i, file = fname_in,form='unformatted',status='old',iostat=ios)
        if(ios .gt. 0) then
          write(*,*) 'Error in reading the data file :', fname_in
          call exit()
        endif
        read(fread_unit+i) Np_tot,itmp(1)         !total particle number, number of scalars
        read(fread_unit+i) itmp(10), t_tmp        !nstep_in, dfloat(nstep_in)*dt
        read(fread_unit+i) Np_frac(i)             !P_list%Ne (number of particles/processor)
        Local_Read_in = Local_Read_in + Np_frac(i)
      enddo


    end if

    ttime = t_tmp
    step = step+1
    if(my_id == 0 .and. mod(step,10) == 0) write(*,*) step, ' ' ,ttime

  enddo  !end while

  call MPI_Finalize(ierr)

  if(my_id == 0) write(*,*) 'MPI finalized ...'

  
contains
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


end Program Parallel_Merge_Lagrangian_Data

