Program Parallel_Merge_Lagrangian_Data
  use mpi
  use netcdf
  implicit none

  type particle_field
    integer :: g_id !global particle id

    real(8), dimension(3) :: Pos,Vel
    real(8), allocatable, dimension(:) :: scalar_var
  end type particle_field

  type(particle_field) :: p_info

  character(100) :: param, path_in, fname_in, fname_out, fname_restart, path_out
  character(5) :: rank

  integer :: my_id, Nproc, ierr

  integer :: i,j,k, iskip, idesti, idx 
  integer :: id_start,id_end, itmp(10), step
  integer :: Np_tot,Np_local,N_infile,N_Local_file, ios,sum_Np
  integer :: file_start, file_end
  integer, parameter :: fread_unit = 100, fread_param = 1000

!  integer, parameter :: N_Scalar = 15
  integer, parameter :: N_Scalar = 11
  character(*), parameter :: restart_name = "Lagrangian_Restart.nc"
  !character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qa","qv","qc","qr","qi","t","tabs","p","b","Fdyn","Fbuoy","Fbuoyd","Fbuoyv","Fbuoyc"/)
  !character(*), dimension(N_Scalar), parameter :: var_unit =(/"kg/m^3","kg/kg","kg/kg","kg/kg","kg/kg","kg/kg","K","K","Pa","m/s2","m/s2","m/s2","m/s2","m/s2","m/s2"/)
  !character(*), dimension(N_Scalar), parameter :: var_lname =(/"Density","Dry air","Water vapor","Liquid water","Rain","ice","Moist Potential Temperature","Temperature","Pressure","Buoyancy","dyn Forcing","buoyant forcing","dry buoy forcing","vapor buoyant forcing","condensate buoyant forcing"/)

!  character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qa","qv","qc","qr","qi","t","tabs","p","b"/)
!  character(*), dimension(N_Scalar), parameter :: var_unit =(/"kg/m^3","kg/kg","kg/kg","kg/kg","kg/kg","kg/kg","K","K","Pa","m/s2"/)
!  character(*), dimension(N_Scalar), parameter :: var_lname =(/"Density","Dry air","Water vapor","Liquid water","Rain","ice","Moist Potential Temperature","Temperature","Pressure","Buoyancy"/)
  character(*), dimension(N_Scalar), parameter :: var_name =(/"rho","qv","qc","qi","t","tabs","Fdyn","Fbuoy","Fbuoyd","Fbuoyv","Fbuoyc"/)
  character(*), dimension(N_Scalar), parameter :: var_unit =(/"kg/m^3","kg/kg","kg/kg","kg/kg","K","K","m/s2","m/s2","m/s2","m/s2","m/s2"/)
  character(*), dimension(N_Scalar), parameter :: var_lname =(/"Density","Water vapor","Liquid water","ice","Moist Potential Temperature","Temperature","Mechanical forcing","buoyant forcing","dry buoyancy","vapor buoyancy","condenate buoyancy"/)

  character(*), dimension(3), parameter :: Coord_var_name =(/"x","y","z"/)
  character(*), dimension(3), parameter :: Coord_var_unit =(/"m","m","m"/)

  character(*), dimension(3), parameter ::   Vel_var_name =(/"u","v","w"/)
  character(*), dimension(3), parameter ::   Vel_var_unit =(/"m/s","m/s","m/s"/)

  integer :: output_ncid, restart_ncid
  integer ::  gid_dimid, gid_varid
  integer :: time_dimid,time_varid
  integer :: dimids(2), Pos_id(3), Vel_id(3)
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
  integer, dimension(:), allocatable :: Int_Read_Buffer,Int_Send_Buffer, Int_Recv_Buffer
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Open Data files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  itmp(1) = N_infile/Nproc

  file_start = my_id*itmp(1)
  file_end   = file_start+itmp(1)-1
  if(my_id == Nproc-1) file_end = N_infile-1
  N_local_file = file_end - file_start+1

  do i = 1, file_end - file_start + 1
    write(rank,'(i5)') file_start + i - 1
    fname_in = trim(path_in)//'LT_data_tmp_'//trim(adjustl(rank))//'.dat'
    open(unit = fread_unit+i, file = fname_in,form='unformatted',status='old', iostat=ios)
    if(ios .gt. 0) then
      write(*,*) 'Error in reading the data file :', fname_in
      call exit()
    endif
    read(fread_unit+i) Np_tot,itmp(1)
  enddo

  Np_local = Np_tot/Nproc
  id_start = Np_local*my_id + 1
  id_end   = id_start -1 + Np_local
  if(my_id == Nproc-1) then
    id_end = Np_tot
    Np_local = id_end - id_start + 1
  endif

  allocate(Pos(Np_local,3),Vel(Np_local,3),Scalar_Var(Np_local,N_scalar))

  call MPI_Allreduce(Np_local,sum_Np,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(my_id == 0) write(*,*) Np_tot, sum_Np 
  if(sum_Np /= Np_tot) then
    write(*,*) 'Something is wrong!',Np_tot,sum_Np
    call MPI_FINALIZE(ierr)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Open Output NCDF files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(rank,'(i5)') my_id
  fname_out = trim(path_out)//'/Lagrangian_data_'//trim(adjustl(rank))//'.nc'

  call check_nc( nf90_create(fname_out    , nf90_noclobber,  output_ncid) )

! Define Coordinate system (gid,time)
  call check_nc( nf90_def_dim(output_ncid, "g_id", Np_local, gid_dimid) )
  call check_nc( nf90_def_dim(output_ncid, "time", nf90_unlimited, time_dimid) )

  call check_nc( nf90_def_var(output_ncid, "g_id", nf90_int , gid_dimid, gid_varid) )
  call check_nc( nf90_put_att(output_ncid, gid_varid, "long_name", "particle global id number") )

  call check_nc( nf90_def_var(output_ncid, "time", data_type,time_dimid,time_varid) )
  call check_nc( nf90_put_att(output_ncid, time_varid, "units", "sec") )
  call check_nc( nf90_put_att(output_ncid, time_varid, "long_name", "time in seconds") )

  dimids = (/gid_dimid,time_dimid/)

! Define Position Vector

  do i = 1,3
    call check_nc( nf90_def_var(output_ncid,Coord_var_name(i),data_type,dimids,Pos_id(i)) )
    call check_nc( nf90_put_att(output_ncid, Pos_id(i), "units", "m") )
  enddo
  call check_nc( nf90_put_att(output_ncid, Pos_id(1), "long_name", "X coordinate - Horizontal axis") )
  call check_nc( nf90_put_att(output_ncid, Pos_id(2), "long_name", "Y coordinate - Horizontal axis") )
  call check_nc( nf90_put_att(output_ncid, Pos_id(3), "long_name", "Z coordinate - Vertical axis") )

! Define Velocity vector
  do i = 1,3
    call check_nc( nf90_def_var(output_ncid, Vel_var_name(i), data_type,dimids,Vel_id(i)) )
    call check_nc( nf90_put_att(output_ncid, Vel_id(i), "units", "m/s") )
  enddo
  call check_nc( nf90_put_att(output_ncid, Vel_id(1), "long_name", "Horizontal velocity in the x-direction") )
  call check_nc( nf90_put_att(output_ncid, Vel_id(2), "long_name", "Horizontal velocity in the y-direction") )
  call check_nc( nf90_put_att(output_ncid, Vel_id(3), "long_name", "Vertical velocity") )

! Define Scalar variables
  do i = 1, N_Scalar
    call check_nc( nf90_def_var(output_ncid, var_name(i), data_type,dimids,Scalar_id(i)) )
    call check_nc( nf90_put_att(output_ncid, Scalar_id(i), "units",var_unit(i)) )
    call check_nc( nf90_put_att(output_ncid, Scalar_id(i), "long_name",var_lname(i)) )
  enddo

  call check_nc( nf90_enddef(output_ncid) )

  do i = 1, Np_local
    call check_nc( nf90_put_var(output_ncid, gid_varid, id_start + i-1 , start=(/i/)) )
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(p_info%scalar_var(N_scalar))

  allocate(p_ios(N_local_file))
  allocate(Np_frac(N_local_file))

  Local_Read_in = 0
  do i = 1,file_end - file_start + 1
    read(fread_unit+i) itmp(10), ttime
    read(fread_unit+i) Np_frac(i)
    Local_Read_in = Local_Read_in + Np_frac(i)
  enddo

  call MPI_Allreduce(Local_Read_in,buffer_size,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

  iskip = N_Scalar+6
  buffer_size = int(real(buffer_size)*1.25)
  allocate(R_Read_Buffer(buffer_size*iskip),R_Send_Buffer(buffer_size*3),R_Recv_Buffer(buffer_size*3))
  allocate(Int_Read_Buffer(buffer_size),Int_Send_Buffer(buffer_size),Int_Recv_Buffer(buffer_size))
  allocate(mpi_desti(buffer_size))


  if(my_id == 0) write(*,*) "Total number of particles:", Np_tot


  step = 1
  do while(.true.)

    call check_nc( nf90_put_var(output_ncid, time_varid, real(ttime,data_bytes), start=(/step/)) )
    itmp(1) = 0
    do i = 1, file_end - file_start + 1
      do j = 1,Np_frac(i)
        itmp(1) = itmp(1) + 1
        idx = (itmp(1)-1)*iskip
        read(fread_unit+i) Int_Read_Buffer( itmp(1) ), R_Read_Buffer(idx+1:idx+iskip)
#if 0
        read(fread_unit+i) p_info%g_id,p_info%Pos(:),p_info%Vel(:),p_info%scalar_var(:)
        Int_Read_Buffer( itmp(1) ) = p_info%g_id
        do k = 1,3
          R_Read_Buffer( idx+k ) = p_info%Pos(k)
        enddo
        do k = 1,3
          R_Read_Buffer( idx+3+k) = p_info%Vel(k)
        enddo
        do k = 1,N_Scalar
          R_Read_Buffer( idx+6+k) = p_info%scalar_var(k)
        enddo
#endif
      enddo !end do over Number of particle
    enddo !end do over N_local_file

    itmp(1) = Np_tot/Nproc
    N_Send = 0; N_Recv = 0; mpi_desti = 0
    do i = 1, Local_Read_in
      idesti = min((Int_Read_Buffer(i)-1)/itmp(1),Nproc-1)
      mpi_desti(i) = idesti
      N_Send(idesti+1) = N_Send(idesti+1) + 1
    enddo

    N_Sum_Send = 0
    do i =2, Nproc
      idx = mpi_id(i-1)+1
      N_Sum_Send(i) = N_Sum_Send(i-1) + N_Send( idx )
    enddo

    N_Recv = 0
    do i = 1, Nproc-1
      call MPI_Irecv(N_Recv( mpi_id(i)+1 ), 1, MPI_INTEGER, mpi_id(i), mpi_id(i), MPI_COMM_WORLD, Recv_Request(i), ierr)
      call MPI_Isend(N_Send( mpi_id(i)+1 ), 1, MPI_INTEGER, mpi_id(i), my_id    , MPI_COMM_WORLD, Send_Request(i), ierr)
    enddo

    Send_Counter = 0
    do i = 1, Local_Read_in
      idesti = mpi_desti(i)
      if(idesti /= my_id) then
        Send_Counter(idesti+1) = Send_Counter(idesti+1) + 1
        if(idesti > my_id) then
          idx = N_Sum_Send(idesti  ) + Send_Counter(idesti+1)
        else
          idx = N_Sum_Send(idesti+1) + Send_Counter(idesti+1)
        endif
        Int_Send_Buffer( idx ) = Int_Read_Buffer(i)
      else
        idx = Int_Read_Buffer(i) - id_start + 1
        do j = 1, 3
          Pos(idx,j) = R_Read_Buffer( (i-1)*iskip + j)
          Vel(idx,j) = R_Read_Buffer( (i-1)*iskip + j + 3)
        enddo
        do j = 1, N_Scalar
          Scalar_Var(idx,j) = R_Read_Buffer( (i-1)*iskip + j + 6)
        enddo
      endif
    enddo

    call MPI_Waitall(Nproc-1, Send_Request, mpi_status,ierr)
    call MPI_Waitall(Nproc-1, Recv_Request, mpi_status,ierr)

    N_Sum_Recv = 0
    do i = 2, Nproc
      idx = mpi_id(i-1)+1
      N_Sum_Recv(i) = N_Sum_Recv(i-1)+N_Recv( idx )
    enddo

    do i = 1, Nproc-1
      idx = N_Sum_Recv(i) 
      itmp(1) = N_Recv( mpi_id(i) + 1 )
      call MPI_Irecv(Int_Recv_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_INTEGER, mpi_id(i), mpi_id(i), MPI_COMM_WORLD, Recv_Request(i), ierr)
      
      idx = N_Sum_Send(i)
      itmp(1) = N_Send( mpi_id(i) + 1 )
      call MPI_Isend(Int_Send_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_INTEGER, mpi_id(i), my_id, MPI_COMM_WORLD, Send_Request(i), ierr)
    enddo

    Send_Counter = 0
    do i = 1, Local_Read_in
      idesti = mpi_desti(i)
      if(idesti /= my_id) then
        Send_Counter(idesti+1) = Send_Counter(idesti+1) + 1
        if(idesti > my_id ) then
          idx = 3*(N_Sum_Send(idesti  ) + Send_Counter(idesti+1) - 1)
        else
          idx = 3*(N_Sum_Send(idesti+1) + Send_Counter(idesti+1) - 1)
        endif
        R_Send_Buffer( idx+1 ) = R_Read_Buffer( (i-1)*iskip+1 )
        R_Send_Buffer( idx+2 ) = R_Read_Buffer( (i-1)*iskip+2 )
        R_Send_Buffer( idx+3 ) = R_Read_Buffer( (i-1)*iskip+3 )
      endif
    enddo

    call MPI_Waitall(Nproc-1, Send_Request, mpi_status,ierr)
    call MPI_Waitall(Nproc-1, Recv_Request, mpi_status,ierr)

    do i = 1, Nproc-1
      idx = 3*N_Sum_Recv(i) 
      itmp(1) = 3*N_Recv( mpi_id(i) + 1 )
      call MPI_Irecv(R_Recv_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), mpi_id(i), MPI_COMM_WORLD, Recv_Request(i), ierr)

      idx = 3*N_Sum_Send(i)
      itmp(1) = 3*N_Send( mpi_id(i) + 1 )
      call MPI_Isend(R_Send_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), my_id, MPI_COMM_WORLD, Send_Request(i), ierr)
    enddo

    call MPI_Waitall(Nproc-1, Send_Request, mpi_status,ierr)

    Send_Counter = 0
    do i = 1, Local_Read_in
      idesti = mpi_desti(i)
      if(idesti /= my_id) then
        Send_Counter(idesti+1) = Send_Counter(idesti+1) + 1
        if(idesti > my_id ) then
          idx = 3*(N_Sum_Send(idesti  ) + Send_Counter(idesti+1) - 1)
        else
          idx = 3*(N_Sum_Send(idesti+1) + Send_Counter(idesti+1) - 1)
        endif
        R_Send_Buffer( idx+1 ) = R_Read_Buffer( (i-1)*iskip+4 )
        R_Send_Buffer( idx+2 ) = R_Read_Buffer( (i-1)*iskip+5 )
        R_Send_Buffer( idx+3 ) = R_Read_Buffer( (i-1)*iskip+6 )
      endif
    enddo

    do i = 1, Nproc-1
      idx = 3*N_Sum_Send(i)
      itmp(1) = 3*N_Send( mpi_id(i) + 1 )
      call MPI_Isend(R_Send_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), my_id, MPI_COMM_WORLD, Send_Request(i), ierr)
    enddo

    call MPI_Waitall(Nproc-1, Recv_Request, mpi_status,ierr)
    do i = 1, N_Sum_Recv(Nproc)
      idx = Int_Recv_Buffer(i) - id_start + 1
      Pos(idx,1) = R_Recv_Buffer( (i-1)*3 + 1)
      Pos(idx,2) = R_Recv_Buffer( (i-1)*3 + 2)
      Pos(idx,3) = R_Recv_Buffer( (i-1)*3 + 3)
    enddo
   
    do i = 1, Nproc-1
      idx = 3*N_Sum_Recv(i) 
      itmp(1) = 3*N_Recv( mpi_id(i) + 1 )
      call MPI_Irecv(R_Recv_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), mpi_id(i), MPI_COMM_WORLD, Recv_Request(i), ierr)
    enddo

    call MPI_Waitall(Nproc-1, Send_Request, mpi_status,ierr)
    Send_Counter = 0
    do i = 1, Local_Read_in
      idesti = mpi_desti(i)
      if(idesti /= my_id) then
        Send_Counter(idesti+1) = Send_Counter(idesti+1) + 1
        if(idesti > my_id ) then
          idx = N_Sum_Send(idesti  ) + Send_Counter(idesti+1) 
        else
          idx = N_Sum_Send(idesti+1) + Send_Counter(idesti+1)
        endif
        R_Send_Buffer( idx ) = R_Read_Buffer( (i-1)*iskip+7 )
      endif
    enddo

    do i = 1, Nproc-1
      idx = N_Sum_Send(i)
      itmp(1) = N_Send( mpi_id(i) + 1 )
      call MPI_Isend(R_Send_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), my_id, MPI_COMM_WORLD, Send_Request(i), ierr)
    enddo

    call MPI_Waitall(Nproc-1, Recv_Request, mpi_status,ierr)
    do i = 1, N_Sum_Recv(Nproc)
      idx = Int_Recv_Buffer(i) - id_start + 1
      Vel(idx,1) = R_Recv_Buffer( (i-1)*3 + 1)
      Vel(idx,2) = R_Recv_Buffer( (i-1)*3 + 2)
      Vel(idx,3) = R_Recv_Buffer( (i-1)*3 + 3)
    enddo
   
    do i = 1, Nproc-1
      idx = N_Sum_Recv(i) 
      itmp(1) = N_Recv( mpi_id(i) + 1 )
      call MPI_Irecv(R_Recv_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), mpi_id(i), MPI_COMM_WORLD, Recv_Request(i), ierr)
    enddo

    call MPI_Waitall(Nproc-1, Send_Request, mpi_status,ierr)
    do k = 1, N_Scalar - 1
      Send_Counter = 0
      do i = 1, Local_Read_in
        idesti = mpi_desti(i)
        if(idesti /= my_id) then
          Send_Counter(idesti+1) = Send_Counter(idesti+1) + 1
          if(idesti > my_id ) then
            idx = N_Sum_Send(idesti  ) + Send_Counter(idesti+1)
          else
            idx = N_Sum_Send(idesti+1) + Send_Counter(idesti+1)
          endif
          R_Send_Buffer( idx ) = R_Read_Buffer( (i-1)*iskip+7+k )
        endif
      enddo

      do i = 1, Nproc-1
        idx = N_Sum_Send(i)
        itmp(1) = N_Send( mpi_id(i) + 1 )
        call MPI_Isend(R_Send_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), my_id, MPI_COMM_WORLD, Send_Request(i), ierr)
      enddo

      call MPI_Waitall(Nproc-1, Recv_Request, mpi_status,ierr)
      do i = 1, N_Sum_Recv(Nproc)
        idx = Int_Recv_Buffer(i) - id_start + 1
        Scalar_Var(idx,k) = R_Recv_Buffer( i )
      enddo
   
      do i = 1, Nproc-1
        idx = N_Sum_Recv(i) 
        itmp(1) = N_Recv( mpi_id(i) + 1 )
        call MPI_Irecv(R_Recv_Buffer(idx+1:idx+itmp(1)), itmp(1), MPI_REAL8, mpi_id(i), mpi_id(i), MPI_COMM_WORLD, Recv_Request(i), ierr)
      enddo

      call MPI_Waitall(Nproc-1, Send_Request, mpi_status,ierr)
    enddo

    call MPI_Waitall(Nproc-1, Recv_Request, mpi_status,ierr)
    do i = 1, N_Sum_Recv(Nproc)
      idx = Int_Recv_Buffer(i) - id_start + 1
      Scalar_Var(idx,N_Scalar) = R_Recv_Buffer( i )
    enddo

    do k = 1,3
      call check_nc( nf90_put_var(output_ncid, Pos_id(k), real(Pos(:,k),data_bytes), start=(/1,step/),count=(/Np_local,1/)) )
      call check_nc( nf90_put_var(output_ncid, Vel_id(k), real(Vel(:,k),data_bytes), start=(/1,step/),count=(/Np_local,1/)) )
    enddo
    do k = 1,N_Scalar
      call check_nc( nf90_put_var(output_ncid, Scalar_id(k), real(Scalar_Var(:,k),data_bytes), start=(/1,step/),count=(/Np_local,1/)) )
    enddo

    Local_Read_in = 0; p_ios = 0
    do i = 1, file_end - file_start + 1
      read(fread_unit+i,iostat = p_ios(i)) itmp(10), t_tmp
      read(fread_unit+i,iostat = p_ios(i)) Np_frac(i)
      Local_Read_in = Local_Read_in + Np_frac(i)
    enddo

    if(sum(p_ios) /= 0) exit

    ttime = t_tmp
    step = step+1
    if(my_id == 0 .and. mod(step,10) == 0) write(*,*) step

  enddo  !end while

#if 0
  if(my_id == 0) then
    Write(*,*) 'Start Writing Restart File'

    N_Recv = 0; N_Sum_Recv = 0
    N_Recv(1:Nproc-2) = Np_tot/Nproc
    N_Recv(  Nproc-1) = Np_tot - (Np_tot/Nproc)*(Nproc-1)

    !fname_restart = trim(path_out)//'data/'//restart_name
    fname_restart = trim(path_out)//restart_name
    call check_nc( nf90_create(fname_restart,nf90_noclobber,restart_ncid) )

   ! Define Coordinate system (gid,time)
    call check_nc( nf90_def_dim(restart_ncid, "g_id", Np_tot, Restart_Dim_id(1)) )
    call check_nc( nf90_def_dim(restart_ncid, "time", 1     , Restart_Dim_id(2)) )

    call check_nc( nf90_def_var(restart_ncid, "g_id", nf90_int , Restart_Dim_id(1), Restart_Dvar_id(1)) )
    call check_nc( nf90_put_att(restart_ncid, Restart_Dvar_id(1), "long_name", "particle global id number") )

    call check_nc( nf90_def_var(restart_ncid, "time", data_type,Restart_Dim_id(2),Restart_Dvar_id(2)) )
    call check_nc( nf90_put_att(restart_ncid, Restart_Dvar_id(2), "units", "sec") )
    call check_nc( nf90_put_att(restart_ncid, Restart_Dvar_id(2), "long_name", "time in seconds") )

    dimids = Restart_Dim_id

  ! Define Position Vector

    do i = 1,3
      call check_nc( nf90_def_var(restart_ncid,Coord_var_name(i),data_type,dimids,Restart_Var_id(i)) )
      call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(i), "units", "m") )
    enddo
    call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(1), "long_name", "X coordinate - Horizontal axis") )
    call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(2), "long_name", "Y coordinate - Horizontal axis") )
    call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(3), "long_name", "Z coordinate - Vertical axis") )

  ! Define Velocity vector
    do i = 1,3
      call check_nc( nf90_def_var(restart_ncid, Vel_var_name(i), data_type,dimids,Restart_Var_id(i+3)) )
      call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(i+3), "units", "m/s") )
    enddo
    call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(4), "long_name", "Horizontal velocity in the x-direction") )
    call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(5), "long_name", "Horizontal velocity in the y-direction") )
    call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(6), "long_name", "Vertical velocity") )

  ! Define Scalar variables
    do i = 1, N_Scalar
      call check_nc( nf90_def_var(restart_ncid, var_name(i), data_type,dimids,Restart_Var_id(i+6)) )
      call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(i+6), "units",var_unit(i)) )
      call check_nc( nf90_put_att(restart_ncid, Restart_Var_id(i+6), "long_name",var_lname(i)) )
    enddo

    call check_nc( nf90_enddef(restart_ncid) )

    call check_nc( nf90_put_var(restart_ncid, Restart_Dvar_id(1), (/(i,i=1,Np_tot)/)) )
    call check_nc( nf90_put_var(restart_ncid, Restart_Dvar_id(2), real(ttime,data_bytes)) )

    do k = 1,3
      call check_nc( nf90_get_var(output_ncid, Pos_id(k), R_Recv_Buffer(1:Np_local), start=(/1,step/),count=(/Np_local,1/)) )
      call check_nc( nf90_put_var(restart_ncid,Restart_Var_id(k), real(R_Recv_Buffer(1:Np_local),data_bytes), start=(/1,1/),count=(/Np_local,1/)) )
    enddo
    do k = 1,3
      call check_nc( nf90_get_var(output_ncid, Vel_id(k), R_Recv_Buffer(1:Np_local), start=(/1,step/),count=(/Np_local,1/)) )
      call check_nc( nf90_put_var(restart_ncid,Restart_Var_id(k+3), real(R_Recv_Buffer(1:Np_local),data_bytes), start=(/1,1/),count=(/Np_local,1/)) )
    enddo
    do k = 1,N_Scalar
      call check_nc( nf90_get_var(output_ncid, Scalar_id(k), R_Recv_Buffer(1:Np_local), start=(/1,step/),count=(/Np_local,1/)) )
      call check_nc( nf90_put_var(restart_ncid,Restart_Var_id(k+6), real(R_Recv_Buffer(1:Np_local),data_bytes), start=(/1,1/),count=(/Np_local,1/)) )
    enddo

  endif

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  do i = 1,Nproc-1
    if(my_id == i) write(*,*) i

    if(my_id == i .or. my_id == 0) then
      do k = 1,3
        if(my_id == i) then
          call check_nc( nf90_get_var(output_ncid,Pos_id(k), R_Send_Buffer(1:Np_local), start=(/1,step/),count=(/Np_local,1/)) )
          call MPI_Send(R_Send_Buffer(1:Np_local), Np_local, MPI_REAL8, 0, i, MPI_COMM_WORLD, Send_Request(1), ierr)
        else
          call MPI_Recv(R_Recv_Buffer(1:N_Recv(i)),N_Recv(i), MPI_REAL8, i, i, MPI_COMM_WORLD, Recv_Request(1), ierr) 
          call check_nc( nf90_put_var(restart_ncid,Restart_Var_id(k), real(R_Recv_Buffer(1:N_Recv(i)),data_bytes), start=(/Np_tot/Nproc*i+1,1/),count=(/N_Recv(i),1/)) )
        endif
        if(my_id == i) then
          call check_nc( nf90_get_var(output_ncid,Vel_id(k), R_Send_Buffer(1:Np_local), start=(/1,step/),count=(/Np_local,1/)) )
          call MPI_Send(R_Send_Buffer(1:Np_local), Np_local, MPI_REAL8, 0, i, MPI_COMM_WORLD, Send_Request(1), ierr)
        else
          call MPI_Recv(R_Recv_Buffer(1:N_Recv(i)),N_Recv(i), MPI_REAL8, i, i, MPI_COMM_WORLD, Recv_Request(1), ierr) 
          call check_nc( nf90_put_var(restart_ncid,Restart_Var_id(k+3), real(R_Recv_Buffer(1:N_Recv(i)),data_bytes), start=(/Np_tot/Nproc*i+1,1/),count=(/N_Recv(i),1/)) )
        endif
      enddo
      do k = 1, N_Scalar
        if(my_id == i) then
          call check_nc( nf90_get_var(output_ncid,Scalar_id(k), R_Send_Buffer(1:Np_local), start=(/1,step/),count=(/Np_local,1/)) )
          call MPI_Send(R_Send_Buffer(1:Np_local), Np_local, MPI_REAL8, 0, i, MPI_COMM_WORLD, Send_Request(1), ierr)
        else
          call MPI_Recv(R_Recv_Buffer(1:N_Recv(i)),N_Recv(i), MPI_REAL8, i, i, MPI_COMM_WORLD, Recv_Request(1), ierr) 
          call check_nc( nf90_put_var(restart_ncid,Restart_Var_id(k+6), real(R_Recv_Buffer(1:N_Recv(i)),data_bytes), start=(/Np_tot/Nproc*i+1,1/),count=(/N_Recv(i),1/)) )
        endif
      enddo
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(my_id == 0) call check_nc( nf90_close(restart_ncid) )
  enddo
#endif

#if 0
  do i = 1, file_end - file_start + 1
    close(fread_unit+i)
  enddo
#endif

  call check_nc( nf90_close(output_ncid) )

  if(my_id == 0) write(*,*) 'NetCDF file closed...'

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

