module mpi_info
  use mpi
  implicit none

  integer, parameter :: max_mpi_buffer = 500000

  logical :: root = .false.
  integer :: my_id, Nproc, mpi_err

  real(8) :: Mpi_Buffer(max_mpi_buffer)

  interface Reduce_Mtx
    module procedure Reduce_1D, Reduce_2D
  end interface Reduce_Mtx

  contains

  subroutine initialize_mpi
    implicit none

    call MPI_Init(mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,mpi_err)

    if(my_id == 0) root = .true.

  end subroutine initialize_mpi

  subroutine finalize_mpi

    call MPI_FINALIZE(mpi_err)

  end subroutine finalize_mpi

  subroutine Reduce_1D(Array)
    implicit none
    integer :: i, data_cnt, istart, iend, dim1
    real(8), dimension(:), intent(inout) :: Array

    dim1 = size(Array)
    
    data_cnt = dim1
    if(data_cnt > max_mpi_buffer) then
      if(root) then
        write(*,*) 'Error in MPI. Increase buffer size'
        write(*,*) 'Data Cnt',data_cnt,': Buffer',max_mpi_buffer
      endif
      call finalize_mpi
    endif

    call MPI_Allreduce(Array,Mpi_Buffer,data_cnt,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

    Array(1:dim1) = Mpi_Buffer(1:dim1)

  end subroutine Reduce_1D


  subroutine Reduce_2D(Mtx)
    implicit none
    integer :: i, data_cnt, istart, iend, dim1, dim2
    real(8), dimension(:,:), intent(inout) :: Mtx

    dim1 = size(Mtx,dim=1)
    dim2 = size(Mtx,dim=2)
    
    data_cnt = dim1*dim2
    if(data_cnt > max_mpi_buffer) then
      if(root) then
        write(*,*) 'Error in MPI. Increase buffer size'
        write(*,*) 'Data Cnt',data_cnt,': Buffer',max_mpi_buffer
      endif
      call finalize_mpi
    endif

    call MPI_Allreduce(Mtx,Mpi_Buffer,data_cnt,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

    do i = 1, dim2
      istart = (i-1)*dim1 + 1
      iend   = (i-1)*dim1 + dim1
      Mtx(:,i) = Mpi_Buffer( istart:iend )
    enddo

  end subroutine Reduce_2D

end module mpi_info
