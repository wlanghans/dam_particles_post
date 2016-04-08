Program Particle_Densities
  use particle_data
  use mpi_info
  use core_info
  implicit none

  integer, parameter :: Nz = 50

  real(8), parameter :: dz = 400.d0, dt = 120.d0
  real(8) :: xl, xr, yl, yr

  integer, parameter :: Nx = 40, Ny = 40
  real(8), parameter :: dx = 200.d0, dy = 200.d0

  real(8) :: dA 

  integer :: i,j,k
  integer :: itmp(10), step, step_cnt, step_out, step_in
  integer :: sum_Np, Nt_tot
  integer :: N_Count, idx, idx_x, idx_y,ierr

  character(140):: fname_out
  character(*), parameter :: output_dir = "./"
  character(16) :: char_dz,char_step,char_dA, param_in,char_icat, char_domain,run_name
  integer :: output_ncid
  integer :: dimids(2)

  integer :: t_varid, z_varid, icat
  integer, dimension(3) :: var_out_id

  real(8), allocatable, dimension(:) :: Number_Density, Specific_Number_Density
  integer, allocatable, dimension(:) :: Particle_Number

  real(8) :: cell_volume

  call initialize_mpi
  call initialize_core_files
  call construct_particle

  allocate( Number_Density(Nz) )
  allocate( Particle_Number(Nz) )
  allocate( Specific_Number_Density(Nz) )

  call getarg(3,param_in)
  step = char2int( trim(adjustl(param_in)) )

  call getarg(4,param_in)
  if (root) print*,'Computing statistics for ',trim(adjustl(param_in)),' domain'
  if (trim(adjustl(param_in)).eq.'vsmall' ) then
     xl = 9.5d3
     xr = 1.05d4
     yl = 9.5d3
     yr = 1.05d4
  elseif (trim(adjustl(param_in)).eq.'small' ) then
     xl = 8.0d3
     xr = 1.2d4
     yl = 8.0d3
     yr = 1.2d4
  elseif (trim(adjustl(param_in)).eq.'full' ) then
     xl = 0.0d0
     xr = real(Nx)*dx
     yl = 0.0d0
     yr = real(Ny)*dy
  else
     if (root) print*,'Error: Domain option not supported...'
     call finalize_mpi
  end if
  char_domain=trim(adjustl(param_in))

  call getarg(5,param_in)
  icat = char2int( trim(adjustl(param_in)) )

  call getarg(6,run_name)
    

  dA = (xr-xl) * (yr-yl)

  call Read_Data(step+1, part_new)
  !get layer mass and particle number in that layer;layer masss is stored in Specific_Number_Density
  call compute_particle_mass_subdomain(Nz,Nx,Ny,dx,dy,dz,xl,xr,yl,yr,Specific_Number_Density,Particle_Number,icat)

  !particles per mass
  do k=1,Nz
  if (Specific_Number_Density(k).ne.0.0d0) then
    Specific_Number_Density(k)=dfloat(Particle_Number(k))/Specific_Number_Density(k)
  end if
  end do
  cell_volume = dz*dA
  Number_Density = dfloat(Particle_Number)/cell_volume
  

  if(root) then
    write(*,*) "Writing Data"

    write(char_step,'(i16)') step
    write(char_dA,'(i16)') int(dA)
    write(char_icat,'(i16)') icat
    fname_out = trim(adjustl(output_dir))//'Lagrangian_Densities_'//trim(adjustl(run_name))//'_dA_'//trim(adjustl(char_domain))//'_tstep_'//trim(adjustl(char_step))//'_'//trim(adjustl(char_icat))//'.nc'
    write(*,*) fname_out
    call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

!  Define Dimensions
    call check_nc( nf90_def_dim(output_ncid,"z",Nz,dimids(1)) )
    call check_nc( nf90_def_dim(output_ncid,"time",1,dimids(2)) )

    call check_nc( nf90_def_var(output_ncid,"z",nf90_double,dimids(1),z_varid) )
    call check_nc( nf90_def_var(output_ncid,"time",nf90_double,dimids(2),t_varid) )

!  Define Variables
    call check_nc( nf90_def_var( output_ncid,"N_Particle"             ,nf90_double,dimids,var_out_id(1) ) )
    call check_nc( nf90_def_var( output_ncid,"Number_Density"         ,nf90_double,dimids,var_out_id(2) ) )
!    call check_nc( nf90_def_var( output_ncid,"Specific_Number_Density",nf90_double,dimids,var_out_id(3) ) )
  !  call check_nc( nf90_put_att(output_ncid, var_out_id(3),'FillValue',NF90_FILL_REAL))

    call check_nc( nf90_enddef(output_ncid) )

!  Put variables
    call check_nc( nf90_put_var(output_ncid, z_varid, (/(dfloat(i)*dz + .5d0*dz ,i=0,Nz    -1)/) ) )
    call check_nc( nf90_put_var(output_ncid, t_varid, step*dt ) )

    call check_nc( nf90_put_var(output_ncid, var_out_id(1), dfloat(Particle_Number)) )
    call check_nc( nf90_put_var(output_ncid, var_out_id(2), Number_Density) )
!    call check_nc( nf90_put_var(output_ncid, var_out_id(3), Specific_Number_Density) )

    call check_nc( nf90_close(output_ncid) )

  endif

  call finalize_particle
  call finalize_core_files
  call finalize_mpi

end Program 

