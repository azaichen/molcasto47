
module h5extractor

contains

subroutine h5molcas(h5file, atoms, charges, coord, overlap, fock,&
  & mo, occ, nbas, prims, primids, basids, density, desym, SYMMETRY)

  ! extract relevant molecule data from Molcas-generated H5 file.  We
  ! use the HDF5 libraries and the Fortran interface (module
  ! hdf5). The assignments of integers in the declarations below as
  ! being of type HID_T vs. platform-native vs. HSIZE_T and
  ! dimension(1) was tested on x86_64 systems running recent versions
  ! of Linux, system-wide installations of HDF5 libraries, and gfortran

  use definitions
  use hdf5

  implicit none

  integer :: ErrorFlag

  character*(LCHARS), intent(IN) :: h5file ! h5file input file

  integer(HID_T) :: file_id

  character(LEN = 1) :: rootname

  integer(HID_T) :: root_id

  character(len=35) :: datatype

  character(len=35) :: data_attribute

  logical :: DMflag, SMflag

  integer(HID_T) :: a_id

  integer(HID_T) :: d_overlap, d_fock, d_density, d_space, d_mo, d_occ ! Datasets

  integer(HID_T) :: d_primitives, d_primids, d_char, d_atoms, d_coord, d_basids

  integer(HID_T) :: a_sid, d_desym

  integer(HID_T), dimension(:), allocatable :: nbas

  real(KREAL), dimension(:), allocatable, intent(OUT) :: prims

  integer, dimension(:), allocatable, intent(OUT) :: primids, basids

  integer(HID_T), dimension(:), allocatable, intent(OUT) :: atoms

  integer(HID_T) :: dmat

  integer(HID_T) :: nprim, nat, nbasids

  integer :: ds_rank, n_ireps

  integer(HSIZE_T), dimension(1) :: ds_dims, maxds_dims, dims, max_dims

  real(KREAL), dimension(:), allocatable, intent(OUT)  :: coord, charges

  real(KREAL), dimension(:), allocatable, intent(OUT)  :: overlap, fock, mo, occ

  real(KREAL), dimension(:), allocatable, intent(OUT)  :: density, desym

  logical, intent(OUT) :: SYMMETRY

  real(KREAL), parameter :: zero=0.0d0

  integer(KINT) :: dbg

  ! ============================================================================

  dbg = 0 ! debug switch

  ! --------------------------------------------
  ! open H5 file and initialize the H5 interface
  ! --------------------------------------------

  write(out,'(/1x,a)') 'Processing H5 file ...'

  call h5open_f(ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error initialising HDF routines'
    stop 'error termination'
  endif

  call h5fopen_f(h5file, H5F_ACC_RDONLY_F, file_id, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) '*** Error opening HDF file'
    stop 'error termination'
  endif

  rootname = '/'

  call h5gopen_f(file_id, rootname, root_id, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error opening root group'
    stop 'error termination'
  endif

  ! -------------------------------------------
  ! retrieve basis set dimension and symmetries
  ! -------------------------------------------

  data_attribute = 'NBAS'

  call h5aopen_name_f(root_id, data_attribute, a_id, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error opening total attribute'
    stop 'error termination'
  endif

  call h5aget_space_f(a_id, a_sid, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error space getting'
    stop 'error termination'
  endif

  call h5sget_simple_extent_dims_f(a_sid, dims, max_dims, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error basis dimensions definition'
    stop 'error termination'
  endif

  n_ireps = INT(max_dims(1),KIND(n_ireps))

  allocate(nbas(n_ireps))

  call h5aread_f(a_id, H5T_STD_I64LE, nbas, dims, ErrorFlag)
  
  write(*,*) nbas

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error basis set size definition'
    stop 'error termination'
  endif

  ! ---------------------------------------------
  ! Extract MCSCF state density matrix if present
  ! ---------------------------------------------

  datatype = 'DENSITY_MATRIX'
  call h5lexists_f(root_id,datatype,DMflag,ErrorFlag)
  if (DMflag) then
    call h5dopen_f(root_id, datatype, d_density, ErrorFlag)
    call h5dget_space_f(d_density, d_space, ErrorFlag)
    call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
    call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
    call h5sget_simple_extent_npoints_f(d_space, dmat, ErrorFlag)
    allocate(density(dmat))
    density = zero
    call h5dread_f(d_density, H5T_IEEE_F64LE, density, ds_dims, ErrorFlag)
    call h5dclose_f(d_density,ErrorFlag)
  else
    write(out,*) ' *** No STATE density matrix found in h5 file.'
  endif

  ! -----------------------------------
  ! Check symmetry and extract the desymmetrization matrix
  ! -----------------------------------

  datatype = 'DESYM_MATRIX'

  call h5lexists_f(root_id,datatype,SMflag,ErrorFlag)

  if (SMflag) then
    SYMMETRY = .true.
    call h5dopen_f(root_id, datatype, d_desym, ErrorFlag)
    call h5dget_space_f(d_desym, d_space, ErrorFlag)
    call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
    call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
    call h5sget_simple_extent_npoints_f(d_space, dmat, ErrorFlag)
    allocate(desym(dmat))
    call h5dread_f(d_desym, H5T_IEEE_F64LE, desym, ds_dims, ErrorFlag)
    call h5dclose_f(d_desym,ErrorFlag)

    write(out,*) ' Calculation appears to have used Symmetry'
  
    ! --------------------------------------
    ! Atomic centers and related information
    ! --------------------------------------

    datatype = 'DESYM_CENTER_ATNUMS'

    call h5dopen_f(root_id, datatype, d_atoms, ErrorFlag)
    
    if (ErrorFlag.lt.0) THEN
         write(*,*) " *** Error geeting atomic numbers"
         return
    endif

    datatype = 'DESYM_CENTER_CHARGES'

    call h5dopen_f(root_id, datatype, d_char, ErrorFlag)
    
    if (ErrorFlag.lt.0) THEN
         write(*,*) " *** Error geeting atomic numbers"
         return
    endif
 

    datatype = 'DESYM_CENTER_COORDINATES'

    call h5dopen_f(root_id, datatype, d_coord, ErrorFlag)

    if (ErrorFlag.lt.0) THEN
        write(*,*) " *** Error geeting atomic coordinates"
        return
    endif

    call h5dget_space_f(d_coord, d_space, ErrorFlag)
    call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
    call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
    call h5sget_simple_extent_npoints_f(d_space, nat, ErrorFlag)


    allocate(coord(nat))
    
    call h5dread_f(d_coord, H5T_IEEE_F64LE, coord, ds_dims, ErrorFlag)


    call h5dget_space_f(d_char, d_space, ErrorFlag)
    call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
    call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
    call h5sget_simple_extent_npoints_f(d_space, nat, ErrorFlag)

    allocate(atoms(nat))
    allocate(charges(nat))

    call h5dread_f(d_char, H5T_IEEE_F64LE, charges, ds_dims, ErrorFlag)
    call h5dread_f(d_atoms, H5T_STD_I64LE, atoms, ds_dims, ErrorFlag)
 
  ! ---------
  ! Basis IDS
  ! ---------

     datatype = 'DESYM_BASIS_FUNCTION_IDS'

     call h5dopen_f(root_id, datatype, d_basids, ErrorFlag)
    
     if (ErrorFlag.lt.0) THEN
         write(*,*) " *** Error geeting basis IDS"
         return
     endif

     call h5dget_space_f(d_basids, d_space, ErrorFlag)
     call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
     call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
     call h5sget_simple_extent_npoints_f(d_space, nbasids, ErrorFlag)

     allocate(basids(nbasids))

     call h5dread_f(d_basids, H5T_NATIVE_INTEGER, basids, ds_dims, ErrorFlag)

  else

    !Nosymmetry case

    write(out,*) ' Symmetry info not found. Assuming C1 symmetry'
    SYMMETRY = .false.
   
     ! --------------------------------------
     ! Atomic centers and related information
     ! --------------------------------------
   
     datatype = 'CENTER_ATNUMS'
   
     call h5dopen_f(root_id, datatype, d_atoms, ErrorFlag)
   
     if (ErrorFlag.lt.0) then
       write(err,*) ' *** Error fetching atomic numbers'
       stop 'error termination'
     endif
   
     datatype = 'CENTER_CHARGES'
   
     call h5dopen_f(root_id, datatype, d_char, ErrorFlag)
   
     if (ErrorFlag.lt.0) then
       write(err,*) ' *** Error geeting atomic numbers'
       stop 'error termination'
     endif
   
     datatype = 'CENTER_COORDINATES'
   
     call h5dopen_f(root_id, datatype, d_coord, ErrorFlag)
   
     if (ErrorFlag.lt.0) then
       write(err,*) ' *** Error geeting atomic coordinates'
       stop 'error termination'
     endif
   
     call h5dget_space_f(d_coord, d_space, ErrorFlag)
     call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
     call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
     call h5sget_simple_extent_npoints_f(d_space, nat, ErrorFlag)
   
     if (dbg>0) write(out,*) 'Coordinates:'
   
     allocate(coord(nat))
   
     call h5dread_f(d_coord, H5T_IEEE_F64LE, coord, ds_dims, ErrorFlag)
   
     if (dbg>0) write(out,'(3F8.4)') coord
   
     ! Atoms:
   
     call h5dget_space_f(d_char, d_space, ErrorFlag)
     call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
     call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
     call h5sget_simple_extent_npoints_f(d_space, nat, ErrorFlag)
   
     if (dbg>0) write(out,*) 'Atoms:'
   
     allocate(atoms(nat))
     allocate(charges(nat))
   
     call h5dread_f(d_char, H5T_IEEE_F64LE, charges, ds_dims, ErrorFlag)
     call h5dread_f(d_atoms, H5T_STD_I64LE, atoms, ds_dims, ErrorFlag)
   
     if (dbg>0) write(out,*) atoms
     if (dbg>0) write(out,*) int(atoms)
   
     ! ---------
     ! Basis IDS
     ! ---------
   
     datatype = 'BASIS_FUNCTION_IDS'
   
     call h5dopen_f(root_id, datatype, d_basids, ErrorFlag)
   
     if (ErrorFlag.lt.0) then
       write(err,*) ' *** Error fetching basis IDS'
       stop 'error termination'
     endif
   
     call h5dget_space_f(d_basids, d_space, ErrorFlag)
     call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
     call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
     call h5sget_simple_extent_npoints_f(d_space, nbasids, ErrorFlag)
   
     if (dbg>0) write(out,*) 'Basis IDS:'
     if (dbg>0) write(out,*) d_space, ds_rank, ds_dims, maxds_dims, nprim
   
     allocate(basids(nbasids))
   
     call h5dread_f(d_basids, H5T_NATIVE_INTEGER, basids, ds_dims, ErrorFlag)
   
     if (dbg>0) write(out,'(4I4)') basids
     
   
  endif !End of the symmetry dependent block
  
  ! --------------------------------------------
  ! Primitives and related information
  ! --------------------------------------------

  datatype = 'PRIMITIVES'

  call h5dopen_f(root_id, datatype, d_primitives, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error geeting primitive basis set'
    stop 'error termination'
  endif

  call h5dget_space_f(d_primitives, d_space, ErrorFlag)
  call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
  call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
  call h5sget_simple_extent_npoints_f(d_space, nprim, ErrorFlag)

  allocate(prims(nprim))

  call h5dread_f(d_primitives, H5T_IEEE_F64LE, prims, ds_dims, ErrorFlag)

  if (dbg>0) write(out,'(2F14.4)') prims

  datatype = 'PRIMITIVE_IDS'

  call h5dopen_f(root_id, datatype, d_primids, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error fetching primitive IDS'
    stop 'error termination'
  endif

  call h5dget_space_f(d_primids, d_space, ErrorFlag)
  call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
  call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
  call h5sget_simple_extent_npoints_f(d_space, nprim, ErrorFlag)

  if (dbg>0) write(out,*) 'Primitives IDS:'
  if (dbg>0) write(out,*) d_space, ds_rank, ds_dims, maxds_dims, nprim

  allocate(primids(nprim))

  call h5dread_f(d_primids, H5T_NATIVE_INTEGER, primids, ds_dims, ErrorFlag)

  if (dbg>0) write(out,'(3I4)') primids

  ! -----------------------------------------
  ! Matrices: MOs, Occupations, Fock, Overlap
  ! -----------------------------------------

  datatype = 'MO_VECTORS'

  call h5dopen_f(root_id, datatype, d_mo, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error fetching MO vector info from h5 file'
    stop 'error termination'
  endif

  datatype = 'AO_FOCKINT_MATRIX'

  call h5dopen_f(root_id, datatype, d_fock, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error fetching Fock matrix info'
    stop 'error termination'
  endif

  datatype = 'AO_OVERLAP_MATRIX'

  call h5dopen_f(root_id, datatype, d_overlap, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error fetching overlap matrix info'
    stop 'error termination'
  endif

  call h5dget_space_f(d_overlap, d_space, ErrorFlag)
  call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
  call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
  call h5sget_simple_extent_npoints_f(d_space, dmat, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error obtaining matrix dimensions!'
    return
  endif

  allocate(overlap(dmat))
  allocate(fock(dmat))
  allocate(mo(dmat))

  call h5dread_f(d_overlap, H5T_IEEE_F64LE, overlap, ds_dims, ErrorFlag)
  call h5dread_f(d_fock, H5T_IEEE_F64LE, fock, ds_dims, ErrorFlag)
  call h5dread_f(d_mo, H5T_IEEE_F64LE, mo, ds_dims, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error fetching matrices from h5 file'
    stop 'error termination'
  endif

  ! Occupation numbers

  datatype = 'MO_OCCUPATIONS'

  call h5dopen_f(root_id, datatype, d_occ, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Error reading MO occupations'
    stop 'error termination'
  endif

  call h5dget_space_f(d_occ, d_space, ErrorFlag)
  call h5sget_simple_extent_ndims_f(d_space, ds_rank, ErrorFlag)
  call h5sget_simple_extent_dims_f(d_space, ds_dims, maxds_dims, ErrorFlag)
  call h5sget_simple_extent_npoints_f(d_space, dmat, ErrorFlag)

  allocate(occ(dmat))

  call h5dread_f(d_occ, H5T_IEEE_F64LE, occ, ds_dims, ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Something is wrong with occupation data in h5 file'
    stop 'error termination'
  endif

  ! --------------------------
  ! Close datasets and h5 file
  ! --------------------------

  call h5dclose_f(d_fock,ErrorFlag)
  call h5dclose_f(d_occ,ErrorFlag)
  call h5dclose_f(d_mo,ErrorFlag)
  call h5aclose_f(a_id, ErrorFlag)
  call h5close_f(ErrorFlag)

  if (ErrorFlag.lt.0) then
    write(err,*) ' *** Something went wrong when closing h5 file'
    stop 'error termination'
  endif

  write(out,'(1x,a/)') '... finished processing H5 file'

  return

end subroutine h5molcas

end module h5extractor
