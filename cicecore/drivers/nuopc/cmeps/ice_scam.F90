module ice_scam

  use ice_kinds_mod

  implicit none

  ! single column control variables (only used for latlon grid)

  logical :: single_column = .false.    ! true => single column mode
  real (kind=dbl_kind) :: scmlat        ! single column latitude (degrees)
  real (kind=dbl_kind) :: scmlon        ! single column longitude (degrees)
  real (kind=dbl_kind) :: scol_frac     ! single column ocn fraction
  real (kind=dbl_kind) :: scol_mask     ! single column ocn mask
  real (kind=dbl_kind) :: scol_area     ! single column ocn area
  integer              :: scol_ni       ! ni size of single column domain file
  integer              :: scol_nj       ! nj size of single column domain file
  logical              :: scol_valid = .false.    ! true => single column mask is 1

end module ice_scam

