! :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
!> \callgraph
!! \callergraph
!!  Take a grid patch with mesh widths **hx**,**hy**, of dimensions **nrow** by
!!  **ncol**,  and set the values of any piece of
!!  of the patch which extends outside the physical domain
!!  using the boundary conditions.
!!
!!
!!   Specific to geoclaw:  extrapolates aux(i,j,1) at boundaries
!!   to constant.
!!
!!  ### Standard boundary condition choices for amr2ez in clawpack
!!
!!  At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
!!
!!  mthbc(k) =
!!  * 0  for user-supplied BC's (must be inserted!)
!!  * 1  for zero-order extrapolation
!!  * 2  for periodic boundary conditions
!!  * 3  for solid walls, assuming this can be implemented
!!                   by reflecting the data about the boundary and then
!!                   negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!!                   component of q.
!!  * 4  for sphere bcs (left half maps to right half of same side, and vice versa), as if domain folded in half
!!
!!  The corners of the grid patch are at
!!     (xlo_patch,ylo_patch)  --  lower left corner
!!     (xhi_patch,yhi_patch) --  upper right corner
!!
!!  The physical domain itself is a rectangle bounded by
!!     (xlower,ylower)  -- lower left corner
!!     (xupper,yupper)  -- upper right corner
!!
!   This figure below does not work with doxygen
!   the picture is the following:
!  ____________________________________________________
!
!                _____________________ (xupper,yupper)
!               |                     |
!           ____|____ (xhi_patch,yhi_patch)
!           |   |    |                |
!           |   |    |                |
!           |   |    |                |
!           |___|____|                |
!  (xlo_patch,ylo_patch) |            |
!               |                     |
!               |_____________________|
!    (xlower,ylower)
!  ____________________________________________________
!!
!!
!>  Any cells that lie outside the physical domain are ghost cells whose
!!  values should be set in this routine.  This is tested for by comparing
!!  xlo_patch with xlower to see if values need to be set at the left
!   as in the figure above,
!
!>  and similarly at the other boundaries.
!!  Patches are guaranteed to have at least 1 row of cells filled
!!  with interior values so it is possible to extrapolate.
!!  Fix [trimbd()](@ref trimbd) if you want more than 1 row pre-set.
!!
!!  Make sure the order the boundaries are specified is correct
!!  so that diagonal corner cells are also properly taken care of.
!!
!!  Periodic boundaries are set before calling this routine, so if the
!!  domain is periodic in one direction only you
!!  can safely extrapolate in the other direction.
!!
!!  Don't overwrite ghost cells in periodic directions!
!!
!! \param val data array for solution \f$q \f$ (cover the whole grid **msrc**)
!! \param aux data array for auxiliary variables
!! \param nrow number of cells in *i* direction on this grid
!! \param ncol number of cells in *j* direction on this grid
!! \param meqn number of equations for the system
!! \param naux number of auxiliary variables
!! \param hx spacing (mesh size) in *i* direction
!! \param hy spacing (mesh size) in *j* direction
!! \param level AMR level of this grid
!! \param time setting ghost cell values at time **time**
!! \param xlo_patch left bound of the input grid
!! \param xhi_patch right bound of the input grid
!! \param ylo_patch lower bound of the input grid
!! \param yhi_patch upper bound of the input grid
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

subroutine bc2amr(val, aux, nrow, ncol, meqn, naux, &
                  hx, hy, level, time, &
                  xlo_patch, xhi_patch, ylo_patch, yhi_patch)

    use amr_module, only: mthbc, xlower, ylower, xupper, yupper
    use amr_module, only: xperdom, yperdom, spheredom

    use helpers, only: q_avac, times, damping, interp_time, interp_bc_space, &
        inflow_mode

    implicit none

! Input/Output
    integer, intent(in) :: nrow, ncol, meqn, naux, level
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), intent(in) :: xlo_patch, xhi_patch
    real(kind=8), intent(in) :: ylo_patch, yhi_patch
    real(kind=8), intent(in out) :: val(meqn, nrow, ncol)
    real(kind=8), intent(in out) :: aux(naux, nrow, ncol)

! Local storage
    integer :: i, j, i_xy, it, ibeg, jbeg, nxl, nxr, nyb, nyt
    real(kind=8) :: hxmarg, hymarg, xc, yc
    real(kind=8), &
    dimension(size(q_avac,2),size(q_avac,3),size(q_avac,4)) :: qt

    hxmarg = hx*.01d0
    hymarg = hy*.01d0
! Use periodic boundary condition specialized code only, if only one
! boundary is periodic we still proceed below
    if (xperdom .and. (yperdom .or. spheredom)) then
        return
    end if

! Interpolate AVAC state in time before interpolating in space
    if (inflow_mode=="bc") then
        qt = interp_time(time, times, q_avac)
    end if

! Each check has an initial check to ensure that the boundary is a real
! boundary condition and otherwise skips the code.  Otherwise
!-------------------------------------------------------
! Left boundary:
!-------------------------------------------------------
    if (xlo_patch < xlower - hxmarg) then
! number of grid cells from this patch lying outside physical domain:
        nxl = int((xlower + hxmarg - xlo_patch)/hx)

        select case (mthbc(1))
        case (0) ! User defined boundary condition

            do j = 1, ncol
                yc = ylo_patch + (j - 0.5d0)*hy
                do i = 1, nxl ! Avalanche
                    val(1:3,i,j) = interp_bc_space( &
                        yc, qt(1,:,2), qt(1,:,3:5) &
                    ) * damping
                end do
                do i = 1, nxl ! Constant bathymetry extrapolation
                    aux(:, i, j) = aux(:, nxl + 1, j)
                end do
            end do

        case (1) ! Zero-order extrapolation
            do j = 1, ncol
                do i = 1, nxl
                    aux(:, i, j) = aux(:, nxl + 1, j)
                    val(:, i, j) = val(:, nxl + 1, j)
                end do
            end do

        case (2) ! Periodic boundary condition
            continue

        case (3) ! Wall boundary conditions
            do j = 1, ncol
                do i = 1, nxl
                    aux(:, i, j) = aux(:, 2*nxl + 1 - i, j)
                    val(:, i, j) = val(:, 2*nxl + 1 - i, j)
                end do
            end do
            do j = 1, ncol
                do i = 1, nxl
                    val(2, i, j) = -val(2, i, j) ! negate the normal velocity
                end do
            end do

        case (4) ! Spherical domain
            continue

        case default
            print *, "Invalid boundary condition requested."
            stop
        end select
    end if

!-------------------------------------------------------
! Right boundary:
!-------------------------------------------------------
    if (xhi_patch > xupper + hxmarg) then

! number of grid cells lying outside physical domain:
        nxr = int((xhi_patch - xupper + hxmarg)/hx)
        ibeg = max(nrow - nxr + 1, 1)

        select case (mthbc(2))
        case (0) ! User defined boundary condition

            do j = 1, ncol
                yc = ylo_patch + (j - 0.5d0)*hy
                do i = ibeg, nrow ! Avalanche
                    val(1:3,i,j) = interp_bc_space( &
                        yc, qt(2,:,2), qt(2,:,3:5) &
                    ) * damping
                end do
                do i = ibeg, nrow ! Zero-order extrapolation
                    aux(:, i, j) = aux(:, ibeg - 1, j)
                end do
            end do

        case (1) ! Zero-order extrapolation
            do i = ibeg, nrow
                do j = 1, ncol
                    aux(:, i, j) = aux(:, ibeg - 1, j)
                    val(:, i, j) = val(:, ibeg - 1, j)
                end do
            end do

        case (2) ! Periodic boundary condition
            continue

        case (3) ! Wall boundary conditions
            do i = ibeg, nrow
                do j = 1, ncol
                    aux(:, i, j) = aux(:, 2*ibeg - 1 - i, j)
                    val(:, i, j) = val(:, 2*ibeg - 1 - i, j)
                end do
            end do
            do i = ibeg, nrow
                do j = 1, ncol
                    val(2, i, j) = -val(2, i, j) ! negate the normal velocity
                end do
            end do

        case (4) ! Spherical domain
            continue

        case default
            print *, "Invalid boundary condition requested."
            stop

        end select
    end if

!-------------------------------------------------------
! Bottom boundary:
!-------------------------------------------------------
    if (ylo_patch < ylower - hymarg) then

! number of grid cells lying outside physical domain:
        nyb = int((ylower + hymarg - ylo_patch)/hy)

        select case (mthbc(3))
        case (0) ! User defined boundary condition

            do i = 1, nrow
                xc = xlo_patch + (i - 0.5d0)*hx
                do j = 1, nyb ! Avalanche
                    val(1:3,i,j) = interp_bc_space( &
                        xc, qt(3,:,1), qt(3,:,3:5) &
                    ) * damping
                end do
                do j = 1, nyb ! Zero-order extrapolation
                    aux(:, i, j) = aux(:, i, nyb + 1)
                end do
            end do

        case (1) ! Zero-order extrapolation
            do j = 1, nyb
                do i = 1, nrow
                    aux(:, i, j) = aux(:, i, nyb + 1)
                    val(:, i, j) = val(:, i, nyb + 1)
                end do
            end do

        case (2) ! Periodic boundary condition
            continue

        case (3) ! Wall boundary conditions
            do j = 1, nyb
                do i = 1, nrow
                    aux(:, i, j) = aux(:, i, 2*nyb + 1 - j)
                    val(:, i, j) = val(:, i, 2*nyb + 1 - j)
                end do
            end do
            do j = 1, nyb
                do i = 1, nrow
                    val(3, i, j) = -val(3, i, j) ! negate the normal velocity
                end do
            end do

        case (4) ! Spherical domain
            continue

        case default
            print *, "Invalid boundary condition requested."
            stop

        end select
    end if

!-------------------------------------------------------
! Top boundary:
!-------------------------------------------------------
    if (yhi_patch > yupper + hymarg) then

! number of grid cells lying outside physical domain:
        nyt = int((yhi_patch - yupper + hymarg)/hy)
        jbeg = max(ncol - nyt + 1, 1)

        select case (mthbc(4))
        case (0) ! User defined boundary condition

            do i = 1, nrow
                xc = xlo_patch + (i - 0.5d0)*hx
                do j = jbeg, ncol ! Avalanche
                    val(1:3,i,j) = interp_bc_space( &
                        xc, qt(4,:,1), qt(4,:,3:5) &
                    ) * damping
                end do
                do j = jbeg, ncol ! Constant bathymetry extrapolation
                    aux(:, i, j) = aux(:, i, jbeg - 1)
                end do
            end do

        case (1) ! Zero-order extrapolation
            do j = jbeg, ncol
                do i = 1, nrow
                    aux(:, i, j) = aux(:, i, jbeg - 1)
                    val(:, i, j) = val(:, i, jbeg - 1)
                end do
            end do

        case (2) ! Periodic boundary condition
            continue

        case (3) ! Wall boundary conditions
            do j = jbeg, ncol
                do i = 1, nrow
                    aux(:, i, j) = aux(:, i, 2*jbeg - 1 - j)
                    val(:, i, j) = val(:, i, 2*jbeg - 1 - j)
                end do
            end do
            do j = jbeg, ncol
                do i = 1, nrow
                    val(3, i, j) = -val(3, i, j) ! negate the normal velocity
                end do
            end do

        case (4) ! Spherical domain
            continue

        case default
            print *, "Invalid boundary condition requested."
            stop

        end select
    end if

    close (1)

end subroutine bc2amr

