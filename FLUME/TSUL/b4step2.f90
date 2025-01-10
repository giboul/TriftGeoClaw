! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux,actualstep)

! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
!
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: g => grav
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: aux_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET

    use storm_module, only: set_storm_fields

    use helpers, only : qa=>q_avac, ts=>times, damping, &
                        inflow_mode, overhang, lake_alt, &
                        fg=>AVAC_fgrid, closest_sup, closest_inf
    ! use helpers, only : fgoutinterp, interp1d4d, AVAC_fgrid, lake_alt
    use fgout_module, only : fgout_grid

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    logical, intent (in) :: actualstep

    ! Local storage
    integer :: index,i,j,k,dummy, inf, sup, w, s
    real(kind=8) :: h,u,v, xw, xe, ys, yn
    real(kind=8) :: xc, yc
    ! real(kind=8) :: qt(size(qa,2),size(qa,3),size(qa,4))
    real(kind=8) :: qt(3,2,2)

    ! ----------------------------------------------------------------
    ! AVAC inflows
    if (trim(inflow_mode) == "src") then
      inf = closest_inf(t, ts)
      sup = closest_sup(t, ts)
      do j = 1, my
        yc = ylower + (j - 0.5d0) * dy
        do i = 1, mx
          xc = xlower + (i - 0.5d0) * dx
          if (lake_alt+overhang<aux(1,i,j)) then
            w = MIN(fg%mx-1, 1+INT((xc-fg%x_low)/(fg%x_hi-fg%x_low)*(fg%mx-1)))
            s = MIN(fg%my-1, 1+INT((yc-fg%y_low)/(fg%y_hi-fg%y_low)*(fg%my-1)))
            if (MAXVAL(qa(inf:sup,1,w:w+1,s:s+1)) > 0) then
              if (ts(inf)<t .and. t<ts(sup) .and. inf<sup) then
                xw = fg%x_low + (w-1)*(fg%x_hi-fg%x_low)/(fg%mx-1)
                ys = fg%y_low + (s-1)*(fg%y_hi-fg%y_low)/(fg%my-1)
                xe = xw + (fg%x_hi-fg%x_low)/(fg%mx-1)
                yn = ys + (fg%y_hi-fg%y_low)/(fg%my-1)
                qt(:,:,:) = qa(inf,1:3,w:w+1,s:s+1) + &
                    (qa(sup,1:3,w:w+1,s:s+1) - &
                     qa(inf,1:3,w:w+1,s:s+1))* &
                    (t-ts(inf)) / (ts(sup)-ts(inf))
                q(1:3,i,j) = (&
                    + (xc-xw) * (yc-ys) * qt(:,2,2) &
                    + (xc-xw) * (yn-yc) * qt(:,2,1) &
                    + (xe-xc) * (yc-ys) * qt(:,1,2) &
                    + (xe-xc) * (yn-yc) * qt(:,1,1) &
                ) / ((xe-xw) * (yn-ys)) * damping
              end if
            end if
          end if
        end do
      end do
    end if
    ! ----------------------------------------------------------------


    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    end forall


    if (aux_finalized < 2 .and. actualstep) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    endif

    ! Set wind and pressure aux variables for this grid
    call set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

end subroutine b4step2
