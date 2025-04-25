module helpers

    use IEEE_ARITHMETIC, only : ieee_value, ieee_positive_inf
    use fgout_module, only : fgout_grid, set_fgout, FGOUT_fgrids, module_setup
    implicit none
    save

    integer :: bc_size
    real(kind=8), allocatable :: q_avac(:,:,:,:)
    real(kind=8), allocatable :: times(:)
    real(kind=8) :: damping, overhang
    real(kind=8) :: lake_alt
    character(len=255) :: AVAC_DIR, inflow_mode
    type(fgout_grid) :: AVAC_fgrid

contains


    SUBROUTINE read_data()

        INTEGER :: unit=2

        call opendatafile(unit, "setprob.data")
            READ(unit,*) inflow_mode
            READ(unit,*) damping
            READ(unit,*) lake_alt
            READ(unit,*) overhang
            READ(unit,*) AVAC_DIR
            READ(unit,*) bc_size
        CLOSE(unit)

        IF (TRIM(inflow_mode) == "None") then
            inflow_mode = "bc"
        END IF

    END SUBROUTINE read_data


    subroutine read_times()
        character(len=255) :: fname
        integer :: io, n, i
        integer :: unit

        unit = 2
        fname = "_bc_inflows/times.txt"
        print "(A,A)", "Reading ", trim(fname)
        open(unit, file=fname)
            n = 0
            do 
                read(unit,*,iostat=io)
                if (io /= 0) then
                    exit
                end if
                n = n + 1
            end do
            allocate(times(n))
        rewind(unit)
            do i = 1, n
                read(unit,*) times(i)
            end do
        close(unit)
    end subroutine read_times


    subroutine init_bc()

        character(len=6), dimension(4) :: sides
        character(len=255) :: fname
        integer :: unit, i, mthbc
    
        unit = 2
        sides = [character(len=6) :: "left", "right", "bottom", "top"]

        call read_times()
 
 
        allocate(q_avac(size(times), 4, bc_size, 5))
        print "(A,I10)", "size(q_avac) = ", size(q_avac)
        print "(A,I10)", "times:     size(q_avac, 1) = ", size(q_avac, 1)
        print "(A,I10)", "sides:     size(q_avac, 2) = ", size(q_avac, 2)
        print "(A,I10)", "bc_size:   size(q_avac, 3) = ", size(q_avac, 3)
        print "(A,I10)", "variables: size(q_avac, 4) = ", size(q_avac, 4)
   
        do mthbc = 1, 4
            do i = 1, size(times)
                fname = "_bc_inflows/"//trim(sides(mthbc))//"_"
                write(fname,"(A,I0.4,A4)") trim(fname), i-1, ".npy"
                open(unit, file=fname, status="unknown", access="stream")
                    read(unit) q_avac(i, mthbc, :, :)
                close(unit)
            end do
        end do

    end subroutine init_bc


    INTEGER FUNCTION closest_inf(value, array)

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: array
        REAL(KIND=8), INTENT(IN) :: value

        IF (ALL(value<array)) THEN
            closest_inf = MINLOC(array, DIM=1)
        ELSE
            closest_inf = MINLOC(ABS(array-value), MASK=array<=value, DIM=1)
        END IF

    END FUNCTION closest_inf


    INTEGER FUNCTION closest_sup(value, array)

        REAL(kind=8), dimension(:), INTENT(IN) :: array
        REAL(kind=8), INTENT(IN) :: value

        IF (ALL(array<=value)) THEN
            closest_sup = MAXLOC(array, DIM=1)
        ELSE
            closest_sup = MINLOC(ABS(array-value), mask=array>value, DIM=1)
        END IF

    END FUNCTION closest_sup


    FUNCTION interp_time(t, ts, q)

        REAL(kind=8), INTENT(IN) :: t, ts(:), q(:,:,:,:)
        REAL(kind=8), & 
        dimension(SIZE(q,2),SIZE(q,3),SIZE(q,4)) :: interp_time
        INTEGER :: i, s

        i = closest_inf(t, ts)
        s = closest_sup(t, ts)
        IF (ts(s)<t) THEN
            interp_time = q(i,:,:,:)*0
        ELSE IF (t<ts(i)) THEN
            interp_time = q(s,:,:,:)*0
        ELSE IF (ts(i)==ts(s)) THEN
            interp_time = q(i,:,:,:)
        ELSE
            interp_time = q(i,:,:,:) + &
                        (q(s,:,:,:) - q(i,:,:,:)) * &
                        (t - ts(i)) / (ts(s) - ts(i))
        END IF

    END FUNCTION interp_time


    FUNCTION interp_bc_space(xnew, x, q)

        REAL(kind=8), INTENT(IN) :: xnew, x(:), q(:,:)
        REAL(kind=8), dimension(SIZE(q,2)) :: interp_bc_space
        INTEGER :: i, s

        ! Finding the closest positions
        i = closest_inf(xnew, x)
        s = closest_sup(xnew, x)

        ! Return 0 if out of bounds
        ! Else return 1D interpolation
        IF (x(s)<xnew) THEN
            interp_bc_space = q(i,:)*0
        ELSE IF (xnew<x(i)) THEN
            interp_bc_space = q(i,:)*0
        ELSE IF (x(i)==x(s)) THEN
            interp_bc_space = q(i,:)
        ELSE
            interp_bc_space = q(i,:) + &
                        (q(s,:) - q(i,:)) * &
                        (MAX(x(i),MIN(x(s),xnew)) - x(i)) / &
                        (x(s) - x(i))
        END IF

    END FUNCTION interp_bc_space


    FUNCTION interp_src(xc, yc, fg, alpha, q)

        REAL(kind=8), intent(in) :: xc, yc, alpha, q(:,:,:,:)
        type(fgout_grid), intent(in) :: fg

        REAL(kind=8), dimension(3) :: interp_src

        INTEGER :: w, s
        REAL(KIND=8), dimension(3,2,2) :: qt
        REAL(kind=8) :: xw, xe, ys, yn
 
        ! finding the indices of position that's closest to (xc, yc)
        ! 'w' for west, 'w+1' for east,
        ! 's' for south and 's+1' for north
        w = MIN(fg%mx-1, 1+INT((xc-fg%x_low)/(fg%x_hi-fg%x_low)*(fg%mx-1)))
        s = MIN(fg%my-1, 1+INT((yc-fg%y_low)/(fg%y_hi-fg%y_low)*(fg%my-1)))

        ! computing the positions of the 4 closest points
        xw = fg%x_low + (w-1)*(fg%x_hi-fg%x_low)/(fg%mx-1)
        ys = fg%y_low + (s-1)*(fg%y_hi-fg%y_low)/(fg%my-1)
        xe = xw + (fg%x_hi-fg%x_low)/(fg%mx-1)
        yn = ys + (fg%y_hi-fg%y_low)/(fg%my-1)

        qt(:,:,:) = q(1,:,w:w+1,s:s+1) + alpha * &
            ( q(2,:,w:w+1,s:s+1) - q(1,:,w:w+1,s:s+1) )

        ! Regular grid interpolation
        interp_src = ( &
            + (xc-xw) * (yc-ys) * qt(:,2,2) &
            + (xc-xw) * (yn-yc) * qt(:,2,1) &
            + (xe-xc) * (yc-ys) * qt(:,1,2) &
            + (xe-xc) * (yn-yc) * qt(:,1,1) &
        ) / ((xe-xw) * (yn-ys))

    END FUNCTION interp_src


    SUBROUTINE init_src_fgout_bin()

        CHARACTER(len=255) :: file, ftemp
        INTEGER :: i

        ftemp = trim(AVAC_DIR) // "/"
        call set_fgout(.false., 4, TRIM(ftemp) // "fgout_grids.data")

        AVAC_fgrid = FGOUT_fgrids(1)
        PRINT *, AVAC_fgrid%mx
        PRINT *, AVAC_fgrid%my
        PRINT *, AVAC_fgrid%x_low
        PRINT *, AVAC_fgrid%x_hi
        PRINT *, AVAC_fgrid%y_low
        PRINT *, AVAC_fgrid%y_hi

        ALLOCATE(times(SIZE(AVAC_fgrid%output_times)))
        times = AVAC_fgrid%output_times

        DEALLOCATE(FGOUT_fgrids)
        module_setup = .false.
        PRINT *, "Allocating q_avac to ", &
            SIZE(times),&
            AVAC_fgrid%nqout,&
            AVAC_fgrid%mx,&
            AVAC_fgrid%my
        ALLOCATE(q_avac(&
            SIZE(times),&
            AVAC_fgrid%nqout,&
            AVAC_fgrid%mx,&
            AVAC_fgrid%my&
        ))

        ftemp = TRIM(ftemp) // "fgout0001."
        DO i = 1, size(times)-1
            WRITE(file,"(A,A,I0.4)") TRIM(ftemp), "b", i
            PRINT *, "READING FGOUT: ", TRIM(file)
            OPEN(2,FILE=file,ACCESS="stream", STATUS="old", ACTION="read")
                READ(2) q_avac(i,:,:,:) ! ti,qi,xi,yj
            CLOSE(2)
            IF (ANY(ISNAN(q_avac(i,:,:,:))))  THEN
                PRINT *, q_avac(i,:,:,:)
            END IF
        END DO

    END SUBROUTINE init_src_fgout_bin


!     REAL(KIND=8) FUNCTION fgoutinterp(fg, q, x, y)

!         REAL(KIND=8), INTENT(IN) :: q(:,:), x, y
!         TYPE(fgout_grid), INTENT(IN) :: fg
!         REAL(KIND=8) :: xw, xe, ys, yn
!         INTEGER :: w, s

!         IF (fg%x_hi<x .or. x<fg%x_low .or. fg%y_hi<y .or. y<fg%y_low) THEN
!             fgoutinterp = 0.d0
!         ELSE

!             w = MIN(fg%mx-1, 1+INT((x-fg%x_low)/(fg%x_hi-fg%x_low)*(fg%mx-1)))
!             s = MIN(fg%my-1, 1+INT((y-fg%y_low)/(fg%y_hi-fg%y_low)*(fg%my-1)))
!             xw = fg%x_low + (w-1)*(fg%x_hi-fg%x_low)/(fg%mx-1)
!             ys = fg%y_low + (s-1)*(fg%y_hi-fg%y_low)/(fg%my-1)
!             xe = xw + (fg%x_hi-fg%x_low)/(fg%mx-1)
!             yn = ys + (fg%y_hi-fg%y_low)/(fg%my-1)

!             fgoutinterp = (&
!                 + (x-xw) * (y-ys) * q(w+1,s+1) &
!                 + (x-xw) * (yn-y) * q(w+1,s) &
!                 + (xe-x) * (y-ys) * q(w,s+1) &
!                 + (xe-x) * (yn-y) * q(w,s) &
!             ) / ((xe-xw) * (yn-ys))
!         END IF

!     END FUNCTION fgoutinterp

end module helpers
