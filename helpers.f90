module helpers

    use IEEE_ARITHMETIC, only : ieee_value, ieee_positive_inf
    use fgout_module, only : fgout_grid, set_fgout, FGOUT_fgrids, module_setup
    implicit none
    save

    integer :: bc_size, avac_fgno
    real(kind=8), allocatable :: q_avac(:,:,:,:)
    real(kind=8), allocatable :: times(:)
    real(kind=8) :: damping, min_alt_avac
    character(len=255) :: AVAC_DIR, mode, input_format
    type(fgout_grid) :: AVAC_fgrid

contains


    SUBROUTINE read_data()

        INTEGER :: unit=2

        call opendatafile(unit, "setprob.data")
            READ(unit,*) mode
            READ(unit,*) damping
            READ(unit,*) min_alt_avac
            READ(unit,*) AVAC_DIR
            READ(unit,*) bc_size
            READ(unit,*) input_format
            READ(unit,*) avac_fgno
        CLOSE(unit)

        IF (TRIM(mode) == "None") then
            mode = "bc"
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
        w = 1+INT((xc-fg%x_low)/(fg%x_hi-fg%x_low)*(fg%mx-1))
        s = 1+INT((yc-fg%y_low)/(fg%y_hi-fg%y_low)*(fg%my-1))
        w = MIN(fg%mx-1, MAX(w, 1))
        s = MIN(fg%my-1, MAX(s, 1))

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


    SUBROUTINE init_src_fgout()

        CHARACTER(len=255) :: fname, ftemp
        INTEGER :: i, j, k, l
        real(kind=4), allocatable :: q_real4(:,:,:)

        ftemp = trim(AVAC_DIR) // "/"
        call set_fgout(.false., 4, TRIM(ftemp) // "fgout_grids.data")

        AVAC_fgrid = FGOUT_fgrids(avac_fgno)
        PRINT *, "mx    ", AVAC_fgrid%mx
        PRINT *, "my    ", AVAC_fgrid%my
        PRINT *, "x_low ", AVAC_fgrid%x_low
        PRINT *, "x_hi  ", AVAC_fgrid%x_hi
        PRINT *, "y_low ", AVAC_fgrid%y_low
        PRINT *, "y_hi  ", AVAC_fgrid%y_hi
        PRINT *, "times ", SIZE(AVAC_fgrid%output_times)

        ALLOCATE(times(SIZE(AVAC_fgrid%output_times)))
        times = AVAC_fgrid%output_times

        DEALLOCATE(FGOUT_fgrids)
        module_setup = .false.
        ALLOCATE(q_avac(&
            SIZE(times),&
            AVAC_fgrid%nqout,&
            AVAC_fgrid%mx,&
            AVAC_fgrid%my&
        ))
        IF (input_format == "binary32") THEN
            ALLOCATE(q_real4(&
                AVAC_fgrid%nqout,&
                AVAC_fgrid%mx,&
                AVAC_fgrid%my&
            ))
        END IF

        !ftemp = TRIM(ftemp) // "fgout0001."
        write(ftemp, "(A,A,I0.4,A)") TRIM(ftemp), "fgout", avac_fgno,"."
        fname = ""
        PRINT *, "Input format: ", TRIM(input_format)
        DO i = 1, size(times)
            PRINT *, "READING FGOUT GRID", i
            IF (TRIM(input_format) == "ascii") THEN  ! ASCII input
                WRITE(fname,"(A,A,I0.4)") TRIM(ftemp), "q", i
                OPEN(2,FILE=fname, STATUS="old", ACTION="read")
                DO l=1,9
                    READ(2,*)
                END DO
                DO l=1,AVAC_fgrid%my
                    DO k=1,AVAC_fgrid%mx
                        READ(2, "(50e26.16)") (q_avac(i,j,k,l), j=1,AVAC_fgrid%nqout)
                    END DO
                    READ(2,*)
                END DO
                CLOSE(2)
            ELSE IF (TRIM(input_format) == "binary32") THEN  ! REAL(4) input
                WRITE(fname,"(A,A,I0.4)") TRIM(ftemp), "b", i
                OPEN(2,FILE=fname,ACCESS="stream", STATUS="old", ACTION="read")
                    READ(2) q_real4(:,:,:)
                CLOSE(2)
                q_avac(i,:,:,:) = REAL(q_real4, kind=8)
            ELSE IF (TRIM(input_format) == "binary64") THEN  ! REAL(8) input
                WRITE(fname,"(A,A,I0.4)") TRIM(ftemp), "b", i
                OPEN(2,FILE=fname,ACCESS="stream", STATUS="old", ACTION="read")
                    READ(2) q_avac(i,:,:,:) ! ti,qi,xi,yj
                CLOSE(2)
            ELSE
                WRITE(*, "(A,I0.1,A)") "Input fgout output_format",&
                    input_format,&
                    "not supported (2 or 3)."
            END IF
            IF (ANY(ISNAN(q_avac(i,:,:,:))))  THEN
                PRINT *, q_avac(i,:,:,:)
            END IF
        END DO

        IF (input_format == "binary32") THEN
            DEALLOCATE(q_real4)
        END IF

    END SUBROUTINE init_src_fgout


end module helpers
