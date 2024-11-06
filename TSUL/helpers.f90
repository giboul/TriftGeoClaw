MODULE helpers
    implicit none
    save

    REAL(kind=8), ALLOCATABLE :: q_avac(:,:,:,:)
    REAL(kind=8), ALLOCATABLE :: x_src(:), y_src(:), q_src(:,:,:,:)
    REAL(kind=8), ALLOCATABLE :: times(:)
    REAL(kind=8) :: damping = 500.d0 / 1e3
    CHARACTER(len=4) :: avid
    CHARACTER(len=3) :: inflow_mode
    INTEGER :: mx, my

CONTAINS


    SUBROUTINE read_data(avid, mode, mx, my)

        CHARACTER(LEN=3), INTENT(OUT) :: mode
        CHARACTER(LEN=4), INTENT(OUT) :: avid
        INTEGER, INTENT(OUT) :: mx, my
        INTEGER :: i, unit=2

        OPEN(unit, FILE="../setprob.data", STATUS='old')
            DO i=1,6
                READ(unit,*)
            END DO
            READ(unit,*) avid
            READ(unit,*) mode
            READ(unit,*) mx 
            READ(unit,*) my 
        CLOSE(unit)

        IF (mode == "None") then
            mode = "bc"
        END IF
        IF (avid == "None") then
            avid = ""
        END IF
        ! ALLOCATE(data(1, my, mx, 5)) ! data(iy, ix, (x, y, h, hu, hv))

    END SUBROUTINE read_data

    INTEGER FUNCTION closest(value, array)

        REAL(kind=8), INTENT(IN) :: value
        REAL(kind=8), dimension(:), INTENT(IN) :: array
        INTEGER, dimension(1) :: tarray

        tarray = MINLOC(ABS(array-value))
        closest = tarray(1)

    END FUNCTION closest


    INTEGER FUNCTION closest_inf(value, array)

        REAL(kind=8), INTENT(IN) :: value
        REAL(kind=8), dimension(:), INTENT(IN) :: array
        INTEGER, dimension(1) :: tarray

        tarray = MINLOC(ABS(array-value), mask=array<value)
        IF (tarray(1) < 1) THEN
            tarray = MINLOC(array)
        ELSE IF (tarray(1) > size(array)) THEN
            tarray = MAXLOC(array)
        END IF
        closest_inf = tarray(1)

    END FUNCTION closest_inf


    INTEGER FUNCTION closest_sup(value, array)

        REAL(kind=8), INTENT(IN) :: value
        REAL(kind=8), dimension(:), INTENT(IN) :: array
        INTEGER, dimension(1) :: tarray

        tarray = MINLOC(ABS(array-value), mask=array>value)
        IF (tarray(1) < 1) THEN
            tarray = MINLOC(array)
        ELSE IF (tarray(1) > size(array)) THEN
            tarray = MAXLOC(array)
        END IF
        closest_sup = tarray(1)

    END FUNCTION closest_sup


    FUNCTION interp1d4d(t, ts, q)

        REAL(kind=8), INTENT(IN) :: t, ts(:), q(:,:,:,:)
        REAL(kind=8), & 
        dimension(SIZE(q,1),SIZE(q,3),SIZE(q,4)) :: interp1d4d
        INTEGER :: i, s

        i = closest_inf(t, ts)
        s = closest_sup(t, ts)
        interp1d4d = q(:,i,:,:) + &
                    (q(:,s,:,:) - q(:,i,:,:)) * &
                    (MAX(ts(i),MIN(ts(s),t)) - ts(i)) / &
                    (ts(s) - ts(i))

    END FUNCTION interp1d4d


    FUNCTION interp1d2d(xnew, x, q)

        REAL(kind=8), INTENT(IN) :: xnew, x(:), q(:,:)
        REAL(kind=8), dimension(SIZE(q,2)) :: interp1d2d
        INTEGER :: i, s

        i = closest_inf(xnew, x)
        s = closest_sup(xnew, x)
        IF (x(i)==x(s)) THEN
            interp1d2d = q(i,:)
        ELSE
            interp1d2d = q(i,:) + &
                        (q(s,:) - q(i,:)) * &
                        (MAX(x(i),MIN(x(s),xnew)) - x(i)) / &
                        (x(s) - x(i))
        END IF

    END FUNCTION interp1d2d


    SUBROUTINE set_q_src_time(q_src, t, times, avid)

        REAL(KIND=8), INTENT(IN) :: t, times(:)
        REAL(KIND=8), INTENT(INOUT) :: q_src(:,:,:,:)
        CHARACTER(LEN=4), INTENT(IN) :: avid
        CHARACTER(LEN=255), PARAMETER :: ftemp = "../../AVAC/_output"
        INTEGER, PARAMETER :: unit = 2
        REAL(KIND=8), ALLOCATABLE :: q1(:,:,:), q2(:,:,:)
        CHARACTER(len=255) :: file
        REAL(KIND=8) :: h, hu, hv, B, xlow, ylow, dx, dy
        INTEGER :: it, mx, my, i, j

        it = closest_inf(t, times)
        write(file,"(A,I0.4)") TRIM(ftemp)//TRIM(avid)//"/fgout0001.q", it
        call read_fgout_ascii(it, avid, q1)
        call read_fgout_ascii(it+1, avid, q2)

        ! write(file,"(A,I0.4)") TRIM(ftemp)//TRIM(avid)//"/fgout0001.q", it+1
        ! OPEN(unit,FILE=file,STATUS="old",FORM="UNFORMATTED",ACCESS="stream")
        !     READ(unit) q2
        ! CLOSE(unit)

        q_src(1,:,:,:) = q1(:,:,:) + (q2(:,:,:)-q1(:,:,:)) & 
            *(t-times(it))/(times(it+1)-times(it))
        print *, MINVAL(q_src), MAXVAL(q_src)
        print *, MINVAL(q1), MAXVAL(q1)
        print *, MINVAL(q2), MAXVAL(q2)

    END SUBROUTINE set_q_src_time

    SUBROUTINE read_fgout(it, avid, q)

        REAL(KIND=8), ALLOCATABLE, INTENT(INOUT) :: q(:,:,:)
        INTEGER, INTENT(IN) :: it
        CHARACTER(LEN=4), INTENT(IN) :: avid
        CHARACTER(LEN=255), PARAMETER :: ftemp = "../../AVAC/_output"
        REAL(KIND=8) :: h, hu, hv, B, xlow, ylow, dx, dy
        CHARACTER(len=255) :: file
        INTEGER :: mx, my, i, j, unit

        WRITE(file,"(A,I0.4)")TRIM(ftemp)//TRIM(avid)//"/fgout0001.q",it
        OPEN(unit,FILE=TRIM(file),STATUS="old")
            READ(unit,*)
            READ(unit,*)
            READ(unit,*) mx
            READ(unit,*) my
            READ(unit,*) xlow
            READ(unit,*) ylow
            READ(unit,*) dx
            READ(unit,*) dy
        CLOSE(unit)
        print *, mx, my, xlow, ylow, dx, dy
        ALLOCATE(q(mx,my,6))
        WRITE(file,"(A,I0.4)") TRIM(ftemp)//TRIM(avid)//"/fgout0001.b", it
        OPEN(unit,FILE=TRIM(file),STATUS="old",ACCESS="stream")
            READ(unit) q(1:20,1:20,3:6)
        CLOSE(unit)
        print *, q(1:20,1:20,3:6)

    END SUBROUTINE read_fgout

    SUBROUTINE read_fgout_ascii(it, avid, q)

        REAL(KIND=8), ALLOCATABLE, INTENT(INOUT) :: q(:,:,:)
        INTEGER, INTENT(IN) :: it
        CHARACTER(LEN=4), INTENT(IN) :: avid
        CHARACTER(LEN=255), PARAMETER :: ftemp = "../../AVAC/_output"
        REAL(KIND=8) :: h, hu, hv, B, xlow, ylow, dx, dy
        CHARACTER(len=255) :: file
        INTEGER :: mx, my, i, j, unit

        write(file,"(A,I0.4)") TRIM(ftemp)//TRIM(avid)//"/fgout0001.q", it
        print "(A,A)", "Reading file ", TRIM(file)
        OPEN(unit,FILE=TRIM(file),STATUS="old")
            DO i = 1, 2
                READ(unit,*)
            END DO
            READ(unit,*) mx
            READ(unit,*) my
            READ(unit,*) xlow
            READ(unit,*) ylow
            READ(unit,*) dx
            READ(unit,*) dy
            ALLOCATE(q(my, mx,5))
            DO j = 1, my
                DO i = 1, mx
                    READ(unit,*) h, hu, hv, B
                    ! print *, h, hu, hv! , B
                    q(j,i,3) = h
                    q(j,i,4) = hu
                    q(j,i,5) = hv
                    ! q(j,i,6) = B
                END DO
                READ(unit,*)
            END DO
        CLOSE(unit)
    END SUBROUTINE read_fgout_ascii

    SUBROUTINE init_src_old(data, avid, num_times)

        REAL(kind=8), ALLOCATABLE, INTENT(INOUT) :: data(:,:,:,:)
        CHARACTER(len=4), INTENT(IN) :: avid
        INTEGER, INTENT(IN) :: num_times
        CHARACTER(len=255) :: ftemp
        CHARACTER(len=255) :: fname
        REAL(kind=8) :: x, y, h, hu, hv
        INTEGER :: unit, io, i, n
        INTEGER :: num_cells

        ftemp = "../_inflows"//TRIM(avid)//"/cut"
        num_cells = 0
        DO i = 1, num_times
            n = 0
            WRITE(fname,"(A,I0.4, A4)") TRIM(ftemp),i-1,".txt"
            PRINT "(A,A)", "Reading ", TRIM(fname)
            OPEN(unit, file=fname, status="old")
                DO
                    READ(unit,*,iostat=io)
                    IF (io /= 0) then
                        exit
                    END IF
                    n = n + 1
                END DO
            CLOSE(unit)
            num_cells = MAX(num_cells, n)
        END DO
    
        ALLOCATE(data(1, num_times, num_cells, 5))
        PRINT "(A,I10)", "SIZE(data)    = ", SIZE(data)
        PRINT "(A,I10)", "SIZE(data, 1) = ", SIZE(data, 1)
        PRINT "(A,I10)", "SIZE(data, 2) = ", SIZE(data, 2)
        PRINT "(A,I10)", "SIZE(data, 3) = ", SIZE(data, 3)
        PRINT "(A,I10)", "SIZE(data, 3) = ", SIZE(data, 4)
   
        DO i = 1, num_times
            WRITE(fname,"(A,I0.4, A4)") TRIM(ftemp),i-1,".txt"
            OPEN(unit, file=fname, status="old")
            DO n = 1, num_cells
                READ(unit,*, iostat=io) x, y, h, hu, hv
                data(1, i, n, 1) = x
                data(1, i, n, 2) = y
                data(1, i, n, 3) = h
                data(1, i, n, 4) = hu
                data(1, i, n, 5) = hv
            END DO
            CLOSE(unit)
        END DO
    END SUBROUTINE init_src_old

    SUBROUTINE init_bc(data, avid, num_times)
        REAL(kind=8), ALLOCATABLE, INTENT(INOUT) :: data(:,:,:,:)
        INTEGER, INTENT(IN) :: num_times
        CHARACTER(len=6), dimension(4) :: sides
        CHARACTER(len=255) :: fdir, ftemp
        CHARACTER(len=255) :: fname
        CHARACTER(len=4) :: avid
        REAL(kind=8) :: x, y, h, hu, hv
        INTEGER :: unit, io, i, n, mthbc
        INTEGER :: num_cells
    
        unit = 2
        sides = [CHARACTER(len=6) :: "left", "right", "bottom", "top"]
 
        PRINT *, "Avalanche id #", TRIM(avid)

        fdir = "../_inflows"//TRIM(avid)//"/"
        ftemp = TRIM(fdir)//TRIM(sides(1))//"_"
        num_cells = 0
        DO mthbc = 1, 4
            ftemp = TRIM(fdir)//TRIM(sides(mthbc))//"_"
            DO i = 1, num_times
                n = 0
                WRITE(fname,"(A,I0.4, A4)") TRIM(ftemp),i-1,".txt"
                PRINT "(A,A)", "Reading ", TRIM(fname)
                OPEN(unit, file=fname, status="old")
                    DO
                        READ(unit,*,iostat=io)
                        IF (io /= 0) then
                            exit
                        END IF
                        n = n + 1
                    END DO
                CLOSE(unit)
                num_cells = MAX(num_cells, n)
            END DO
        END DO
    
        ALLOCATE(data(4, num_times, num_cells, 5))
        PRINT "(A,I10)", "SIZE(data)    = ", SIZE(data)
        PRINT "(A,I10)", "SIZE(data, 1) = ", SIZE(data, 1)
        PRINT "(A,I10)", "SIZE(data, 2) = ", SIZE(data, 2)
        PRINT "(A,I10)", "SIZE(data, 3) = ", SIZE(data, 3)
        PRINT "(A,I10)", "SIZE(data, 3) = ", SIZE(data, 4)
   
        DO mthbc = 1, 4
            ftemp = TRIM(fdir)//TRIM(sides(mthbc))//"_"
            DO i = 1, num_times
                WRITE(fname,"(A,I0.4, A4)") TRIM(ftemp),i-1,".txt"
                OPEN(unit, file=fname, status="old")
                DO n = 1, num_cells
                    READ(unit,*, iostat=io) x, y, h, hu, hv
                    data(mthbc, i, n, 1) = x
                    data(mthbc, i, n, 2) = y
                    data(mthbc, i, n, 3) = h
                    data(mthbc, i, n, 4) = hu
                    data(mthbc, i, n, 5) = hv
                END DO
                CLOSE(unit)
            END DO
        END DO
    END SUBROUTINE init_bc

    SUBROUTINE read_times(times, avid)
        CHARACTER(len=4), INTENT(IN) :: avid
        CHARACTER(len=255) :: fname
        INTEGER :: io, n, i
        INTEGER :: unit
        REAL(kind=8), ALLOCATABLE :: times(:)

        unit = 2
        fname = "../_inflows"//TRIM(avid)//"/timing.txt"
        PRINT "(A,A)", "Reading ", TRIM(fname)
        OPEN(unit, file=fname)
            n = 0
            DO 
                READ(unit,*,iostat=io)
                IF (io /= 0) then
                    exit
                END IF
                n = n + 1
            END DO
            ALLOCATE(times(n))
        rewind(unit)
            DO i = 1, n
                READ(unit,*) times(i)
            END DO
        CLOSE(unit)
    END SUBROUTINE read_times

    subroutine gridinterp(x, y, Z, xnew, ynew, Znew)
        real(kind=8), intent(in) :: x(:), y(:), Z(:,:)
        real(kind=8), intent(out) :: Znew(:,:)
        real(kind=8), intent(inout) :: xnew(:), ynew(:)
        integer :: i, j, w, e, s, n

        do j = 1, size(ynew)
            do i = 1, size(xnew)
                w = closest_inf(xnew(i), x)
                s = closest_inf(ynew(j), y)
                e = w + 1
                n = s + 1
                Znew(j, i) = (&
                    + (xnew(i)-x(w)) * (ynew(j)-y(s)) * Z(n,e) &
                    + (xnew(i)-x(w)) * (y(n)-ynew(j)) * Z(s,e) &
                    + (x(e)-xnew(i)) * (ynew(j)-y(s)) * Z(n,w) &
                    + (x(e)-xnew(i)) * (y(n)-ynew(j)) * Z(s,w) &
                ) / ((x(e)-x(w)) * (y(n)-y(s)))
            end do
        end do
    end subroutine gridinterp

END module helpers
