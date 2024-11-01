MODULE helpers
    implicit none
    save

    REAL(kind=8), ALLOCATABLE :: q_avac(:,:,:,:)
    REAL(kind=8), ALLOCATABLE :: times(:)
    REAL(kind=8) :: damping = 500.d0 / 1e3
    CHARACTER(len=4) :: avid
    CHARACTER(len=3) :: inflow_mode

CONTAINS

    SUBROUTINE read_avid(avid)
        CHARACTER(len=4), INTENT(OUT) :: avid
        OPEN(2, file="../avac.data", status='old')
            READ(2,*) avid
        CLOSE(2)
        IF (avid == "None") then
            avid = ""
        END IF
    END SUBROUTINE read_avid

    SUBROUTINE read_inflow_mode(mode)
        CHARACTER(len=3), INTENT(OUT) :: mode
        OPEN(2, file="../inflow.data", status='old')
            READ(2,*) mode
        CLOSE(2)
        IF (mode == "None") then
            mode = "bc"
        END IF
    END SUBROUTINE read_inflow_mode

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

    END function


    SUBROUTINE init_src(data, avid, num_times)

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
    END SUBROUTINE init_src

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
 
        OPEN(unit, file="../avac.data", status='old')
            READ(unit,*) avid
        CLOSE(unit)
        IF (avid == "None") then
            avid = "    "
        END IF
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

END module helpers
