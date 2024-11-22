module helpers

    use IEEE_ARITHMETIC, only : ieee_value, ieee_positive_inf
    use fgout_module, only : fgout_grid, set_fgout, FGOUT_fgrids
    implicit none
    save

    real(kind=8), allocatable :: q_avac(:,:,:,:)
    real(kind=8), allocatable :: times(:)
    real(kind=8) :: damping = 500.d0 / 1e3
    character(len=4) :: avid
    character(len=3) :: inflow_mode

contains

    SUBROUTINE read_data(avid, mode)

        CHARACTER(LEN=3), INTENT(OUT) :: mode
        CHARACTER(LEN=4), INTENT(OUT) :: avid
        INTEGER :: i, unit=2

        OPEN(unit, FILE="../setprob.data", STATUS='old')
            DO i=1,6
                READ(unit,*)
            END DO
            READ(unit,*) avid
            READ(unit,*) mode
        CLOSE(unit)

        IF (TRIM(mode) == "None") then
            mode = "bc"
        END IF
        IF (TRIM(avid) == "None") then
            avid = ""
        END IF

    END SUBROUTINE read_data

    subroutine read_times(times)
        character(len=255) :: fname
        integer :: io, n, i
        integer :: unit
        real(kind=8), allocatable :: times(:)

        unit = 2
        fname = "times.txt"
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

    integer function closest(value, array)

        real(kind=8), intent(in) :: value
        real(kind=8), dimension(:), intent(in) :: array

        closest = MINLOC(ABS(array-value), DIM=1)

    end function closest

    INTEGER FUNCTION closest_inf(value, array)

        REAL(KIND=8), INTENT(IN) :: value
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: array

        closest_inf = MINLOC(ABS(array-value), MASK=array<value, DIM=1)
        IF (closest_inf < 1) THEN
            closest_inf = MINLOC(array, DIM=1)
        END IF
        closest_inf = MIN(SIZE(array)-1, closest_inf)

    END FUNCTION closest_inf

    INTEGER FUNCTION next_closest(i1, x, y, xarray, yarray)

        INTEGER, INTENT(IN) :: i1
        REAL(KIND=8), INTENT(IN) :: x, y
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xarray, yarray
        REAL(KIND=8) :: a1, a2, b1, b2, c1, c2
        INTEGER :: i2

        i2 = i1 + 1
        a1 = xarray(i2) - xarray(i1)
        a2 = yarray(i2) - yarray(i1)
        b1 = xarray(i1) - x
        b2 = yarray(i1) - y
        c1 = xarray(i2) - x
        c2 = yarray(i2) - y
        IF (a1*(b1+c1)+a2*(b2+c2)<0) THEN
            next_closest = i2
        ELSE
            next_closest = i1 -1
        END IF

    END FUNCTION next_closest

    FUNCTION interp2contours(x, y, x1, y1, q1, x2, y2, q2)

        REAL(KIND=8), DIMENSION(3) :: interp2contours
        REAL(KIND=8), INTENT(IN) :: x, y
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: x1, y1, x2, y2
        REAL(KIND=8), INTENT(IN), DIMENSION(:,:) :: q1, q2
        INTEGER :: c1, n1
        INTEGER :: c2, n2
        REAL(KIND=8) :: h, hu, hv

        c1 = closest(0.d0, (x-x1)**2+(y-y1**2))
        c2 = closest(0.d0, (x-x2)**2+(y-y2**2))
        n1 = next_closest(c1, x, y, x1, y1)
        n2 = next_closest(c2, x, y, x2, y2)

        h = quadrangular_interp(x, y, &
            x1(c1), y1(c1), q1(c1,1), &
            x1(n1), y1(n1), q1(n1,1), &
            x2(c2), y2(c2), q2(c2,1), &
            x2(n2), y2(n2), q2(n2,1) &
        )
        hu = quadrangular_interp(x, y, &
            x1(c1), y1(c1), q1(c1,2), &
            x1(n1), y1(n1), q1(n1,2), &
            x2(c2), y2(c2), q2(c2,2), &
            x2(n2), y2(n2), q2(n2,2) &
        )
        hv = quadrangular_interp(x, y, &
            x1(c1), y1(c1), q1(c1,3), &
            x1(n1), y1(n1), q1(n1,3), &
            x2(c2), y2(c2), q2(c2,3), &
            x2(n2), y2(n2), q2(n2,3) &
        )

        interp2contours = [h, hu, hv]

    END FUNCTION interp2contours

    REAL(KIND=8) FUNCTION quadrangular_interp(x, y, &
                                              x1, y1, z1, &
                                              x2, y2, z2, &
                                              x3, y3, z3, &
                                              x4, y4, z4)
        ! doi:10.3390/atmos10030123 
        ! An Alternative Bilinear Interpolation Method Between Spherical Grids
        REAL(KIND=8), INTENT(IN) :: x, y, &
                                    x1, x2, x3, x4, &
                                    y1, y2, y3, y4, &
                                    z1, z2, z3, z4
        REAL(KIND=8), DIMENSION(4,4) :: D0, D1, D2, D3, D4
        REAL(KIND=8) :: a, b, c, d, f

        D0(1,:) = [1.d0, x1, y1, x1*y1]
        D0(2,:) = [1.d0, x2, y2, x2*y2]
        D0(3,:) = [1.d0, x3, y3, x3*y3]
        D0(4,:) = [1.d0, x4, y4, x4*y4]

        D1(1,:) = [z1, x1, y1, x1*y1]
        D1(2,:) = [z2, x2, y2, x2*y2]
        D1(3,:) = [z3, x3, y3, x3*y3]
        D1(4,:) = [z4, x4, y4, x4*y4]

        D2(1,:) = [1.d0, z1, y1, x1*y1]
        D2(2,:) = [1.d0, z2, y2, x2*y2]
        D2(3,:) = [1.d0, z3, y3, x3*y3]
        D2(4,:) = [1.d0, z4, y4, x4*y4]

        D3(1,:) = [1.d0, x1, z1, x1*y1]
        D3(2,:) = [1.d0, x2, z2, x2*y2]
        D3(3,:) = [1.d0, x3, z3, x3*y3]
        D3(4,:) = [1.d0, x4, z4, x4*y4]

        D4(1,:) = [1.d0, x1, y1, z1]
        D4(2,:) = [1.d0, x2, y2, z2]
        D4(3,:) = [1.d0, x3, y3, z3]
        D4(4,:) = [1.d0, x4, y4, z4]

        call det(D1, 4, a)
        call det(D2, 4, b)
        call det(D3, 4, c)
        call det(D4, 4, d)
        call det(D0, 4, f)
        quadrangular_interp = (a + b*x + c*y + d*x*y)/f
        IF (ISNAN(quadrangular_interp)) THEN
            STOP
        END IF

    END FUNCTION quadrangular_interp

    SUBROUTINE det(A, m, D)                                              

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: m
        REAL(KIND=8), DIMENSION(m,m), INTENT(IN) :: A
        INTEGER :: i
        INTEGER, DIMENSION(m) :: ipiv
        REAL(KIND=8), INTENT(OUT) :: D

        D = 0.d0
        CALL dgetrf(m, m, A, m, ipiv, i)

        D = 1.0d0
        DO i = 1, m
            D = D * A(i,i)
            IF (ipiv(i) /= i) THEN
                D = -D
            END IF
        END DO

    END SUBROUTINE det

    REAL(KIND=8) FUNCTION det4x4(D)
        ! Because numpy and lapack use LU decomposition, I wanted a direct method
        REAL(KIND=8), DIMENSION(4,4), INTENT(IN) :: D
        print *, D(1,:)
        print *, D(2,:)
        print *, D(3,:)
        print *, D(4,:)
        det4x4 = &
         + D(1,1)*D(2,2)*D(3,3)*D(4,4) + D(1,1)*D(2,3)*D(3,4)*D(4,2) + D(1,1)*D(2,4)*D(3,2)*D(4,3) &
         - D(1,1)*D(2,4)*D(3,3)*D(4,2) - D(1,1)*D(2,3)*D(3,2)*D(4,4) - D(1,1)*D(2,2)*D(3,4)*D(4,3) &
         - D(1,2)*D(2,1)*D(3,3)*D(4,4) - D(1,3)*D(2,1)*D(3,4)*D(4,2) - D(1,4)*D(2,1)*D(3,2)*D(4,3) &
         + D(1,4)*D(2,1)*D(3,3)*D(4,2) + D(1,3)*D(2,1)*D(3,2)*D(4,4) + D(1,2)*D(2,1)*D(3,4)*D(4,3) &
         + D(1,2)*D(2,3)*D(3,1)*D(4,4) + D(1,3)*D(2,4)*D(3,1)*D(4,2) + D(1,4)*D(2,2)*D(3,1)*D(4,3) &
         - D(1,4)*D(2,3)*D(3,1)*D(4,2) - D(1,3)*D(2,2)*D(3,1)*D(4,4) - D(1,2)*D(2,4)*D(3,1)*D(4,3) &
         - D(1,2)*D(2,3)*D(3,4)*D(4,1) - D(1,3)*D(2,4)*D(3,2)*D(4,1) - D(1,4)*D(2,2)*D(3,3)*D(4,1) &
         + D(1,4)*D(2,3)*D(3,2)*D(4,1) + D(1,3)*D(2,2)*D(3,4)*D(4,1) + D(1,2)*D(2,4)*D(3,3)*D(4,1) 
        print *, "==>", det4x4
    END FUNCTION det4x4

    real(kind=8) function interp(newx, x, y)
        real(kind=8), intent(in) :: x(:), y(:), newx
        real(kind=8) :: xp, yp
        integer :: ix

        ix = closest_inf(newx, x)
        if (ix == 0) then
            interp = y(1)
        else if (ix == size(x)) then
            interp = y(size(y))
        else
            xp = x(ix)
            yp = y(ix)
            interp = yp + (newx-xp)/(x(ix+1)-xp) * (y(ix+1)-yp)
        end if
    end function interp

    subroutine init_src(data, num_times)
        real(kind=8), allocatable, intent(inout) :: data(:,:,:,:)
        integer, intent(in) :: num_times
        character(len=255) :: fname
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n, c
        integer :: num_cells
    
        num_cells = 0
        do c = 2, 3
            do i = 1, num_times
                n = 0
                write(fname,"(A,I0.1,A,I0.4, A4)") "cut",c,"_",i-1,".npy"
                print "(A,A)", "Reading ", trim(fname)
                open(unit, file=fname, status="old")
                    do
                        read(unit,*,iostat=io)
                        if (io /= 0) then
                            exit
                        end if
                        n = n + 1
                    end do
                close(unit)
                num_cells = max(num_cells, n)
            end do
        end do
 
        allocate(data(2, num_times, num_cells, 5))
        print "(A,I10)", "size(data)    = ", size(data)
        print "(A,I10)", "size(data, 1) = ", size(data, 1)
        print "(A,I10)", "size(data, 2) = ", size(data, 2)
        print "(A,I10)", "size(data, 3) = ", size(data, 3)
        print "(A,I10)", "size(data, 3) = ", size(data, 4)
   
        do c = 2, 3
            do i = 1, num_times
                write(fname,"(A,I0.1,A,I0.4, A4)") "cut",c,"_",i-1,".npy"
                open(unit, file=fname, status="old")
                do n = 1, num_cells
                    read(unit,*, iostat=io) x, y, h, hu, hv
                    if (io == 0) then
                        data(c-1, i, n, 1) = x
                        data(c-1, i, n, 2) = y
                        data(c-1, i, n, 3) = h
                        data(c-1, i, n, 4) = hu
                        data(c-1, i, n, 5) = hv
                    else
                        data(c-1, i, n, :) = ieee_value(x, ieee_positive_inf)
                    end if
                end do
                close(unit)
            end do
        end do
    end subroutine init_src

    subroutine init_bc(data, num_times)
        real(kind=8), allocatable, intent(inout) :: data(:,:,:,:)
        integer, intent(in) :: num_times
        character(len=6), dimension(4) :: sides
        character(len=255) :: ftemp
        character(len=255) :: fname
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n, mthbc
        integer :: num_cells
    
        unit = 2
        sides = [character(len=6) :: "left", "right", "bottom", "top"]
 
        num_cells = 0
        do mthbc = 1, 4
            do i = 1, num_times
                n = 0
                write(fname,"(A,I0.4, A4)") trim(sides(mthbc))//"_",i-1,".npy"
                print "(A,A)", "Reading ", trim(fname)
                open(unit, file=fname, status="old")
                    do
                        read(unit,*,iostat=io)
                        if (io /= 0) then
                            exit
                        end if
                        n = n + 1
                    end do
                close(unit)
                num_cells = max(num_cells, n)
            end do
        end do
    
        allocate(data(4, num_times, num_cells, 5))
        print "(A,I10)", "size(data)    = ", size(data)
        print "(A,I10)", "size(data, 1) = ", size(data, 1)
        print "(A,I10)", "size(data, 2) = ", size(data, 2)
        print "(A,I10)", "size(data, 3) = ", size(data, 3)
        print "(A,I10)", "size(data, 3) = ", size(data, 4)
   
        do mthbc = 1, 4
            do i = 1, num_times
                write(fname,"(A,I0.4, A4)") trim(sides(mthbc))//"_",i-1,".npy"
                open(unit, file=fname, status="old")
                do n = 1, num_cells
                    read(unit,*, iostat=io) x, y, h, hu, hv
                    data(mthbc, i, n, 1) = x
                    data(mthbc, i, n, 2) = y
                    data(mthbc, i, n, 3) = h
                    data(mthbc, i, n, 4) = hu
                    data(mthbc, i, n, 5) = hv
                end do
                close(unit)
            end do
        end do
    end subroutine init_bc


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

    ! SUBROUTINE set_q_src_time(q_src, t, times, avid)

    !     REAL(KIND=8), INTENT(IN) :: t, times(:)
    !     REAL(KIND=8), INTENT(INOUT) :: q_src(:,:,:,:)
    !     CHARACTER(LEN=4), INTENT(IN) :: avid
    !     CHARACTER(LEN=255), PARAMETER :: ftemp = "../../AVAC/_output"
    !     INTEGER, PARAMETER :: unit = 2
    !     REAL(KIND=8), ALLOCATABLE :: q1(:,:,:), q2(:,:,:)
    !     CHARACTER(len=255) :: file
    !     REAL(KIND=8) :: h, hu, hv, B, xlow, ylow, dx, dy
    !     INTEGER :: it, mx, my, i, j

    !     it = closest_inf(t, times)
    !     write(file,"(A,I0.4)") TRIM(ftemp)//TRIM(avid)//"/fgout0001.q", it
    !     call read_fgout(it, avid, q1)
    !     call read_fgout(it+1, avid, q2)

    !     ! write(file,"(A,I0.4)") TRIM(ftemp)//TRIM(avid)//"/fgout0001.q", it+1
    !     ! OPEN(unit,FILE=file,STATUS="old",FORM="UNFORMATTED",ACCESS="stream")
    !     !     READ(unit) q2
    !     ! CLOSE(unit)

    !     q_src(1,:,:,:) = q1(:,:,:) + (q2(:,:,:)-q1(:,:,:)) & 
    !         *(t-times(it))/(times(it+1)-times(it))
    !     print *, MINVAL(q_src), MAXVAL(q_src)
    !     print *, MINVAL(q1), MAXVAL(q1)
    !     print *, MINVAL(q2), MAXVAL(q2)

    ! END SUBROUTINE set_q_src_time

    SUBROUTINE init_src_fgout_bin(avid, num_times, q_avac, FGOUT_fgrids)

        CHARACTER(LEN=4), INTENT(IN) :: avid
        INTEGER, INTENT(IN) :: num_times
        REAL(KIND=8), ALLOCATABLE, INTENT(INOUT) :: q_avac(:,:,:,:)
        TYPE(fgout_grid), ALLOCATABLE :: FGOUT_fgrids(:)
        CHARACTER(len=255) :: file, ftemp
        INTEGER :: i, j, k

        ftemp = "../../AVAC/_output" // trim(avid) // "/"
        call set_fgout(.false., 4, TRIM(ftemp) // "fgout_grids.data")
        ftemp = TRIM(ftemp) // "fgout0001.b"

        PRINT *, FGOUT_fgrids(1)%mx
        PRINT *, FGOUT_fgrids(1)%my
        PRINT *, FGOUT_fgrids(1)%x_low
        PRINT *, FGOUT_fgrids(1)%x_hi
        PRINT *, FGOUT_fgrids(1)%y_low
        PRINT *, FGOUT_fgrids(1)%y_hi

        ALLOCATE(q_avac(&
            size(times),&
            FGOUT_fgrids(1)%nqout,&
            FGOUT_fgrids(1)%mx,&
            FGOUT_fgrids(1)%my&
        ))

        DO i = 1, num_times-1
            WRITE(file,"(A,I0.4)") TRIM(ftemp), i
            PRINT *, "READING FGOUT: ", TRIM(file)
            OPEN(2,FILE=file,ACCESS="stream", STATUS="old", ACTION="read")
                READ(2) q_avac(i,:,:,:) ! ti,qi,xi,yj
            CLOSE(2)
        END DO

        ! OPEN(2,file="/home/giboul/Downloads/data.npy")
        !     DO i = 1, SIZE(q_avac,3)
        !         DO j = 1, SIZE(q_avac,4)
        !             WRITE(2,"(E20.10E5)",advance="no") q_avac(1,3,i,j)
        !             WRITE(2,"(A1)",advance="no") " "
        !         END DO
        !         WRITE(2,*) ""
        !     END DO
        ! CLOSE(2)

    END SUBROUTINE init_src_fgout_bin

    SUBROUTINE read_fgout_ascii(it, avid, q)

        REAL(KIND=8), ALLOCATABLE, INTENT(INOUT) :: q(:,:,:)
        INTEGER, INTENT(IN) :: it
        CHARACTER(LEN=4), INTENT(IN) :: avid
        CHARACTER(LEN=255), PARAMETER :: ftemp = "_output"
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

    real(kind=8) function gridinterp(x, y, Z, xnew, ynew)

        real(kind=8), intent(in) :: x(:), y(:), Z(:,:), xnew, ynew
        integer :: i, j, w, e, s, n

        w = closest_inf(xnew, x)
        s = closest_inf(ynew, y)
        e = w + 1
        n = s + 1
        gridinterp = (&
            + (xnew-x(w)) * (ynew-y(s)) * Z(n,e) &
            + (xnew-x(w)) * (y(n)-ynew) * Z(s,e) &
            + (x(e)-xnew) * (ynew-y(s)) * Z(n,w) &
            + (x(e)-xnew) * (y(n)-ynew) * Z(s,w) &
        ) / ((x(e)-x(w)) * (y(n)-y(s)))

    end function gridinterp

    real(kind=8) function fgoutinterp(fg, q, xnew, ynew)

        real(kind=8), intent(in) :: q(:,:), xnew, ynew
        type(fgout_grid), intent(in) :: fg
        REAL(KIND=8) :: x, y, xw, xe, ys, yn
        INTEGER :: w, e, s, n

        x = MAX(fg%x_low,MIN(fg%x_hi,xnew))
        y = MAX(fg%y_low,MIN(fg%y_hi,ynew))
        ! print *, xnew, x
        ! print *, ynew, y

        w = INT((x-fg%x_low)/(fg%x_hi-fg%x_low)*fg%mx)
        s = INT((y-fg%y_low)/(fg%y_hi-fg%y_low)*fg%my)
        ! print *, &
        !     q(w+1,s+1), &
        !     q(w+1,s), &
        !     q(w,s+1), &
        !     q(w,s)
        xw = fg%x_low + w*(fg%x_hi-fg%x_low)/fg%mx
        ys = fg%y_low + s*(fg%y_hi-fg%y_low)/fg%my
        xe = xw + (fg%x_hi-fg%x_low)/fg%mx
        yn = ys + (fg%y_hi-fg%y_low)/fg%my
        fgoutinterp = (&
            + (x-xw) * (y-ys) * q(w+1,s+1) &
            + (x-xw) * (yn-y) * q(w+1,s) &
            + (xe-x) * (y-ys) * q(w,s+1) &
            + (xe-x) * (yn-y) * q(w,s) &
        ) / ((xe-xw) * (yn-ys))
        if (fgoutinterp < 0 .and. .false.) then
            print *, fgoutinterp
            print *, xw, x, xe
            print *, ys, y, yn
            print *, (x<xw) .or. (xe<x) .or. (y<ys) .or. (yn<y)
            print *, ANY(q<0)
            print *, "--------------"
        end if

    end function fgoutinterp

end module helpers
