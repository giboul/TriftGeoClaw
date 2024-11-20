module helpers

    use IEEE_ARITHMETIC, only : ieee_value, ieee_positive_inf
    implicit none
    save

    real(kind=8), allocatable :: q_avac(:,:,:,:)
    real(kind=8), allocatable :: times(:)
    real(kind=8) :: damping = 500.d0 / 1e3
    character(len=4) :: avid
    character(len=3) :: inflow_mode
    INTEGER :: mx, my

contains

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

    INTEGER FUNCTION next_closest(i_closest, value, array)

        INTEGER, INTENT(IN) :: i_closest
        REAL(KIND=8), INTENT(IN) :: value
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: array
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ma

        ma = [array(:(i_closest-1)), &
              ieee_value(array(i_closest),  ieee_positive_inf), &
              array((i_closest+1):)]
        next_closest = closest(value, ma)

    END FUNCTION next_closest

    FUNCTION interp2contours(x, y, x1, y1, q1, x2, y2, q2)

        REAL(KIND=8), DIMENSION(3) :: interp2contours
        REAL(KIND=8), INTENT(IN) :: x, y
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: x1, y1, x2, y2
        REAL(KIND=8), INTENT(IN), DIMENSION(:,:) :: q1, q2
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dist2
        INTEGER :: c1, n1
        INTEGER :: c2, n2
        REAL(KIND=8) :: h, hu, hv

        dist2 = (x-x1)**2 + (y-y1)**2
        c1 = closest(0.d0, dist2)
        n1 = next_closest(c1, 0.d0, dist2)
        dist2 = (x-x2)**2 + (y-y2)**2
        c2 = closest(0.d0, dist2)
        n2 = next_closest(c2, 0.d0, dist2)

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
        ! An alternative bilinear interpolation method between spherical grids
        REAL(KIND=8), INTENT(IN) :: x, y, &
                                    x1, x2, x3, x4, &
                                    y1, y2, y3, y4, &
                                    z1, z2, z3, z4
        REAL(KIND=8), DIMENSION(4,4) :: D0, D1, D2, D3, D4
        REAL(KIND=8) :: a, b, c, d, e

        IF ((x<x1 .and. x<x2 .and. x<x3 .and. x<x4) .or. &
             (x1<x .and. x2<x .and. x3<x .and. x4<x) .or. &
             (y<y1 .and. y<y2 .and. y<y3 .and. y<y4) .or. &
             (y1<y .and. y2<y .and. y3<y .and. y4<y)) THEN
             quadrangular_interp = 0.d0
             RETURN
         END IF 

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
         call det(D0, 4, e)

         IF (e>1.0E-8) THEN
             quadrangular_interp = (a + b*x + c*y + d*x*y)/e
             ! PRINT *, quadrangular_interp
         ELSE
             quadrangular_interp = (z1+z2+z3+z4)/4
         END IF
         quadrangular_interp = MAX(0.d0, quadrangular_interp)

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
        do c = 1, 2
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
   
        do c = 1, 2
            do i = 1, num_times
                write(fname,"(A,I0.1,A,I0.4, A4)") "cut",c,"_",i-1,".npy"
                print "(A,A)", "Reading ", trim(fname)
                open(unit, file=fname, status="old")
                do n = 1, num_cells
                    read(unit,*, iostat=io) x, y, h, hu, hv
                    if (io == 0) then
                        data(c, i, n, 1) = x
                        data(c, i, n, 2) = y
                        data(c, i, n, 3) = h
                        data(c, i, n, 4) = hu
                        data(c, i, n, 5) = hv
                    else
                        data(c, i, n, :) = ieee_value(x, ieee_positive_inf)
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
            ftemp = trim(sides(mthbc))//"_"
            do i = 1, num_times
                n = 0
                write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".npy"
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
            ftemp = trim(sides(mthbc))//"_"
            do i = 1, num_times
                write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".npy"
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

end module helpers
