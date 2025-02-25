module helpers

    use IEEE_ARITHMETIC, only : ieee_value, ieee_positive_inf
    use fgout_module, only : fgout_grid, set_fgout, FGOUT_fgrids, module_setup
    implicit none
    save

    real(kind=8), allocatable :: q_avac(:,:,:,:)
    real(kind=8), allocatable :: times(:)
    real(kind=8) :: damping, overhang
    real(kind=8) :: lake_alt
    character(len=255) :: inflow_mode
    type(fgout_grid) :: AVAC_fgrid

contains


    SUBROUTINE read_data()

        INTEGER :: unit=2

        call opendatafile(unit, "setprob.data")
            READ(unit,*) inflow_mode
            READ(unit,*) damping
            READ(unit,*) lake_alt
            READ(unit,*) overhang
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


    subroutine init_bc()

        character(len=6), dimension(4) :: sides
        character(len=255) :: ftemp
        character(len=255) :: fname
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n, mthbc
        integer :: num_cells
    
        unit = 2
        sides = [character(len=6) :: "left", "right", "bottom", "top"]

        call read_times()
 
        num_cells = 0
        do mthbc = 1, 4
            do i = 1, size(times)
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
 
        allocate(q_avac(size(times), 4, num_cells, 5))
        print "(A,I10)", "size(q_avac)    = ", size(q_avac)
        print "(A,I10)", "size(q_avac, 1) = ", size(q_avac, 1)
        print "(A,I10)", "size(q_avac, 2) = ", size(q_avac, 2)
        print "(A,I10)", "size(q_avac, 3) = ", size(q_avac, 3)
        print "(A,I10)", "size(q_avac, 3) = ", size(q_avac, 4)
   
        do mthbc = 1, 4
            do i = 1, size(times)
                write(fname,"(A,I0.4, A4)") trim(sides(mthbc))//"_",i-1,".npy"
                open(unit, file=fname, status="old")
                do n = 1, num_cells
                    read(unit,*, iostat=io) x, y, h, hu, hv
                    q_avac(i, mthbc, n, 1) = x
                    q_avac(i, mthbc, n, 2) = y
                    q_avac(i, mthbc, n, 3) = h
                    q_avac(i, mthbc, n, 4) = hu
                    q_avac(i, mthbc, n, 5) = hv
                end do
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


    FUNCTION interp1d4d(t, ts, q)

        REAL(kind=8), INTENT(IN) :: t, ts(:), q(:,:,:,:)
        REAL(kind=8), & 
        dimension(SIZE(q,2),SIZE(q,3),SIZE(q,4)) :: interp1d4d
        INTEGER :: i, s

        i = closest_inf(t, ts)
        s = closest_sup(t, ts)
        IF (ts(s)<t) THEN
            interp1d4d = q(i,:,:,:)*0
        ELSE IF (t<ts(i)) THEN
            interp1d4d = q(s,:,:,:)*0
        ELSE
            interp1d4d = q(i,:,:,:) + &
                        (q(s,:,:,:) - q(i,:,:,:)) * &
                        (t - ts(i)) / (ts(s) - ts(i))
        END IF

    END FUNCTION interp1d4d


    FUNCTION interp1d2d(xnew, x, q)

        REAL(kind=8), INTENT(IN) :: xnew, x(:), q(:,:)
        REAL(kind=8), dimension(SIZE(q,2)) :: interp1d2d
        INTEGER :: i, s

        i = closest_inf(xnew, x)
        s = closest_sup(xnew, x)
        IF (.not.(x(i)<=xnew.and.xnew<=x(s))) THEN
            IF (xnew>=MINVAL(x) .and. MAXVAL(x)>=xnew) THEN
                PRINT *, "####", MINVAL(x), x(i), xnew, x(s), MAXVAL(x)
                stop
            END IF
        END IF
        IF (x(i)==x(s)) THEN
            interp1d2d = q(i,:)
        ELSE
            interp1d2d = q(i,:) + &
                        (q(s,:) - q(i,:)) * &
                        (MAX(x(i),MIN(x(s),xnew)) - x(i)) / &
                        (x(s) - x(i))
        END IF

    END FUNCTION interp1d2d


    SUBROUTINE init_src_fgout_bin()

        CHARACTER(len=255) :: file, ftemp
        INTEGER :: i

        ftemp = "../../avac/_output/"
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


    SUBROUTINE read_fgout_ascii(it, q)

        REAL(KIND=8), ALLOCATABLE, INTENT(INOUT) :: q(:,:,:)
        INTEGER, INTENT(IN) :: it
        CHARACTER(LEN=255), PARAMETER :: ftemp = "_output"
        REAL(KIND=8) :: h, hu, hv, B, xlow, ylow, dx, dy
        CHARACTER(len=255) :: file
        INTEGER :: mx, my, i, j, unit

        write(file,"(A,I0.4)") TRIM(ftemp)//"/fgout0001.q", it
        print "(A,A)", "Reading file ", TRIM(file)
        OPEN(unit,FILE=TRIM(file),STATUS="old")
            DO i = 1, 2
                READ(unit,*)
            END DO
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


!     real(kind=8) function gridinterp(x, y, Z, xnew, ynew)
! 
!         real(kind=8), intent(in) :: x(:), y(:), Z(:,:), xnew, ynew
!         integer :: i, j, w, e, s, n
! 
!         w = closest_inf(xnew, x)
!         s = closest_inf(ynew, y)
!         e = w + 1
!         n = s + 1
!         gridinterp = (&
!             + (xnew-x(w)) * (ynew-y(s)) * Z(n,e) &
!             + (xnew-x(w)) * (y(n)-ynew) * Z(s,e) &
!             + (x(e)-xnew) * (ynew-y(s)) * Z(n,w) &
!             + (x(e)-xnew) * (y(n)-ynew) * Z(s,w) &
!         ) / ((x(e)-x(w)) * (y(n)-y(s)))
! 
!     end function gridinterp

    REAL(KIND=8) FUNCTION fgoutinterp(fg, q, x, y)

        REAL(KIND=8), INTENT(IN) :: q(:,:), x, y
        TYPE(fgout_grid), INTENT(IN) :: fg
        REAL(KIND=8) :: xw, xe, ys, yn
        INTEGER :: w, s

        IF (fg%x_hi<x .or. x<fg%x_low .or. fg%y_hi<y .or. y<fg%y_low) THEN
            fgoutinterp = 0.d0
        ELSE

            w = MIN(fg%mx-1, 1+INT((x-fg%x_low)/(fg%x_hi-fg%x_low)*(fg%mx-1)))
            s = MIN(fg%my-1, 1+INT((y-fg%y_low)/(fg%y_hi-fg%y_low)*(fg%my-1)))
            xw = fg%x_low + (w-1)*(fg%x_hi-fg%x_low)/(fg%mx-1)
            ys = fg%y_low + (s-1)*(fg%y_hi-fg%y_low)/(fg%my-1)
            xe = xw + (fg%x_hi-fg%x_low)/(fg%mx-1)
            yn = ys + (fg%y_hi-fg%y_low)/(fg%my-1)

            fgoutinterp = (&
                + (x-xw) * (y-ys) * q(w+1,s+1) &
                + (x-xw) * (yn-y) * q(w+1,s) &
                + (xe-x) * (y-ys) * q(w,s+1) &
                + (xe-x) * (yn-y) * q(w,s) &
            ) / ((xe-xw) * (yn-ys))
        END IF

    END FUNCTION fgoutinterp

end module helpers
