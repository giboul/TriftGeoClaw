module helpers
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

        IF (mode == "None") then
            mode = "bc"
        END IF
        IF (avid == "None") then
            avid = ""
        END IF

    END SUBROUTINE read_data

    subroutine read_times(times, avid)
        character(len=4), intent(in) :: avid
        character(len=255) :: fname
        integer :: io, n, i
        integer :: unit
        real(kind=8), allocatable :: times(:)

        unit = 2
        fname = "../_inflows"//trim(avid)//"/timing.txt"
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
        integer :: i
        real(kind=8), intent(in) :: value
        real(kind=8), dimension(:), intent(in) :: array
        closest = MINLOC(ABS(array-value), DIM=1)
    end function closest

    INTEGER FUNCTION closest_inf(value, array)

        REAL(KIND=8), INTENT(IN) :: value
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: array
        LOGICAL, DIMENSION(SIZE(array)) :: mask

        closest_inf = MINLOC(ABS(array-value), MASK=array<value, DIM=1)
        IF (closest_inf < 1) THEN
            closest_inf = MINLOC(array, DIM=1)
        END IF
        closest_inf = MIN(SIZE(array)-1, closest_inf)

    END FUNCTION closest_inf

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

    subroutine init_src(data, avid, num_times)
        real(kind=8), allocatable, intent(inout) :: data(:,:,:,:)
        character(len=4), intent(in) :: avid
        integer, intent(in) :: num_times
        character(len=255) :: ftemp
        character(len=255) :: fname
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n, mthbc
        integer :: num_cells
        logical :: res
    
        ftemp = "../_inflows"//trim(avid)//"/cut"
        num_cells = 0
        do i = 1, num_times
            n = 0
            write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".txt"
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
    
        allocate(data(1, num_times, num_cells, 5))
        print "(A,I10)", "size(data)    = ", size(data)
        print "(A,I10)", "size(data, 1) = ", size(data, 1)
        print "(A,I10)", "size(data, 2) = ", size(data, 2)
        print "(A,I10)", "size(data, 3) = ", size(data, 3)
        print "(A,I10)", "size(data, 3) = ", size(data, 4)
   
        do i = 1, num_times
            write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".txt"
            open(unit, file=fname, status="old")
            do n = 1, num_cells
                read(unit,*, iostat=io) x, y, h, hu, hv
                data(1, i, n, 1) = x
                data(1, i, n, 2) = y
                data(1, i, n, 3) = h
                data(1, i, n, 4) = hu
                data(1, i, n, 5) = hv
            end do
            close(unit)
        end do
    end subroutine init_src

    subroutine init_bc(data, avid, num_times)
        real(kind=8), allocatable, intent(inout) :: data(:,:,:,:)
        integer, intent(in) :: num_times
        character(len=6), dimension(4) :: sides
        character(len=255) :: fdir, ftemp
        character(len=255) :: fname
        character(len=4) :: avid
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n, mthbc
        integer :: num_cells
        logical :: res
    
        unit = 2
        sides = [character(len=6) :: "left", "right", "bottom", "top"]
 
        open(unit, file="../avac.data", status='old')
            read(unit,*) avid
        close(unit)
        if (avid == "None") then
            avid = "    "
        end if
        print *, "Avalanche id #", trim(avid)

        fdir = "../_inflows"//trim(avid)//"/"
        ftemp = trim(fdir)//trim(sides(1))//"_"
        num_cells = 0
        do mthbc = 1, 4
            ftemp = trim(fdir)//trim(sides(mthbc))//"_"
            do i = 1, num_times
                n = 0
                write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".txt"
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
            ftemp = trim(fdir)//trim(sides(mthbc))//"_"
            do i = 1, num_times
                write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".txt"
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
