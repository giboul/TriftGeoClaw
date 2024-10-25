module helpers
    implicit none
    save

    real(kind=8), allocatable :: q_avac(:,:,:)
    real(kind=8), allocatable :: times(:)
    real(kind=8) :: damping = 500.d0 / 1e3
    character(len=4) :: avid

contains

    subroutine read_avid(avid)
        character(len=4), intent(out) :: avid
        open(2, file="../avac.data", status='old')
            read(2,*) avid
        close(2)
        if (avid == "None") then
            avid = ""
        end if
    end subroutine read_avid

    integer function closest(value, array)
        integer :: i
        real(kind=8), intent(in) :: value
        real(kind=8), dimension(:), intent(in) :: array

        closest = 1
        do i = 1, size(array)
            if (abs(value-array(i))<abs(value-array(closest))) then
                closest = i
            end if
        end do
    end function closest

    real(kind=8) function interp(newx, x, y)
        real(kind=8), intent(in) :: x(:), y(:), newx
        real(kind=8) :: xp, yp
        integer :: ix

        ix = closest(newx, x) ! closest_back
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


    subroutine init_inflows(data, avid, num_times)
        real(kind=8), allocatable, intent(inout) :: data(:,:,:)
        character(len=4), intent(in) :: avid
        integer, intent(in) :: num_times
        character(len=255) :: ftemp
        character(len=255) :: fname
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n, mthbc
        integer :: num_cells
        logical :: res
    
        ftemp = "../../AVAC/_cut_output"//trim(avid)//"/cut"
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
    
        allocate(data(num_times, num_cells, 5))
        print "(A,I10)", "size(data)    = ", size(data)
        print "(A,I10)", "size(data, 1) = ", size(data, 1)
        print "(A,I10)", "size(data, 2) = ", size(data, 2)
        print "(A,I10)", "size(data, 3) = ", size(data, 3)
   
        do i = 1, num_times
            write(fname,"(A,I0.4, A4)") trim(ftemp),i-1,".txt"
            open(unit, file=fname, status="old")
            do n = 1, num_cells
                read(unit,*, iostat=io) x, y, h, hu, hv
                data(i, n, 1) = x
                data(i, n, 2) = y
                data(i, n, 3) = h
                data(i, n, 4) = hu
                data(i, n, 5) = hv
            end do
            close(unit)
        end do
    end subroutine init_inflows

    subroutine read_times(times, avid)
        character(len=4), intent(in) :: avid
        character(len=255) :: fname
        integer :: io, n, i
        integer :: unit
        real(kind=8), allocatable :: times(:)

        unit = 2
        fname = "../../AVAC/_cut_output"//trim(avid)//"/timing.txt"
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

end module helpers
