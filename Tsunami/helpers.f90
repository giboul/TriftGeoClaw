module helpers
    implicit none
    save
contains
    subroutine allocate_space(data, mthbc)
        real(kind=8), allocatable, intent(inout) :: data(:,:,:)
        integer, intent(in) :: mthbc
    
        character(len=6), dimension(4) :: sides
        character(len=32) :: ftemp
        character(len=40) :: fname
        real(kind=8) :: x, y, h, hu, hv
        integer :: unit, io, i, n
        integer :: num_cells, num_times
        logical :: res
    
        unit = 2
        sides = [character(len=6) :: "left", "right", "bottom", "top"]
    
        num_times = 0
        ftemp = "../../AVAC/_cut_output/" // trim(sides(mthbc)) // "_"
        do
            write(fname,"(A,I0.4, A4)") trim(ftemp),num_times+1,".txt"
            inquire(file=trim(fname),exist=res)
            if (.not. res) then
                exit
            end if
            num_times = num_times + 1
        end do
        close(unit)
        ! print "(A1,A,A1)", "<", trim(fname), ">"
    
        num_cells = 0
        ftemp = "../../AVAC/_cut_output/" // trim(sides(mthbc)) // "_"
        do i = 1, num_times
            n = 0
            write(fname,"(A,I0.4, A4)") trim(ftemp),i,".txt"
            print *, fname
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
        print "(A,I5)", "num_cells = ", num_cells
    
        allocate(data(num_times, num_cells, 5))
        print "(A,I10)", "size(data) = ", size(data)
   
        do i = 1, num_times
            write(fname,"(A,I0.4, A4)") trim(ftemp),i,".txt"
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
    end subroutine allocate_space

    function read_times()
        integer :: io, n, i
        integer :: unit
        real(kind=8), allocatable :: read_times(:)

        unit = 2
        open(unit, file="../../AVAC/_cut_output/timing.txt")
            n = 0
            do 
                read(unit,*,iostat=io)
                if (io /= 0) then
                    exit
                end if
                n = n + 1
            end do
            allocate(read_times(n))
        rewind(unit)
            do i = 1, n
                read(unit,*) read_times(i)
            end do
        close(unit)

    end function read_times
end module helpers
