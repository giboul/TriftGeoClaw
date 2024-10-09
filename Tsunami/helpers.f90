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
    
        allocate(data(num_times, num_cells, 6))
        print "(A,I10)", size(data)
    
    end subroutine allocate_space
end module helpers
