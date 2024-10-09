subroutine setprob
    implicit none
    save
    ! common /params/ damping, q_left, q_top, q_right, q_bottom

    real(kind=8) :: damping = 500.d0 / 1e3
    real(kind=8), allocatable :: q_left(:,:,:), q_right(:,:,:)
    real(kind=8), allocatable :: q_top(:,:,:), q_bottom(:,:,:)

    call allocate_space(q_left, 1)
    call allocate_space(q_right, 2)
    call allocate_space(q_bottom, 3)
    call allocate_space(q_top, 4)
    print *, "setprob end"

end subroutine setprob

subroutine allocate_space(data, mthbc)
    real(kind=8), allocatable, intent(in) :: data(:,:,:)
    integer, intent(in) :: mthbc

    character(len=6), dimension(4) :: sides
    character(len=32) :: ftemp
    character(len=40) :: fname
    integer :: unit, io, i, n
    integer :: num_cells, num_times
    logical :: res
    real(kind=8) :: x, y, h, hu, hv

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


! function read_q_t(filename)
!     real(kind=8), allocatable :: read_q(:,:), data(:,:)
!     integer :: n, i, io, unit
!     real(kind=8) :: x, y, h, hu, hv
!     character(*), intent(in) :: filename
!     unit = 2
!     n = 0
!     open(unit, file=filename, status="replace")
!         do
!             read(unit,*,iostat=io) t, x, h, hu, hv
!             if (io /= 0) then
!                 exit
!             end if
!             n = n + 1
!         end do
!         allocate(data(n, 5))
!         rewind(unit)
!         do i = 1, n
!             read(unit,*,iostat=io) t, x
!             data(i,1) = x
!             data(i,2) = y
!             data(i,3) = h
!             data(i,4) = hu
!             data(i,5) = hv
!         end do
!     close(unit)
!     allocate(read_q(shape(data)))
!     read_q = data
! end function read_q_t
! 
! function locate(xx, x)
!     !from fortran recipes
!     implicit none
!     real(kind=8), dimension(:), intent(in) :: xx
!     real(kind=8) :: x
!     integer :: locate
!     !Given an array xx(1:N), and given a value x, returns a value j such that x is betwe>
!     !xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
!     !j = N is returned to indicate that x is out of range.
!     integer :: n, jl, jm, ju
!     Logical :: ascnd
!     n = size(xx)
!     ascnd = (xx(n) >= xx(1)) !True if ascending order of table, false otherwise.
!     jl = 0                     !Initialize lower
!     ju = n + 1                   !and upper limits.
!     do
!         if (ju - jl <= 1) exit !Repeat until this condition is satisfied.
!         jm = (ju + jl)/2 !Compute a midpoint,
!         if (ascnd .eqv. (x >= xx(jm))) then
!            jl = jm !and replace either the lower limit
!         else
!            ju = jm !or the upper limit, as appropriate.
!         end if
!     end do
!     if (x == xx(1)) then !Then set the output, being careful with the endpoints.
!         locate = 1
!     else if (x == xx(n)) then
!         locate = n - 1
!     else
!         locate = jl
!     end if
! end function locate
! 
! 
! real(kind=8) function interp1d(xt, yt, x)
!     real(kind=8), intent(in) :: xt, yt
!     real(kind=8) :: x
!     interp1d = x
! end function interp1d
