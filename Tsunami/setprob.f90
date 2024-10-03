subroutine setprob
    implicit none
    real(kind=8) :: damping
    character(len=255) :: path
    common /params/ damping

    damping = 500 / 1e3 ! snow density over water density

end subroutine setprob
