subroutine setprob
    use helpers
    implicit none
    save

    real(kind=8) :: damping = 500.d0 / 1e3

    call init_inflows(q_avac)
    call read_times(times)

end subroutine setprob

