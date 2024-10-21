subroutine setprob
    use helpers, only: q_avac, times, init_inflows, read_times
    implicit none
    save

    call init_inflows(q_avac)
    call read_times(times)

end subroutine setprob

