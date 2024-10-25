subroutine setprob
    use helpers
    implicit none
    save

    call read_avid(avid)
    call read_times(times, avid)
    call init_inflows(q_avac, avid, size(times))

end subroutine setprob

