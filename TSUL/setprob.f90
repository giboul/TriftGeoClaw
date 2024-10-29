subroutine setprob
    use helpers
    implicit none
    save

    call read_inflow_mode(inflow_mode)
    call read_avid(avid)
    call read_times(times, avid)
    if (trim(inflow_mode) == "bc") then
        call init_bc(q_avac, avid, size(times))
    else
        call init_src(q_avac, avid, size(times))
    end if

end subroutine setprob

