subroutine setprob
    use helpers
    implicit none
    save

    call read_data(avid, inflow_mode, mx, my)
    call read_times(times)
    if (trim(inflow_mode) == "bc") then
        call init_bc(q_avac, size(times))
    else if (trim(inflow_mode) == "src") then
        call init_src(q_avac, size(times))
    end if

end subroutine setprob

