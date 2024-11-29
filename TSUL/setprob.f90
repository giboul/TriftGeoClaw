subroutine setprob
    use helpers
    implicit none
    save

    call read_data(avid, inflow_mode)
    call read_times(times)
    if (trim(inflow_mode) == "bc") then
        call init_bc(q_avac, size(times))
    else if (trim(inflow_mode) == "src") then
        call init_src_fgout_bin(avid, size(times), q_avac, FGOUT_fgrids, FGOUT_num_grids)
    end if

end subroutine setprob
