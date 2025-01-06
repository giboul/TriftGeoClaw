subroutine setprob
    use helpers
    implicit none
    save

    call read_data()
    ! call read_times()
    if (trim(inflow_mode) == "bc") then
        call init_bc()
    else if (trim(inflow_mode) == "src") then
        call init_src_fgout_bin()
    end if

end subroutine setprob
