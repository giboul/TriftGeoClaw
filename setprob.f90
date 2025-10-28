subroutine setprob
    use helpers
    implicit none
    save

    call read_data()
    if (trim(mode) == "bc") then
        call init_bc()
    else if (trim(mode) == "src") then
        call init_src_fgout()
    end if

end subroutine setprob
