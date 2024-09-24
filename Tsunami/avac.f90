module avac
    implicit none
    save

    real(kind=8) :: h0, hu0, hv0

    OPEN(unit=1,file='bc.data')
        READ(1,*) h0
        READ(1,*) hu0
        READ(1,*) hv0
    CLOSE(1)
    print *,"h0, hu0, hv0 = ", h0, hu0, hv0

end module avac
