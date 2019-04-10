  subroutine RA_DEC_sgs(RA, DEC,&
       RA_hr, RA_min, RA_sec, DEC_hr, DEC_min, DEC_sec)

    implicit none
    double precision :: RA
    double precision :: DEC
    integer :: RA_hr
    integer :: RA_min
    double precision :: RA_sec
    integer :: DEC_hr
    integer :: DEC_min
    double precision :: DEC_sec

    RA_hr = int(RA)
    DEC_hr = int(DEC)

    RA_min = int((RA - int(RA))*60.d0)
    DEC_min = abs(int((DEC - int(DEC))*60.d0))

    RA_sec = (RA*60.d0 - int(RA*60.d0))*60.d0
    DEC_sec = abs((DEC*60.d0 - int(DEC*60.d0))*60.d0)

    return
  end subroutine RA_DEC_sgs
