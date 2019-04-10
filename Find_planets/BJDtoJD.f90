  subroutine BJDtoJD(ep_nr, BJD_T0, JD_T0, BJD_T1, JD_T1, BJD_T4, JD_T4, ra, dec)

    implicit none
    
    integer :: i
    integer :: ep_nr
    
    double precision, dimension(ep_nr) :: BJD_T0
    double precision, dimension(ep_nr) :: JD_T0
    double precision, dimension(ep_nr) :: BJD_T1
    double precision, dimension(ep_nr) :: JD_T1
    double precision, dimension(ep_nr) :: BJD_T4
    double precision, dimension(ep_nr) :: JD_T4
    double precision :: ra
    double precision :: dec

    character(len=8), parameter :: script = './dia.sh'   
    character(len=4) :: input_file
    character(len=3) :: output_file
    character(len=100) :: comando

    ! ------------------------------------------------------------- !
    ! Writing basic information in the script.                      !
    ! ------------------------------------------------------------- !

    ! The if statement should take care of any possible DEC and RA, considering now K2.

    open(75, file = "dia.sh", status = "unknown")
    write(75,'(A,a300)')"#!/bin/bash"
    write(75,'(A,a300)')"touch $2"
    write(75,'(A,a300)')"for jd in `cat $1`"
    write(75,'(A,a300)')"do"
    write(75,760)"  curl --referer http://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html -d jds=$jd \\"
    if (int(ra).le.9) then
       write(75,7611)"  -d ra=",ra," \\"
    elseif (int(ra).ge.10) then
       write(75,7612)"  -d ra=",ra," \\"
    endif
    if (int(dec).le.9.and.int(dec).ge.0) then
       write(75,764)"  -d raunits=hours -d dec=",dec," \\"
    elseif (int(dec).ge.10) then
       write(75,762)"  -d raunits=hours -d dec=",dec," \\"
    elseif (int(dec).lt.0.and.int(dec).ge.-9) then
       write(75,762)"  -d raunits=hours -d dec=",dec," \\"
    elseif (int(dec).le.-10) then
       write(75,765)"  -d raunits=hours -d dec=",dec," \\"
    endif
    write(75,763)"  -d spaceobs=none ""http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.php"" \\"
    write(75,*)"  2>/dev/null | sed -n '/^[0-9]/ p' >> $2"
    write(75,'(A,a4)')"done"
    close(75)
760 format(a90)
7611 format(a8,f10.8,a2)
7612 format(a8,f11.8,a2)
762 format(a26,f11.8,a2)
764 format(a26,f10.8,a2)
765 format(a26,f12.8,a2)
763 format(a82)

    input_file = "BJDs"
    output_file = "JDs"

    ! ------------------------------------------------------------- !
    ! BJD to JD for BJD_T0.                                         !
    ! ------------------------------------------------------------- !

    open(45, file = input_file, status="unknown")

    do i=1, ep_nr
       write(45,*)BJD_T0(i)
    enddo
    close(45)

    open(45, file=input_file, status="old")
    open(55, file = output_file, status="unknown")
    comando = script // ' ' // input_file // ' ' // output_file
    
    call system(comando)

    close(45)
    close(55)
    
    open(55, file = output_file, status="unknown")
    do i=1, ep_nr
       read(55,*)JD_T0(i)
    enddo
    close(55)

    call system('rm BJDs')
    call system('rm JDs')

    ! ------------------------------------------------------------- !
    ! BJD to JD for BJD_T1.                                         !
    ! ------------------------------------------------------------- !

    open(45, file = input_file, status="unknown")

    do i=1, ep_nr
       write(45,*)BJD_T1(i)
    enddo
    close(45)

    open(45, file=input_file, status="old")
    open(55, file = output_file, status="unknown")
    comando = script // ' ' // input_file // ' ' // output_file
    
    call system(comando)

    close(45)
    close(55)
    
    open(55, file = output_file, status="unknown")
    do i=1, ep_nr
       read(55,*)JD_T1(i)
    enddo
    close(55)

    call system('rm BJDs')
    call system('rm JDs')

    ! ------------------------------------------------------------- !
    ! BJD to JD for BJD_T0.                                         !
    ! ------------------------------------------------------------- !

    open(45, file = input_file, status="unknown")

    do i=1, ep_nr
       write(45,*)BJD_T4(i)
    enddo
    close(45)

    open(45, file=input_file, status="old")
    open(55, file = output_file, status="unknown")
    comando = script // ' ' // input_file // ' ' // output_file
    
    call system(comando)

    close(45)
    close(55)
    
    open(55, file = output_file, status="unknown")
    do i=1, ep_nr
       read(55,*)JD_T4(i)
    enddo
    close(55)

    call system('rm BJDs')
    call system('rm JDs')

    return
  end subroutine BJDtoJD
