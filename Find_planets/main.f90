   program bla
     implicit none
     integer :: i,j
     integer :: ep_nr                  ! Number of transits when Kepler is visible (a bit more)
     integer, parameter :: KOI_nr =           43
     character(50), dimension(KOI_nr) :: file_input_BJD
     character(50), dimension(KOI_nr) :: file_input_STARS
     character(30), dimension(KOI_nr) :: Name_KOI
     character(30), dimension(KOI_nr) :: file_data
     double precision, dimension(KOI_nr) :: Kmag

     character(40), parameter :: output_1 = "observable_eclipses"
     double precision :: orb_per       ! Best-fit orbital period of the KOIs
     double precision :: err_amp       ! Error on the timings (scatter of the residual O-C) 
     double precision :: amp_OC        ! amplitude OC diagram
     double precision, allocatable, dimension(:) :: BJD_T0      ! Array with predictions in BJD
     double precision, allocatable, dimension(:) :: BJD_T1      ! First contact in BJD
     double precision, allocatable, dimension(:) :: BJD_T4      ! Fourth contact in BJD
     double precision, allocatable, dimension(:) :: JD_T0       ! Array with predictions in JD
     double precision, allocatable, dimension(:) :: JD_T1       ! First contact in JD
     double precision, allocatable, dimension(:) :: JD_T4       ! Fourth contact in JD
     double precision, dimension(KOI_nr) :: Tdur                ! Transit duration per KOI in days
     double precision, dimension(KOI_nr) :: RA                  ! Right ascention
     double precision, dimension(KOI_nr) :: DEC                 ! Declination
     integer :: RA_hr, RA_min                                   ! To put RA in sexagesimal
     double precision :: RA_sec                                 ! To put RA in sexagesimal
     integer :: DEC_hr, DEC_min                                 ! To put DEC in sexagesimal
     double precision :: DEC_sec                                ! To put DEC in sexagesimal
     double precision, dimension(KOI_nr) :: RA_prec             ! Right ascention, precessed
     double precision, dimension(KOI_nr) :: DEC_prec            ! Declination, precessed
     double precision :: UTRising
     double precision :: UTSetting
     double precision :: UT_T0                                  ! JD to UT variable mid-transit.
     double precision :: UT_T1                                  ! JD to UT variable first contact.
     double precision :: UT_T4                                  ! JD to UT variable fourth contact.
     double precision :: StarH_T0                               ! Star altitude mid-transit time.
     double precision :: StarH_T1                               ! Star altitude first contact.
     double precision :: StarH_T4                               ! Star altitude fourth contact.
     character(10) :: cal_date_T0                               ! Calender date of JD-mid-transit time.
     character(5) :: cal_date_T1                                ! Calender date of JD-ingress.
     character(5) :: cal_date_T4                                ! Calender date of JD-egress.
     character(5) :: UT0                                        ! UT_T0 in HH:MM.
     character(5) :: UT1                                        ! UT_T1 in HH:MM.
     character(5) :: UT4                                        ! UT_T4 in HH:MM.

     ! Change here to addapt to the different observatories we have for TESS.
     integer, parameter :: obs_nr = 1
     integer :: ii
     double precision, dimension(1), parameter :: LAT = (/28.29833/)
     double precision, dimension(1), parameter :: LON = (/343.4906/)
     double precision, dimension(1), parameter :: H = (/2346/)

     integer, dimension(KOI_nr) :: flag
     double precision :: aux
     double precision :: aux1, aux3, aux4, aux5, aux6
     character(30) :: aux2, aux7

     open(901, file = output_1, status = "unknown")
     
     !Create a loop reading on input data that fulfills this. RA/DEC
     !have to be in concordance with TESS CVZ. KOI_nr is the variable
     !that assigns how many targets are in the input list.

     open(501, file = 'file_input_STARS', status = "unknown")
     do i=1, KOI_nr
        read(501,*)file_data(i)
     enddo
     close(501)

     do ii=1, obs_nr
     write(*,*)ii, LAT(ii), LON(ii), H(ii)

     do i=1, KOI_nr
        open(502, file = file_data(i), status = "unknown")
        read(502,*)aux1, aux2, aux3, aux4, aux5, aux6, aux7
        close(502)
        flag(i) = aux1
        Name_KOI(i) = aux2
        Tdur(i) = aux3/24.d0
        RA(i) = aux4
        DEC(i) = aux5
        Kmag(i) = aux6
        file_input_BJD(i) = aux7
     enddo

     do i=1, KOI_nr
 
        open(45, file = file_input_BJD(i), status = "unknown")
        ep_nr = 1
        do while(.true.)
           read(45,*,end=46)
           ep_nr = ep_nr + 1
        enddo
 46     close(45)
        ep_nr = ep_nr - 1
 
        allocate( BJD_T0(ep_nr), BJD_T1(ep_nr), BJD_T4(ep_nr) )
        allocate( JD_T0(ep_nr), JD_T1(ep_nr), JD_T4(ep_nr) )
        open(45, file = file_input_BJD(i), status = "unknown")
        do j=1, ep_nr
           read(45,*)aux, BJD_T0(j)
           BJD_T1(j) = BJD_T0(j) - (Tdur(i)/2.d0)
           BJD_T4(j) = BJD_T0(j) + (Tdur(i)/2.d0)
        enddo
        close(45)
 
!        call BJDtoJD(ep_nr, BJD_T0, JD_T0, BJD_T1, JD_T1, BJD_T4, JD_T4, RA(i), DEC(i))
	JD_T0 = BJD_T0
	JD_T1 = BJD_T1
	JD_T4 = BJD_T4 

        call RA_DEC_sgs(RA(i), DEC(i), RA_hr, RA_min, RA_sec, DEC_hr, DEC_min, DEC_sec)

        do j=1, ep_nr
           call precession(JD_T0(j), RA(i), DEC(i), RA_prec(i), DEC_prec(i))
           call day_or_night(LON(ii), LAT(ii), JD_T0(j), UTRising, UTSetting)
 
           UT_T0 = JD_T0(j) - int(JD_T0(j)) - 0.5d0
           if (UT_T0.lt.0.d0) then
              UT_T0 = UT_T0 + 1.d0
           endif
 
           UT_T1 = JD_T1(j) - int(JD_T1(j)) - 0.5d0
           if (UT_T1.lt.0.d0) then
              UT_T1 = UT_T1 + 1.d0
           endif
 
           UT_T4 = JD_T4(j) - int(JD_T4(j)) - 0.5d0
           if (UT_T4.lt.0.d0) then
              UT_T4 = UT_T4 + 1.d0
           endif
 
           call jd2date_2(JD_T0(j), cal_date_T0, JD_T1(j), cal_date_T1,&
                JD_T4(j), cal_date_T4)
 
           call ut2hhmm(UT_T0, UT_T1, UT_T4, UT0, UT1, UT4)
 
           if (UTSetting.lt.UTRising) then
              if ( (UT_T1.gt.UTSetting).and.(UT_T4.lt.UTRising).and.&
                   (UT_T4.gt.UT_T1) ) then
 
              call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                   StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
 
              if ( StarH_T1.gt.20.d0.and.StarH_T0.gt.20.d0.and.StarH_T4.gt.20.d0 ) then
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec),RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec),DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "C",flag(i)
 
 911             format(a10,2x,i2.2,a1,i2.2,a1,i2.2,f0.2,2x,&
                      i3.2,a1,i2.2,a1,i2.2,f0.2,2x,&
                      f6.3,2x,f12.4,&
                      5x,a5,1x,a5,1x,"(",i2,"°)",&
                      5x,a10,1x,a5,1x,"(",i2,"°)",&
                      5x,a5,1x,a5,1x,"(",i2,"°)",&
                      2x,a1,2x,i2)
 912             format(a10,2x,i2.2,a1,i2.2,a1,i2.2,f0.2,2x,&
                      i3.2,a1,i2.2,a1,i2.2,f0.2,2x,&
                      f6.3,2x,f12.4,&
                      5x,a5,1x,a5,1x,"(",i2,"°)",&
                      5x,a10,1x,a5,1x,"(",i2,"°)",&
                      5x,a5,1x,a5,1x,"(",i3,"°)",&
                      2x,a1,2x,i2)
 913             format(a10,2x,i2.2,a1,i2.2,a1,i2.2,f0.2,2x,&
                      i3.2,a1,i2.2,a1,i2.2,f0.2,2x,&
                      f6.3,2x,f12.4,&
                      5x,a5,1x,a5,1x,"(",i2,"°)",&
                      5x,a10,1x,a5,1x,"(",i2,"°)",&
                      5x,a5,1x,a5,1x,"(",i2,"°)",&
                      2x,a1,2x,i2)
              endif
 
           elseif ( (UT_T1.lt.UTRising).and.(UT_T1.ge.UTSetting).and.&
                (UT_T4.gt.UTRising) ) then
 
              call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                   StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
              if ( StarH_T1.gt.20.d0.and.StarH_T0.gt.20.d0 ) then
                 if ( StarH_T4.le.-10 ) then
                 write(901,912)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "I",flag(i)
 
                 else
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "I",flag(i)
 
                 endif
              endif
 
           elseif ( (UT_T1.lt.UTSetting).and.(UT_T4.gt.UTSetting).and.&
                (UT_T4.lt.UTRising) ) then
 
              call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                   StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
              if ( StarH_T0.gt.20.d0.and.StarH_T4.gt.20.d0 ) then
                 if ( StarH_T1.le.-10 ) then
                 write(901,913)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "E",flag(i)
 
                 else
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "E",flag(i)
 
                 endif
              endif
 
           elseif ( (UT_T4.gt.UTSetting).and.(UT_T4.lt.UTRising).and.&
                (UT_T1.lt.UTRising).and.(UT_T1.gt.UTSetting).and.(UT_T1.gt.UT_T4) ) then  
 
              call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                   StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
              if ( StarH_T1.gt.30.d0.and.StarH_T4.gt.30.d0 ) then
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "S",flag(i)
 
              endif
           endif
           elseif (UTRising.lt.UTSetting) then
 
              if ( (UT_T1.gt.UTSetting.and.UT_T4.gt.UTSetting.and.UT_T4.gt.UT_T1).or.&
                   (UT_T1.lt.UTRising.and.UT_T4.lt.UTRising.and.UT_T4.gt.UT_T1).or.&
                   (UT_T1.gt.UTSetting.and.UT_T4.lt.UTRising.and.UT_T1.gt.UT_T4) ) then
                 call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                      StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
                 
                 if ( StarH_T1.gt.20.d0.and.StarH_T0.gt.20.d0.and.StarH_T4.gt.20.d0 ) then
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "C",flag(i)
 
                 endif
              elseif ( (UT_T1.lt.UTRising.and.UT_T4.gt.UTRising.and.UT_T4.lt.UTSetting).or.&
                   (UT_T1.gt.UTSetting.and.UT_T4.gt.UTRising.and.UT_T1.gt.UT_T4) ) then
                 call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                      StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
                 
                 if ( StarH_T1.gt.20.d0.and.StarH_T0.gt.20.d0 ) then
                    if ( StarH_T4.le.-10 ) then
                 write(901,912)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "I",flag(i)
 
                    else
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "I",flag(i)
 
                    endif
                 endif
              elseif ( (UT_T4.gt.UTSetting.and.UT_T1.lt.UTSetting.and.UT_T1.gt.UTRising).or.&
                   (UT_T4.lt.UTRising.and.UT_T1.gt.UTRising.and.UT_T1.lt.UTSetting.and.&
                   UT_T1.gt.UT_T4) ) then
                 call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                      StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
                
                 if ( StarH_T0.gt.20.d0.and.StarH_T4.gt.20.d0 ) then
                    if ( StarH_T1.le.-10 ) then
                 write(901,913)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "E",flag(i)
 
                    else
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "E",flag(i)
 
                    endif
                 endif
              elseif ( (UT_T1.lt.UTRising.and.UT_T4.gt.UTSetting.and.UT_T4.gt.UT_T1) ) then
                 call star_altitude(LAT(ii), LON(ii), H(ii), JD_T1(j), JD_T0(j), JD_T4(j),&
                      StarH_T1, StarH_T0, StarH_T4, RA_prec(i), DEC_prec(i))
                 
                 if ( StarH_T1.gt.30.d0.and.StarH_T4.gt.30.d0 ) then
                 write(901,911)Name_KOI(i), RA_hr,":",RA_min,":",int(RA_sec), RA_sec - int(RA_sec),&
                      DEC_hr,":",DEC_min,":",int(DEC_sec), DEC_sec - int(DEC_sec),&
                      Kmag(i), BJD_T0(j),&
                      cal_date_T1, UT1, nint(StarH_T1),&
                      cal_date_T0, UT0, nint(StarH_T0),&
                      cal_date_T4, UT4, nint(StarH_T4),&
                      "S",flag(i)
 
                 endif
              endif
           endif
        enddo
 
        deallocate( BJD_T0, BJD_T1, BJD_T4 )
        deallocate( JD_T0, JD_T1, JD_T4 )
 
     enddo
     enddo
     close(901)
   end program bla
