  subroutine jd2date_2(JD_T0, cal_date_T0, JD_T1, cal_date_T1,&
               JD_T4, cal_date_T4)

    integer, dimension(12), parameter :: MonthDays = (/ 31, 28, 31, 30, 31, &
         30, 31, 31, 30, 31, 30, 31 /)
    integer, dimension(12), parameter :: MonthDaysLeap = (/ 31, 29, 31, 30, &
         31, 30, 31, 31, 30, 31, 30, 31 /)
    integer, dimension(12) :: MonthDaysFinal

    integer, dimension(10), parameter :: days_year_2013 = (/ 365, 365, 365, &
         366, 365, 365, 365, 366, 365, 365 /)

    double precision, parameter :: Jul0 = 2456293.5d0 ! 01/01/2013, 00 hs
    double precision, parameter :: Year0 = 365.d0

    integer :: YearNr
    integer :: MonthNr
    integer :: DayNr
    integer :: MoonAge

    integer :: Day
    double precision :: Year
    integer :: leapyear

    integer :: j, aux

    double precision :: Toi

    double precision :: JD_T0
    double precision :: JD_T1
    double precision :: JD_T4

    character(10) :: cal_date_T0
    character(5) :: cal_date_T1
    character(5) :: cal_date_T4
    character(2) :: char1_aux
    character(2) :: char2_aux
    character(4) :: char3_aux

    ! ----------------------------------------------------- !
    ! Convert JD-mid-transit into DD.MM                     !
    ! ----------------------------------------------------- !

    Toi = JD_T0
    YearNr = int((Toi - Jul0)/Year0)
    Year = (Toi - Jul0)/Year0
    
    DayNr = 0
    do j =1, YearNr
       DayNr = DayNr + days_year_2013(j)
    enddo
    DayNr = floor(Toi - (Jul0 + DayNr))

    leapyear = mod(YearNr + 2013, 4)
    if (leapyear.le.0) then
       MonthDaysFinal = MonthDaysLeap
    elseif (leapyear.gt.0) then
       MonthDaysFinal = MonthDays
    endif

    aux = 0
    MonthNr = 0
    do while (aux.le.DayNr)
       aux = aux + MonthDaysFinal(MonthNr + 1)
       MonthNr = MonthNr + 1
    enddo
    MonthNr = MonthNr - 1
    
    aux = 0
    do j=1, MonthNr
       aux = aux + MonthDaysFinal(j)
    enddo

    write(char1_aux,"(i2.2)") DayNr - aux + 1
    write(char2_aux,"(i2.2)") MonthNr+1
    write(char3_aux,"(i4.4)") YearNr + 2013
    cal_date_T0 = char1_aux//'.'//char2_aux//'.'//char3_aux

    ! ----------------------------------------------------- !
    ! Convert JD-ingress into DD.MM                         !
    ! ----------------------------------------------------- !

    Toi = JD_T1
    YearNr = int((Toi - Jul0)/Year0)
    Year = (Toi - Jul0)/Year0

    DayNr = 0
    do j =1, YearNr
       DayNr = DayNr + days_year_2013(j)
    enddo
    DayNr = floor(Toi - (Jul0 + DayNr))

    aux = 0
    MonthNr = 0
    do while (aux.le.DayNr)
       aux = aux + MonthDaysFinal(MonthNr + 1)
       MonthNr = MonthNr + 1
    enddo
    MonthNr = MonthNr - 1
    
    aux = 0
    do j=1, MonthNr
       aux = aux + MonthDaysFinal(j)
    enddo

    write(char1_aux,"(i2.2)") DayNr - aux + 1
    write(char2_aux,"(i2.2)") MonthNr+1
    cal_date_T1 = char1_aux//'.'//char2_aux


    ! ----------------------------------------------------- !
    ! Convert JD-mid-egress into DD.MM                      !
    ! ----------------------------------------------------- !

    Toi = JD_T4
    YearNr = int((Toi - Jul0)/Year0)
    Year = (Toi - Jul0)/Year0

    DayNr = 0
    do j =1, YearNr
       DayNr = DayNr + days_year_2013(j)
    enddo
    DayNr = floor(Toi - (Jul0 + DayNr))

    aux = 0
    MonthNr = 0
    do while (aux.le.DayNr)
       aux = aux + MonthDaysFinal(MonthNr + 1)
       MonthNr = MonthNr + 1
    enddo
    MonthNr = MonthNr - 1
    
    aux = 0
    do j=1, MonthNr
       aux = aux + MonthDaysFinal(j)
    enddo

    write(char1_aux,"(i2.2)") DayNr - aux + 1
    write(char2_aux,"(i2.2)") MonthNr+1
    cal_date_T4 = char1_aux//'.'//char2_aux

    return
  end subroutine jd2date_2
