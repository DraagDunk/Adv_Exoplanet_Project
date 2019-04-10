  subroutine precession(Toi, StarRAo, StarDECo, StarRA, StarDEC)

    implicit none

    !Subroutine variables:
    double precision :: Toi
    double precision :: StarRAo
    double precision :: StarDECo
    double precision :: StarRA
    double precision :: StarDEC

    ! Internal variables:
    double precision, parameter :: pi = 4.d0*atan(1.d0)
    double precision :: StarRAm
    double precision :: StarDECm
    double precision, parameter :: JDo = 2451545.d0
    double precision, parameter :: Yearo = 36525.d0
    double precision :: T
    double precision :: M
    double precision :: N
    double precision, parameter :: M1 = 1.2811556689d0
    double precision, parameter :: M2 = 0.00038655131d0
    double precision, parameter :: M3 = 0.000010079625d0
    double precision, parameter :: M4 = 9.60194e-9
    double precision, parameter :: M5 = 1.68806e-10
    double precision, parameter :: N1 = 0.5567199731d0
    double precision, parameter :: N2 = 0.00011930372d0
    double precision, parameter :: N3 = 0.000011617400d0
    double precision, parameter :: N4 = 1.96917e-9
    double precision, parameter :: N5 = 3.5389e-11

    T = (Toi - JDo)/Yearo

    M = M1*T + M2*T**2 + M3*T**3 - M4*T**4 - M5*T**5
    N = N1*T - N2*T**2 - N3*T**3 - N4*T**4 - N5*T**5

    StarRAm = StarRAo + 0.5d0*(M/15.d0 + N/15.d0*sin(pi*15.d0*StarRAo/180.d0)*&
         tan(pi*StarDECo/180.d0))
    StarDECm = StarDECo + 0.5d0*N*cos(pi*15.d0*StarRAm/180.d0)

    StarRA = StarRAo + M/15.d0 + N/15.d0*sin(pi*15.d0*StarRAm/180.d0)*&
         tan(pi*StarDECm/180.d0)
    StarDEC = StarDECo + N*cos(pi*15.d0*StarRAm/180.d0)

    return

  end subroutine precession
