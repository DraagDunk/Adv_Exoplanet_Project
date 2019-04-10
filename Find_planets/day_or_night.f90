subroutine day_or_night(lamb, phi, Toi, UTRising, UTSetting)
  
  implicit none

  ! Subroutine variables.
  double precision :: lamb
  double precision :: phi
  double precision :: UTRising 
  double precision :: UTSetting
  double precision :: Toi

  ! Internal variables:
  integer :: DiaAnno
  integer, parameter :: Ene = 31
  integer, parameter :: Feb = 28
  integer, parameter :: Mar = 21 ! solsticio de invierno.
  double precision, parameter :: JDo = 2451545.0d0      ! CE 2000 January 01 12:00:00.0 UT
  double precision, parameter :: anno = 365.25d0
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: eps = 23.439d0*pi/180.d0
  double precision :: lngHour
  double precision :: tRising, tSetting
  double precision :: MRising, MSetting
  double precision :: LRising, LSetting
  double precision :: RARising, RASetting
  double precision :: Lquadrant, RAquadrant
  double precision :: sinDecRising, cosDecRising
  double precision :: sinDecSetting, cosDecSetting
  double precision :: cosHRising, cosHSetting
  double precision :: HRising, HSetting
  double precision :: T_Rising, T_Setting 
  double precision, parameter :: zenith = (90.d0 + 16.d0)*pi/180.d0
  double precision :: lambRAD
  double precision :: phiRAD

  lambRAD = lamb*(pi/180.d0)
  phiRAD = phi*(pi/180.d0)

  ! -----------------------------------------------------------------!
  ! Calculateion of sunrise and sunset as a function of the UT,      !
  ! which will be compared to the JD - 12 hs.                        !
  ! -----------------------------------------------------------------!

  ! First calculate the day of the year
  DiaAnno = int(((Toi - JDo)/anno-int((Toi - JDo)/anno))*anno)

  ! -----------------------------------------------------------------!
  ! Convert the longitude to hour value and calculate an approximate !
  ! time. lngHour in hours.                                          !
  ! -----------------------------------------------------------------!

  lngHour = 180.d0/pi*lambRAD/15.d0                         
  tRising = float(DiaAnno) + ((6.d0 - lngHour)/24.d0)
  tSetting = float(DiaAnno) + ((18.d0 - lngHour)/24.d0)

  ! -----------------------------------------------------------------!  
  ! Calculate the Sun's mean anomaly.                                !
  ! -----------------------------------------------------------------!

  MRising = (0.9856d0 * tRising) - 3.289d0
  MSetting = (0.9856d0 * tSetting) - 3.289d0

  ! -----------------------------------------------------------------! 
  ! Calculate the Sun's true longitude.                              !
  ! -----------------------------------------------------------------!

  LRising = MRising + (1.916d0 * sin(pi/180.d0*MRising)) + (0.020d0 * &
       sin(2.d0 * pi/180.d0*MRising)) + 282.634d0
  LSetting = MSetting + (1.916d0 * sin(pi/180.d0*MSetting)) + (0.020d0 * &
       sin(2.d0 * pi/180.d0*MSetting)) + 282.634d0

  if (LRising.le.0.d0) then
     do while (LRising.lt.0.d0)
        LRising = LRising + 360.d0
     enddo
  else
     do while (LRising.gt.360.d0)
        LRising = LRising - 360.d0
     enddo
  endif

  if (LSetting.le.0.d0) then
     do while (LSetting.lt.0.d0)
        LSetting = LSetting + 360.d0
     enddo
  else
     do while (LSetting.gt.360.d0)
        LSetting = LSetting - 360.d0
     enddo
  endif

  ! -----------------------------------------------------------------!
  ! NOTE: L potentially needs to be adjusted into the range [0,360)  !
  ! by adding/subtracting 360                                        !
  ! -----------------------------------------------------------------!

  ! -----------------------------------------------------------------!
  ! Calculate the Sun's right ascension.                             !
  ! -----------------------------------------------------------------!

  RARising = 180.d0/pi*atan(0.91764d0 * tan(pi/180.d0*LRising))
  RASetting = 180.d0/pi*atan(0.91764d0 * tan(pi/180.d0*LSetting))

  if (RARising.le.0.d0) then
     do while (RARising.lt.0.d0)
        RARising = RARising + 360.d0
     enddo
  else
     do while (RARising.gt.360.d0)
        RARising = RARising - 360.d0
     enddo
  endif

  if (RASetting.le.0.d0) then
     do while (RASetting.lt.0.d0)
        RASetting = RASetting + 360.d0
     enddo
  else
     do while (RASetting.gt.360.d0)
        RASetting = RASetting - 360.d0
     enddo
  endif

  ! -----------------------------------------------------------------!
  ! NOTE: RA potentially needs to be adjusted into the range [0,360) !
  ! by adding/subtracting 360.                                       !
  ! -----------------------------------------------------------------!

  ! -----------------------------------------------------------------!
  ! Right ascension value needs to be in the same quadrant as L.     !
  ! -----------------------------------------------------------------!

  Lquadrant  = (floor(LRising/90.d0)) * 90.d0
  RAquadrant = (floor(RARising/90.d0)) * 90.d0
  RARising = RARising + (Lquadrant - RAquadrant)    
  
  Lquadrant  = (floor( LSetting/90.d0)) * 90.d0
  RAquadrant = (floor(RASetting/90.d0)) * 90.d0
  RASetting = RASetting + (Lquadrant - RAquadrant)

  ! -----------------------------------------------------------------!
  ! Right ascension value needs to be converted into hours.          !
  ! -----------------------------------------------------------------!

  RARising = RARising / 15.d0
  RASetting = RASetting / 15.d0

  ! -----------------------------------------------------------------!  
  ! Calculate the Sun's declination.                                 !
  ! -----------------------------------------------------------------!

  sinDecRising = 0.39782d0 * sin(pi/180.d0*LRising)
  cosDecRising = cos(asin(sinDecRising))

  sinDecSetting = 0.39782d0 * sin(pi/180.d0*LSetting)
  cosDecSetting = cos(asin(sinDecSetting))
  
  ! -----------------------------------------------------------------!
  ! Calculate the Sun's local hour angle.                            !
  ! -----------------------------------------------------------------!

  cosHRising = (cos(zenith) - (sinDecRising * sin(phiRAD))) / (cosDecRising * &
       cos(phiRAD))
  cosHSetting = (cos(zenith) - (sinDecSetting * sin(phiRAD))) / (cosDecSetting * &
       cos(phiRAD))
  
  ! -----------------------------------------------------------------!
  ! if (cosH >  1)                                                   !
  ! the sun never rises on this location (on the specified date)     !
  ! if (cosH < -1)                                                   !
  ! the sun never sets on this location (on the specified date)      !
  ! -----------------------------------------------------------------!

  ! -----------------------------------------------------------------!
  ! Finish calculating H and convert into hours.                     !
  ! -----------------------------------------------------------------!

  HRising = 360.d0 - 180.d0/pi*acos(cosHRising)
  HSetting = 180.d0/pi*acos(cosHSetting)
  HRising = HRising / 15.d0
  HSetting = HSetting / 15.d0

  ! -----------------------------------------------------------------!
  ! Calculate local mean time of rising/setting.                     !
  ! -----------------------------------------------------------------!

  T_Rising = HRising + RARising - (0.06571d0 * tRising) - 6.622d0
  T_Setting = HSetting + RASetting - (0.06571d0 * tSetting) - 6.622d0

  ! -----------------------------------------------------------------!
  ! Adjust back to UTC.                                              !
  ! -----------------------------------------------------------------!

  UTRising = T_Rising - lngHour
  UTSetting = T_Setting - lngHour

  if (UTRising.le.0.d0) then
     do while (UTRising.lt.0.d0)
        UTRising = UTRising + 24.d0
     enddo
  else
     do while (UTRising.gt.24.d0)
        UTRising = UTRising - 24.d0
     enddo
  endif
  if (UTSetting.le.0.d0) then
     do while (UTSetting.lt.0.d0)
        UTSetting = UTSetting + 24.d0
     enddo
  else     
     do while (UTSetting.gt.24.d0)
        UTSetting = UTSetting - 24.d0
     enddo
  endif

  ! -----------------------------------------------------------------!
  ! NOTE: UT potentially needs to be adjusted into the range [0,24)  !
  ! by adding/subtracting 24.                                        !
  ! -----------------------------------------------------------------!

  ! -----------------------------------------------------------------!
  ! UTRising y UTSetting en fraccion de dia.                         !
  ! -----------------------------------------------------------------!

  UTRising = UTRising/24.d0 
  UTSetting = UTSetting/24.d0

  return
end subroutine day_or_night
