subroutine star_altitude(phiDEG, lambDEG, H, jd_1, jd_0, jd_4,&
     StarH_T1, StarH_T0, StarH_T4, StarRA, StarDEC)

  implicit none

  ! Subroutine variables:
  double precision :: phiDEG
  double precision :: lambDEG
  double precision :: StarDEC
  double precision :: StarRA
  double precision :: jd_1
  double precision :: jd_0
  double precision :: jd_4
  double precision :: StarH_T1
  double precision :: StarH_T0
  double precision :: StarH_T4
  double precision :: Toi
  double precision :: H
  double precision :: H_aux

  ! Internal variables:
  double precision :: aux1
  double precision :: aux2
  double precision :: T
  double precision :: GMST
  double precision :: HA
  double precision :: phi
  double precision :: lamb
  double precision, parameter :: JDo = 2451545.0d0
  double precision, parameter :: anno100 = 36525.d0
  double precision, parameter :: c1 = 280.4606183d0
  double precision, parameter :: c2 = 360.98564736629d0
  double precision, parameter :: c3 = 0.000387933
  double precision, parameter :: c4 = 38710000.0d0
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: R = 6371.d0
  double precision :: alpha

  ! ---------------------------------------------------------- !
  ! This subroutine calculates the altitude of the star, given !
  ! the JULIAN DATE of the mid-transit time and first and      !
  ! fourth contacts.                                           !
  ! HA = hour angle                                            !
  ! GMST = Greenwich mean sidereal time                        !
  ! GMST is needed to go from JD to hour angle.                !
  ! ---------------------------------------------------------- !

  ! ---------------------------------------------------------- !
  ! Convert geographic coordinates to radians.                 !
  ! ---------------------------------------------------------- !

  phi = phiDEG*pi/180.d0
  lamb = lambDEG*pi/180.d0

  ! ---------------------------------------------------------- !
  ! Mid-transit time.                                          !
  ! ---------------------------------------------------------- !

  Toi = jd_0

  aux1 = StarRA
  aux2 = StarDEC

  aux1 = aux1*15.d0*pi/180.d0
  aux2 = aux2*pi/180.d0

  T = (Toi - JDo)/anno100  
  GMST = c1 + c2*(Toi - JDo) + c3*T**2 - T**3/c4 

  if (GMST.ge.0.d0) then
     do while (GMST.gt.360.d0)
        GMST = GMST - 360.d0
     enddo
  elseif (GMST.lt.0.d0) then
     do while (GMST.lt.0.d0)
        GMST = GMST + 360.d0
     enddo
  endif
  
  GMST = GMST*pi/180.d0
  HA = GMST + lamb - aux1
  
  aux2 = abs(aux2)
  phi = abs(phi)

  StarH_T0 = asin(sin(phi)*sin(aux2) + cos(phi)*cos(aux2)*cos(HA))
  StarH_T0 = StarH_T0*180.d0/pi

  ! ---------------------------------------------------------- !
  ! First contact time.                                        !
  ! ---------------------------------------------------------- !

  Toi = jd_1

  aux1 = StarRA
  aux2 = StarDEC

  aux1 = aux1*15.d0*pi/180.d0
  aux2 = aux2*pi/180.d0

  T = (Toi - JDo)/anno100  
  GMST = c1 + c2*(Toi - JDo) + c3*T**2 - T**3/c4 

  if (GMST.ge.0.d0) then
     do while (GMST.gt.360.d0)
        GMST = GMST - 360.d0
     enddo
  elseif (GMST.lt.0.d0) then
     do while (GMST.lt.0.d0)
        GMST = GMST + 360.d0
     enddo
  endif
  
  GMST = GMST*pi/180.d0
  HA = GMST + lamb - aux1

  aux2 = abs(aux2)
  phi = abs(phi)

  StarH_T1 = asin(sin(phi)*sin(aux2) + cos(phi)*cos(aux2)*cos(HA))
  StarH_T1 = StarH_T1*180.d0/pi

  ! ---------------------------------------------------------- !
  ! Fourth contact time.                                       !
  ! ---------------------------------------------------------- !

  Toi = jd_4

  aux1 = StarRA
  aux2 = StarDEC

  aux1 = aux1*15.d0*pi/180.d0
  aux2 = aux2*pi/180.d0

  T = (Toi - JDo)/anno100  
  GMST = c1 + c2*(Toi - JDo) + c3*T**2 - T**3/c4 

  if (GMST.ge.0.d0) then
     do while (GMST.gt.360.d0)
        GMST = GMST - 360.d0
     enddo
  elseif (GMST.lt.0.d0) then
     do while (GMST.lt.0.d0)
        GMST = GMST + 360.d0
     enddo
  endif
  
  GMST = GMST*pi/180.d0
  HA = GMST + lamb - aux1

  aux2 = abs(aux2)
  phi = abs(phi)

  StarH_T4 = asin(sin(phi)*sin(aux2) + cos(phi)*cos(aux2)*cos(HA))
  StarH_T4 = StarH_T4*180.d0/pi

  H_aux = H/1000.d0 ! Elevation in km.
  
  alpha = asin( (2.d0*H_aux*R + H_aux**2)/&
       ((R + H_aux)*sqrt(H_aux**2 + 2.d0*H_aux*R)) )
  alpha = alpha*180.d0/pi

  StarH_T0 = StarH_T0 + alpha
  StarH_T1 = StarH_T1 + alpha
  StarH_T4 = StarH_T4 + alpha

  return
end subroutine star_altitude

