#!/bin/bash

rm BJDs JDs
gfortran -o main main.f90 BJDtoJD.f90 precession.f90 day_or_night.f90 star_altitude.f90 jd2date_2.f90 ut2hhmm.f90 RA_DEC_sgs.f90

./main
