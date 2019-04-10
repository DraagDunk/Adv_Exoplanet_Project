  subroutine ut2hhmm(UT_T0, UT_T1, UT_T4, UT0, UT1, UT4)

    implicit none

    double precision :: UT_T0
    double precision :: UT_T1
    double precision :: UT_T4
    character(5) :: UT0
    character(5) :: UT1
    character(5) :: UT4

    character(2) :: char1_aux
    character(2) :: char2_aux

    
    write(char1_aux,"(i2.2)") int(UT_T0*24.d0)
    write(char2_aux,"(i2.2)") int((UT_T0*24.d0 - int(UT_T0*24.d0))*60.d0)
    UT0 = char1_aux//":"//char2_aux

    write(char1_aux,"(i2.2)") int(UT_T1*24.d0)
    write(char2_aux,"(i2.2)") int((UT_T1*24.d0 - int(UT_T1*24.d0))*60.d0)
    UT1 = char1_aux//":"//char2_aux

    write(char1_aux,"(i2.2)") int(UT_T4*24.d0)
    write(char2_aux,"(i2.2)") int((UT_T4*24.d0 - int(UT_T4*24.d0))*60.d0)
    UT4 = char1_aux//":"//char2_aux

    return
  end subroutine ut2hhmm

