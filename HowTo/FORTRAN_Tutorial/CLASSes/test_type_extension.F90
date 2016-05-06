PROGRAM test_type_extension
    !purpose: illustrate use of extension of derived data types
    IMPLICIT NONE

    !declare type point
    TYPE :: point
        REAL::x,y
    END TYPE

    !declare inherited extension type
    TYPE, EXTENDS(point) :: point3d
        REAL :: z
    END TYPE

    !declare a variable and print out
    TYPE (point3d)::my_point
    !initialize
    my_point%x = 1.
    my_point%y = 2.
    my_point%z = 3.

    !print out
    WRITE(*,*)'my_point = ',my_point
END PROGRAM test_type_extension
