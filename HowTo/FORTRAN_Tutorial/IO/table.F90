PROGRAM table
!purpose: illustrate formatted write statements.  This programs gives a table containing square roots, squares, and cubes of all integers between 1 and 10.  The table includes a title and column headings.
IMPLICIT NONE

INTEGER::cube,i,square
REAL::square_root

!Print the title of the table on a new page
WRITE(*,100)
100 FORMAT ('1',T3,'Table of Squazre Roots, Squares, and Cubes')

!Print the column headings after skipping one line.
WRITE(*,110)
110 FORMAT ('0',T4,'Number',T13,'Square Root',T29,'Square',T39,'Cubes')
WRITE(*,120)
120 FORMAT (1x, T4,'======',T13,'===========',T29,'======',T39,'====='///) !slash / does the output to buffer.  This is similar to <<endl<<endl<<endl; in c++
! 1x means 1 space, T4 means tab to 4th column

!Generate the required values, and print them out
DO i=1,10
    square_root=sqrt(real(i))
    square=i**2
    cube=i**3
    WRITE(*,130)i,square_root,square,cube
    130 FORMAT (1x, T4, I4, T13, F10.6, T27, I6, T37, I6)
! see this for more formatting stuff http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
! it will explain the F10.6 and I4 and stuff
END DO

WRITE(*,*) 'testing stuff and spaces'
WRITE(*,*) 'next line'//
WRITE(*,*) 'double space above?'

END PROGRAM table
