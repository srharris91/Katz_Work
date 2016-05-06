PROGRAM test_vectors
    !Purpose: test the adding and subtracting of 2d vectors
    USE vector_module
    IMPLICIT NONE
    !enter points
    TYPE (vector) :: v1,v2
    
    !get first and second vector
    WRITE(*,*) 'Enter the first sector (x,y):'
    READ(*,*) v1%x,v1%y
    WRITE(*,*) 'Enter the second vector (x,y):'
    READ(*,*) v2%x, v2%y
    ! add the points
    WRITE(*,1000) vector_add(v1,v2)
    1000 FORMAT (1x,'The sum o fthe points is (',F8.2,',',F8.2')')

    !subtract the points
    WRITE(*,1010) vector_sub(v1,v2)
    1010 FORMAT (1x,'The difference of the poiunts is (',F8.2,',',F8.2,')')
END PROGRAM test_vectors
