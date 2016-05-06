MODULE vector_module
    !purpose: derived data type for 2d vectors, plus addition and subtraction operations

    IMPLICIT NONE
    !declare data type
    TYPE :: vector
        REAL :: x,y
    END TYPE vector

    !add procedures
    CONTAINS
        TYPE (vector) FUNCTION vector_add (v1,v2)
            !purpose: add two vectors
            IMPLICIT NONE
            !data dictionary: parameters
            TYPE (vector),INTENT(IN)::v1,v2
            
            !add the vectors
            vector_add%x = v1%x + v2%x
            vector_add%y = v1%y + v2%y
        END FUNCTION vector_add

        TYPE (vector) FUNCTION vector_sub(v1,v2)
            !purpose: subtract two vectors
            IMPLICIT NONE
            !data dictionary: parameters
            TYPE (vector),INTENT(IN)::v1,v2
            
            !subtract the points
            vector_sub%x = v1%x - v2%x
            vector_sub%y = v1%y - v2%y
        END FUNCTION vector_sub
            
END MODULE vector_module
