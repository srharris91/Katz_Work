MODULE vector_module_2
    !purpose: derived data type for 2d vectors, plus addition and subtraction operations

    IMPLICIT NONE
    !declare data type
    TYPE :: vector
        REAL :: x,y
        CONTAINS
            PROCEDURE,PASS :: vector_add
            PROCEDURE,PASS :: vector_sub
    END TYPE vector

    !add procedures
    CONTAINS
        TYPE (vector) FUNCTION vector_add (this,v2)
            !purpose: add two vectors
            IMPLICIT NONE
            !data dictionary: parameters
            CLASS (vector),INTENT(IN)::this,v2 !this is the vector class it is addressed with
            
            !add the vectors
            vector_add%x = this%x + v2%x
            vector_add%y = this%y + v2%y
        END FUNCTION vector_add

        TYPE (vector) FUNCTION vector_sub(this,v2)
            !purpose: subtract two vectors
            IMPLICIT NONE
            !data dictionary: parameters
            CLASS (vector),INTENT(IN)::this,v2 !this is the vector class it is addressed with
            
            !subtract the points
            vector_sub%x = this%x - v2%x
            vector_sub%y = this%y - v2%y
        END FUNCTION vector_sub
            
END MODULE vector_module_2
