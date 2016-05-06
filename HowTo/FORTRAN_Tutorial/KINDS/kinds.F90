PROGRAM kinds
    !purpose: determine the kinds of single and double precision real values in a particular computer
    IMPLICIT NONE
    !Write out the kinds of single and double precision values
    WRITE(*,'(" The KIND for a single precision is",I2)')KIND(0.0)
    WRITE(*,'(" The KIND for a double precision is",I2)')KIND(0.0D0)
END PROGRAM kinds


