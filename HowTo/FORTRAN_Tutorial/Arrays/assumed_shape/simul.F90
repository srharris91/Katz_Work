SUBROUTINE simul (a,b,ndim,n,error)
    !purpose: solve a set of n linear equations in n unknowns using Gaussian elimination and the max pivot technique
    IMPLICIT NONE
    !data dictionary: parameters
    INTEGER,INTENT(IN) :: ndim
    REAL,INTENT(INOUT),DIMENSION(ndim,ndim):: a  !matrix
    REAL,INTENT(INOUT),DIMENSION(ndim):: b       ! rhs
    INTEGER,INTENT(IN):: n                       !# of equations
    INTEGER,INTENT(OUT):: error                  !0 = no error, 1-singular equations
    !data dictionary: constants
    REAL,PARAMETER::EPSILON=1.0E-6
    !DATA dictionary: variables
    REAL::factor
    INTEGER::irow,ipeak,jrow,kcol
    REAL::temp

    !process n times to get all equations
    mainloop: DO irow=1,n
        !find peak pivot for column irows
        ipeak=irow
        max_pivot: DO jrow=irow+1,n
            IF (ABS(a(jrow,irow)) > ABS(a(ipeak,irow))) THEN
                ipeak=jrow
            END IF
        END DO max_pivot

        ! check for singular equations
        singular: IF (ABS(a(ipeak,irow)) < EPSILON) THEN
            error=1
            RETURN
        END IF singular

        ! otherwise swap equations irow&ipeak if ipeak/=irow
        swap_eqn: IF (ipeak/=irow) THEN
            DO kcol=1,n
                temp    = a(ipeak,kcol)
                a(ipeak,kcol)=a(irow,kcol)
                a(irow,kcol) =temp
            END DO
            temp    = b(ipeak)
            b(ipeak)= b(irow)
            b(irow) = temp
        END IF swap_eqn

        !multiply equation by factor
        eliminate: DO jrow=1,n
            IF (jrow/=irow) THEN
                factor=-a(jrow,irow)/a(irow,irow)
                DO kcol=1,n
                    a(jrow,kcol)=a(irow,kcol)*factor + a(jrow,kcol)
                END DO
                b(jrow)=b(irow)*factor+b(jrow)
            END IF
        END DO eliminate
    END DO mainloop
    !all diagonal terms are now zero, onw must divide each equation by coefficient of its on-diagonal term
    divide: DO irow=1,n
        b(irow)     = b(irow)/a(irow,irow)
        a(irow,irow)= 1.
    END DO divide

    ! set error flag to 0 and return
    error = 0
END SUBROUTINE simul
