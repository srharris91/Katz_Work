MODULE types
    !purpose: define the derived data type used for hte customer database
    IMPLICIT NONE
    !declare personal info
    TYPE::personal_info
        CHARACTER(len=12)   ::  first   !first name
        CHARACTER           ::  mi      !middle initial
        CHARACTER(len=12)   ::  last
        CHARACTER(len=26)   ::  street  !street address
        CHARACTER(len=12)   ::  city
        CHARACTER(len=2)    ::  state
        INTEGER             ::  zip
    END TYPE personal_info

END MODULE types

PROGRAM customer_database
    !purpose: read in a character input data set, sort it and use selection algorithm.  
    USE types
    IMPLICIT NONE

    !data dictionary: constants
    INTEGER,PARAMETER::MAX_SIZE=100
    LOGICAL,EXTERNAL::lt_last
    LOGICAL,EXTERNAL::lt_city
    LOGICAL,EXTERNAL::lt_zip

    !Data dictionary: variables
    TYPE(personal_info),DIMENSION(MAX_SIZE)::customers

    INTEGER::choice
    LOGICAL::exceed=.FALSE.
    CHARACTER(len=20)::filename
    INTEGER::i,nvals=0,status
    TYPE(personal_info)::temp

    ! get file containing the input data
    WRITE(*,*) 'Enter the file name with customer database: '
    READ(*,'(A20)') filename

    !open input file
    OPEN (UNIT=9,FILE=filename,STATUS='OLD',IOSTAT=status)
    !successful open then
    fileopen: IF (status==0) THEN
        DO
            READ(9,1010,IOSTAT=status) temp
            1010 FORMAT (A12,1x,A1,1x,A12,1x,A26,1x,A12,1x,A2,1x,I5)
            IF (status/=0) EXIT
            nvals=nvals+1
            size: IF (nvals<=MAX_SIZE) THEN
                customers(nvals)=temp
            ELSE
                exceed=.TRUE.
            END IF size
        END DO

        !was the array size exceeded... 
        toobig: IF (exceed) THEN
            WRITE(*,1020) nvals, MAX_SIZE
            1020 FORMAT (' Maximum array size exceeded: ',I6,' > ',I6)
        ELSE 
            WRITE(*,1030)
            1030 FORMAT (1x,'Enter way to sort database:',/, &
                1x,'   1 -- By last name ',/, &
                1x,'   2 -- By city ',/, &
                1x,'   3 -- By zip code ')
            READ(*,*) choice
            !sort database
            SELECT CASE (choice)
            CASE(1)
                CALL sort_database(customers,nvals,lt_last)
            CASE(2)
                CALL sort_database(customers,nvals,lt_city)
            CASE(3)
                CALL sort_database(customers,nvals,lt_zip)
            CASE DEFAULT
                WRITE(*,*) 'Invalid choice entered!'
            END SELECT

            ! now write out sorted data
            WRITE(*,'(A)') ' The sorted database values are: '
            WRITE(*,1040) ( customers(i),i=1,nvals)
            1040 FORMAT (1x,A12,1x,A1,1x,A12,1x,A26,1x,A12,1x,A2,1x,I5)
        END IF toobig
    ELSE fileopen
    WRITE(*,'(A,I6)') ' File open error: IOSTAT = ',status
END IF fileopen

END PROGRAM customer_database

SUBROUTINE sort_database(array,n,lt_fun)
    !purpose: sort array into ascending order using selection sort where array is a derived data type "personal_info"
    USE types
    IMPLICIT NONE

    !data dictionary parameters
    INTEGER,INTENT(IN)::n
    TYPE(personal_info),DIMENSION(n),INTENT(INOUT)::array
    LOGICAL,EXTERNAL::lt_fun

    !data dictionary: variables
    INTEGER::i,iptr,j
    TYPE(personal_info)::temp

    !sor tarray
    outer: DO i=1,n-1
        !find min value
        iptr=i
        inner: DO j=i+1,n
            minval: IF (lt_fun(array(j),array(iptr)))THEN
                iptr=j
            END IF minval
        END DO inner

        ! iptr points to min value so swap array(iptr with array(i)
        swap: IF (i/=iptr) THEN
            temp    =   array(i)
            array(i)=   array(iptr)
            array(iptr)=temp
        END IF swap

    END DO outer
END SUBROUTINE sort_database

LOGICAL FUNCTION lt_last(a,b)
    !purpose: compare a and b and find which one is smaller (alphabetical)
    USE types
    IMPLICIT NONE
    !parameters
    TYPE (personal_info),INTENT(IN)::a,b

    !make comparison
    lt_last = LLT(a%last,b%last)
END FUNCTION lt_last

LOGICAL FUNCTION lt_city(a,b)
    !purpose: compare a and b and find which one is smaller (alphabetical)
    USE types
    IMPLICIT NONE
    !parameters
    TYPE (personal_info),INTENT(IN)::a,b

    !make comparison
    lt_city = LLT(a%city,b%city)
END FUNCTION lt_city

LOGICAL FUNCTION lt_zip(a,b)
    !purpose: compare a and b and find which one is smaller (numeric)
    USE types
    IMPLICIT NONE
    !parameters
    TYPE (personal_info),INTENT(IN)::a,b

    !make comparison
    lt_zip = a%zip < b%zip
END FUNCTION lt_zip

