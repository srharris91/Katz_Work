PROGRAM test_io
    !purpose: illustrate I/O of variables with derived data types
    IMPLICIT NONE
    !Declare type person
    TYPE::person
        CHARACTER(len=14)::first_name
        CHARACTER::middle_initial
        CHARACTER(len=14)::last_name
        CHARACTER(len=14)::phone
        INTEGER::age
        CHARACTER::sex
        CHARACTER(len=11)::ssn
    END TYPE person

    !declare variable
    TYPE (person) :: john
    !initialize variable
    john=person('John','R','Jones','323-6439',21,'M','123-45-6789')

    !output variable using free format I/O
    WRITE(*,*) 'Free format: ',john

    !output using formatted I/O
    WRITE(*,1000)john
    1000 FORMAT (' Formatted I/O:' ,/,4(1x,A,/),1x,I4,/,1x,A,/,1x,A)
END PROGRAM test_io
