MODULE set3d_mod
    !Shauns c hange for fun
    IMPLICIT NONE
    TYPE::SET
        !REAL::phi,phiS,phiN,gradPhiMag,gradPhi,grad2Phi,gradMixPhi,gridX
        REAL:: x,y,z                        ! x,y,z coordinates for each node
        REAL:: phi,phiN                     ! phi old and new values.
        REAL:: gradPhiX,gradPhiY, gradPhiZ  ! gradient of phi in x,y, and z direction
        REAL:: gradPhiXX                ! second gradient of phi in x and then x direction
        REAL:: gradPhiYY                ! second gradient of phi in y and then y direction
        REAL:: gradPhiZZ                ! second gradient of phi in z and then z direction
        REAL:: gradPhiXY                ! second gradient of phi in x and then y direction
        REAL:: gradPhiXZ                ! second gradient of phi in x and then z direction
        REAL:: gradPhiYZ                ! second gradient of phi in y and then z direction
        REAL:: gradPhiMag               ! Magnitude of the gradient of phi
        REAL:: gausCurv                 ! gausCurve 
        REAL:: curv                     ! curvature
        REAL:: gradPhiDir               ! gradient direction (positive or negative)
        REAL:: principleCurvMax         ! Max and Min principle curvature 
        REAL:: principleCurvMin
    END TYPE SET
    
    TYPE::SET_smaller
        REAL::phi,gradPhiX,gradPhiY,gradPhiZ
!        REAL::phi1,phi2
    END TYPE SET_smaller
END MODULE set3d_mod
