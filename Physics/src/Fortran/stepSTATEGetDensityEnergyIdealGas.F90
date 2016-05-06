!> \brief
!!   This subroutine returns density and internal energy given pressure
!!   and temperature for and ideal gas.
!! \param npts Number of points to operate on.
!! \param rGas Ideas gas constant.
!! \param gamma Ratio of specific heats.
!! \param p Pressure.
!! \param t Temperature.
!! \param rho Density.
!! \param e Internal Energy.
!! \author:
!!   Aaron Katz
!! \version:
!!   1.0
!! \date
!!   2012-9-13
!! \par Further Documentation:
!! \par Source code:
!! \include src/Fortran/stepSTATEGetDensityEnergyIdealGas.F90


SUBROUTINE stepstategetdensityenergyidealgas(&
     npts, &
     rGas, &
     gamma, &
     p, &
     t, &
     rho, &
     e)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: npts
REAL,   INTENT(IN   ) :: rGas,gamma
REAL,   INTENT(IN   ),DIMENSION(npts) :: p,t
REAL,   INTENT(  OUT),DIMENSION(npts) :: rho,e
INTEGER :: n
REAL :: gm1


gm1 = gamma-1.

DO n=1,npts
   rho(n) = p(n)/(rgas*t(n))
   e(n)   = rgas*t(n)/gm1
END DO


END SUBROUTINE stepstategetdensityenergyidealgas
