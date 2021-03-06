/**
 * \brief
 * Class Sutherland holds the data and specifies the operations for
 * transport properties using Sutherland's Law.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-13
 */


#ifndef included_Sutherland
#define included_Sutherland

#include "PHYSICS_defs.h"
#include "TransportType.h"


class Sutherland: public TransportType
{

 public:

  /**
   * \brief
   * Constructor for the Sutherland class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C Sutherland
   */
  Sutherland();

  /**
   * \brief
   * Destructor for the Sutherland class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C ~Sutherland
   */
  ~Sutherland();

  /**
   * \brief
   * Reads inputs and initializes the Sutherland class objects.
   * \param comp Component number of this instance of Sutherland's Law.
   * \param inputFile Name of Transport input file.
   * \param state Pointer to the current State instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C initialize
   */
  void initialize(const int& comp,
		  const string& inputFile,
		  State* state);

  /**
   * \brief
   * Returns the Prandtl number.
   * \param Prn Prandtl number.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-06-9
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getPrn
   */
  void getPrn(double* Prn);

  /**
   * \brief
   * Returns the turbulent Prandtl number.
   * \param PrnT Turbulent Prandtl number.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getPrnT
   */
  void getPrnT(double* PrnT);

  /**
   * \brief
   * Returns dynamic viscosity given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param mu Dynamic viscosity.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getViscosity
   */
  void getViscosity(const int& npts,
		    const double* p,
		    const double* t,
		    double* mu);

  /**
   * \brief
   * Returns the derivative of dynamic viscosity with respect to x
   * given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param px Pressure x-derivative.
   * \param tx Temperature x-derivative.
   * \param mux Dynamic viscosity x-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getViscosityX
   */
  void getViscosityX(const int& npts,
		     const double* p,
		     const double* t,
		     const double* px,
		     const double* tx,
		     double* mux);

  /**
   * \brief
   * Returns the derivative of dynamic viscosity with respect to y
   * given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param py Pressure y-derivative.
   * \param ty Temperature y-derivative.
   * \param muy Dynamic viscosity y-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getViscosityY
   */
  void getViscosityY(const int& npts,
		     const double* p,
		     const double* t,
		     const double* py,
		     const double* ty,
		     double* muy);

  /**
   * \brief
   * Returns the derivative of dynamic viscosity with respect to z
   * given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param pz Pressure z-derivative.
   * \param tz Temperature z-derivative.
   * \param muz Dynamic viscosity z-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getViscosityZ
   */
  void getViscosityZ(const int& npts,
		     const double* p,
		     const double* t,
		     const double* pz,
		     const double* tz,
		     double* muz);

  /**
   * \brief
   * Returns the derivative of dynamic viscosity wrt temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param dmudT Derivative of dynamic viscosity wrt temperature.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getDViscosityDT
   */
  void getDViscosityDT(const int& npts,
		       const double* p,
		       const double* t,
		       double* dmudT);

  /**
   * \brief
   * Returns thermal conductivity given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param k Thermal conductivity.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getConductivity
   */
  void getConductivity(const int& npts,
		       const double* p,
		       const double* t,
		       double* k);

  /**
   * \brief
   * Returns the derivative of thermal conductivity with respect to x
   * given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param px Pressure x-derivative.
   * \param tx Temperature x-derivative.
   * \param kx Thermal conductivity x-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getConductivityX
   */
  void getConductivityX(const int& npts,
			const double* p,
			const double* t,
			const double* px,
			const double* tx,
			double* kx);

  /**
   * \brief
   * Returns the derivative of thermal conductivity with respect to y
   * given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param py Pressure y-derivative.
   * \param ty Temperature y-derivative.
   * \param ky Thermal conductivity y-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getConductivityY
   */
  void getConductivityY(const int& npts,
			const double* p,
			const double* t,
			const double* py,
			const double* ty,
			double* ky);

  /**
   * \brief
   * Returns the derivative of thermal conductivity with respect to z
   * given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param pz Pressure z-derivative.
   * \param tz Temperature z-derivative.
   * \param kz Thermal conductivity z-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getConductivityZ
   */
  void getConductivityZ(const int& npts,
			const double* p,
			const double* t,
			const double* pz,
			const double* tz,
			double* kz);

  /**
   * \brief
   * Returns the derivative of thermal conductivity wrt temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param dkdT Derivative of thermal conductivity wrt temperature.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C getDConductivityDT
   */
  void getDConductivityDT(const int& npts,
			  const double* p,
			  const double* t,
			  double* dkdT);

  /**
   * \brief
   * Deallocates any memory for Sutherland class objects.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Transport/Sutherland.C finalize
   */
  void finalize();


 private:

  double mu0;
  double t0;
  double s;
  double prn;
  double prnT;

};
#endif
