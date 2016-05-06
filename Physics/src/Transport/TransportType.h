/**
 * \brief
 * Class TransportType is the base class for various
 * equations of transport, such as Sutherland's Law.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-13
 */


#ifndef included_TransportType
#define included_TransportType

#include "PHYSICS_defs.h"
#include "State.h"


class TransportType
{
 public:

  /**
   * \brief
   * Constructor for the TransportType class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/TransportType/TransportType.C TransportType
   */
  TransportType();

  /**
   * \brief
   * Destructor for the TransportType class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/TransportType/TransportType.C ~TransportType
   */
  virtual ~TransportType();

  /**
   * \brief
   * Reads inputs and allocates/initializes the TransportType class objects.
   * \param icomp Fluid component number.
   * \param inputFile Name of TransportType input file.
   * \param state Pointer to the current State instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   */
  virtual void initialize(const int& icomp,
			  const string& inputFile,
			  State* state)=0;

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
   * \snippet src/Transport/TransportType.C getPrn
   */
  virtual void getPrn(double* Prn);

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
   * \snippet src/Transport/TransportType.C getPrnT
   */
  virtual void getPrnT(double* PrnT);

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
   */
  virtual void getViscosity(const int& npts,
			    const double* p,
			    const double* t,
			    double* mu)=0;

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
   */
  virtual void getViscosityX(const int& npts,
			     const double* p,
			     const double* t,
			     const double* px,
			     const double* tx,
			     double* mux)=0;

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
   */
  virtual void getViscosityY(const int& npts,
			     const double* p,
			     const double* t,
			     const double* py,
			     const double* ty,
			     double* muy)=0;

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
   */
  virtual void getViscosityZ(const int& npts,
			     const double* p,
			     const double* t,
			     const double* pz,
			     const double* tz,
			     double* muz)=0;

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
   */
  virtual void getDViscosityDT(const int& npts,
			       const double* p,
			       const double* t,
			       double* dmudT)=0;

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
   */
  virtual void getConductivity(const int& npts,
			       const double* p,
			       const double* t,
			       double* k)=0;

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
   */
  virtual void getConductivityX(const int& npts,
				const double* p,
				const double* t,
				const double* px,
				const double* tx,
				double* kx)=0;

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
   */
  virtual void getConductivityY(const int& npts,
				const double* p,
				const double* t,
				const double* py,
				const double* ty,
				double* ky)=0;

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
   */
  virtual void getConductivityZ(const int& npts,
				const double* p,
				const double* t,
				const double* pz,
				const double* tz,
				double* kz)=0;

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
   */
  virtual void getDConductivityDT(const int& npts,
				  const double* p,
				  const double* t,
				  double* dkdT)=0;

  /**
   * \brief
   * Deletes memory.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   */
  virtual void finalize()=0;


 protected:

  State* state;


 private:

};
#endif
