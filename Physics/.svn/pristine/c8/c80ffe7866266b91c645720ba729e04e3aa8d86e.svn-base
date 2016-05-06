/**
 * \brief
 * Class StateType is the base class for various
 * equations of state, such as ideal gas.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-13
 */


#ifndef included_StateType
#define included_StateType

#include "PHYSICS_defs.h"


class StateType
{
 public:

  /**
   * \brief
   * Constructor for the StateType class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StateType/StateType.C StateType
   */
  StateType();

  /**
   * \brief
   * Destructor for the StateType class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StateType/StateType.C ~StateType
   */
  virtual ~StateType();

  /**
   * \brief
   * Reads inputs and allocates/initializes the StateType class objects.
   * \param icomp Fluid component number.
   * \param inputFile Name of StateType input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   */
  virtual void initialize(const int& icomp,
			  const string& inputFile)=0;

  /**
   * \brief
   * Returns Gamma, given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param gamma Ratio of specific heats.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   */
  virtual void getGamma(const int& npts,
			const double* p,
			const double* t,
			double* gamma)=0;

  /**
   * \brief
   * Returns Gamma for constant gamma equations of state.
   * \param gamma Ratio of specific heats.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/State/StateType.C getGamma2
   */
  virtual void getGamma(double* gamma);

  /**
   * \brief
   * Returns the ideal gas constant for ideal gases.
   * \param rGas Ideal gas constant.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/State/StateType.C getRGas
   */
  virtual void getRGas(double* rGas);

  /**
   * \brief
   * Returns Cp, given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param Cp Constant pressure specific heat.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   */
  virtual void getCp(const int& npts,
		     const double* p,
		     const double* t,
		     double* Cp)=0;

  /**
   * \brief
   * Returns sound speed given pressure and temperature for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param c Sound Speed.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-12
   * \par Further Documentation:
   */
  virtual void getSoundSpeed(const int& npts,
			     const double* p,
			     const double* t,
			     double* c)=0;

  /**
   * \brief
   * Returns the derivative of density with respect to pressure at constant
   * temperature for separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param rp d rho/ d p.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   */
  virtual void getRhoP(const int& npts,
		       const double* p,
		       const double* t,
		       double* rp)=0;

  /**
   * \brief
   * Returns the derivative of density with respect to temperature at constant
   * pressure for separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param rt d rho/ d t.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   */
  virtual void getRhoT(const int& npts,
		       const double* p,
		       const double* t,
		       double* rt)=0;

  /**
   * \brief
   * Returns the derivative of enthalpy with respect to pressure at constant
   * temperature for separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param hp d h/ d p.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   */
  virtual void getHP(const int& npts,
		     const double* p,
		     const double* t,
		     double* hp)=0;

  /**
   * \brief
   * Returns the derivative of enthalpy with respect to temperature at constant
   * pressure for separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param ht d h/ d t.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   */
  virtual void getHT(const int& npts,
		     const double* p,
		     const double* t,
		     double* ht)=0;

  /**
   * \brief
   * Returns density given pressure and temperature.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param rho Density.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   */
  virtual void getDensity(const int& npts,
			  const double* p,
			  const double* t,
			  double* rho)=0;

  /**
   * \brief
   * Returns the derivative of density with respect to x,
   * given pressure and temperature and their derivatives for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param px Pressure x-derivative.
   * \param tx Temperature x-derivative.
   * \param rhox Density x-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-03
   * \par Further Documentation:
   */
  virtual void getDensityX(const int& npts,
			   const double* p,
			   const double* t,
			   const double* px,
			   const double* tx,
			   double* rhox)=0;

  /**
   * \brief
   * Returns the derivative of density with respect to y,
   * given pressure and temperature and their derivatives for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param py Pressure y-derivative.
   * \param ty Temperature y-derivative.
   * \param rhoy Density y-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-03
   * \par Further Documentation:
   */
  virtual void getDensityY(const int& npts,
			   const double* p,
			   const double* t,
			   const double* py,
			   const double* ty,
			   double* rhoy)=0;

  /**
   * \brief
   * Returns the derivative of density with respect to z,
   * given pressure and temperature and their derivatives for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param pz Pressure z-derivative.
   * \param tz Temperature z-derivative.
   * \param rhoz Density z-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-03
   * \par Further Documentation:
   */
  virtual void getDensityZ(const int& npts,
			   const double* p,
			   const double* t,
			   const double* pz,
			   const double* tz,
			   double* rhoz)=0;

  /**
   * \brief
   * Returns enthalpy given pressure and temperature.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param h Enthalpy.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   */
  virtual void getEnthalpy(const int& npts,
			  const double* p,
			  const double* t,
			  double* h)=0;

  /**
   * \brief
   * Returns the derivative of enthalpy with respect to x,
   * given pressure and temperature and their derivatives for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param px Pressure x-derivative.
   * \param tx Temperature x-derivative.
   * \param hx Enthalpy x-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-03
   * \par Further Documentation:
   */
  virtual void getEnthalpyX(const int& npts,
			    const double* p,
			    const double* t,
			    const double* px,
			    const double* tx,
			    double* hx)=0;

  /**
   * \brief
   * Returns the derivative of enthalpy with respect to y,
   * given pressure and temperature and their derivatives for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param py Pressure y-derivative.
   * \param ty Temperature y-derivative.
   * \param hy Enthalpy y-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-03
   * \par Further Documentation:
   */
  virtual void getEnthalpyY(const int& npts,
			    const double* p,
			    const double* t,
			    const double* py,
			    const double* ty,
			    double* hy)=0;

  /**
   * \brief
   * Returns the derivative of enthalpy with respect to z,
   * given pressure and temperature and their derivatives for
   * separate, non-mixture components.
   * \param npts Number of points to operate on.
   * \param p Pressure.
   * \param t Temperature.
   * \param pz Pressure z-derivative.
   * \param tz Temperature z-derivative.
   * \param hz Enthalpy z-derivative.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-03
   * \par Further Documentation:
   */
  virtual void getEnthalpyZ(const int& npts,
			    const double* p,
			    const double* t,
			    const double* pz,
			    const double* tz,
			    double* hz)=0;

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

 private:

};
#endif
