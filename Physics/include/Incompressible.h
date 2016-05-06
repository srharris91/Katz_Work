/**
 * \brief
 * Class Incompressible holds the data and specifies the operations for
 * an incompressible fluid.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-13
 */


#ifndef included_Incompressible
#define included_Incompressible

#include "PHYSICS_defs.h"
#include "StateType.h"


class Incompressible: public StateType
{

 public:

  /**
   * \brief
   * Constructor for the Incompressible class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/State/Incompressible.C Incompressible
   */
  Incompressible();

  /**
   * \brief
   * Destructor for the Incompressible class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/State/Incompressible.C ~Incompressible
   */
  ~Incompressible();

  /**
   * \brief
   * Reads inputs and initializes the Incompressible class objects.
   * \param comp Component number of this instance of ideal gas.
   * \param inputFile Name of State input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/State/Incompressible.C initialize
   */
  void initialize(const int& comp,
		  const string& inputFile);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getGamma
   */
  void getGamma(const int& npts,
		const double* p,
		const double* t,
		double* gamma);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getCp
   */
  void getCp(const int& npts,
	     const double* p,
	     const double* t,
	     double* Cp);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getSoundSpeed
   */
  void getSoundSpeed(const int& npts,
		     const double* p,
		     const double* t,
		     double* c);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getRhoP
   */
  void getRhoP(const int& npts,
	       const double* p,
	       const double* t,
	       double* rp);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getRhoT
   */
  void getRhoT(const int& npts,
	       const double* p,
	       const double* t,
	       double* rt);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getHP
   */
  void getHP(const int& npts,
	     const double* p,
	     const double* t,
	     double* hp);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getHT
   */
  void getHT(const int& npts,
	     const double* p,
	     const double* t,
	     double* ht);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getDensity
   */
  void getDensity(const int& npts,
		  const double* p,
		  const double* t,
		  double* rho);

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
  void getDensityX(const int& npts,
		   const double* p,
		   const double* t,
		   const double* px,
		   const double* tx,
		   double* rhox);

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
  void getDensityY(const int& npts,
		   const double* p,
		   const double* t,
		   const double* py,
		   const double* ty,
		   double* rhoy);

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
  void getDensityZ(const int& npts,
		   const double* p,
		   const double* t,
		   const double* pz,
		   const double* tz,
		   double* rhoz);

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
   * \par Source Code:
   * \snippet src/State/Incompressible.C getEnthalpy
   */
  void getEnthalpy(const int& npts,
		   const double* p,
		   const double* t,
		   double* h);

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
  void getEnthalpyX(const int& npts,
		    const double* p,
		    const double* t,
		    const double* px,
		    const double* tx,
		    double* hx);

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
  void getEnthalpyY(const int& npts,
		    const double* p,
		    const double* t,
		    const double* py,
		    const double* ty,
		    double* hy);

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
  void getEnthalpyZ(const int& npts,
		    const double* p,
		    const double* t,
		    const double* pz,
		    const double* tz,
		    double* hz);

  /**
   * \brief
   * Deallocates any memory for Incompressible class objects.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/State/Incompressible.C finalize
   */
  void finalize();


 private:

  double rho;
  double Cp;
  double Cv;
  double c;
};
#endif
