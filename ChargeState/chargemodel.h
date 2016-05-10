#ifndef INCChargeModel
#define INCChargeModel

#include <iostream>
#include <vector>

#include "toi.h"

/** \file chargemodel.h
 * \brief File in which the charge distribution class are declared.
 */

//! This the mother of ALL charge State model.
/*!
 All the model
calcuating charge state distribution must derive from this
function.

 Three virtual function are declared, they must be impelemented
in the class describing the model.
*/
class ChargeStateDist
{
 private :
   //! Define the ion for which you want to calculate the Charge state fraction.
   Isotope ion;
   //! Check if the ion is defined.
   int IonSet;
   //! Vector containing all the critical value for the calculation - Will be filled the first time fq is called
   std::vector<float> criticaldata;
   //! Vector containing the fq for all charge state
   std::vector<double> fqvector;
   //! Return the unormalized fraction of the charge z at energy E. To be implemented in the class describing the model.
   virtual double GetUnNormalizefq(double z,double E)=0;
   //! Member which return the critical value of the model. To be implemented in the class describing the model.
   virtual const std::vector<float> GetCriticalData()=0;
   //!Last energy used to compute fq vector
   double energy;
   //!Member which fill the criticaldata vector with the data from GetCriticalData()
   bool FillCritical();
   //!Member which compare criticaldata with the current status of the model
   bool CompareCritical();
 public :
   //! Constructor
   ChargeStateDist(){IonSet=0; criticaldata.reserve(4); energy=-1.;fqvector.reserve(10);}
   //! Constructor (this is the most obvious choice).
   ChargeStateDist(Isotope a){ion=a; IonSet=1; criticaldata.reserve(4);energy=-1.;fqvector.reserve(a.GetZ());}
   //!vitual destructor not implemented. If you use pointed you could be in trouble.
   virtual ~ChargeStateDist(){};
   //! Check if the ion is defined.
   int GetCheck(){return IonSet;}
   //! Return the ion
   Isotope GetIon(){return ion;}
   //! Set the isotope == ion.
   void SetIsotope(Isotope a){ion=a; IonSet=1;fqvector.reserve(a.GetZ());}
   //! Set the isotope == ion.
   void SetIon(Isotope a){ion=a; IonSet=1;fqvector.reserve(a.GetZ());}
   //! Member which return the charge state fraction. To be implemented in the class describing the model.
   virtual double GetqBar(double E)=0;
   //! Member returning the fraction of charge z at an energy E.
   double Getfq(double z,double E);
   //! Member returning the most abundant charge state at the energy E
   double GetMaxq(double E);
   //! Print the reference of the charge state model. To be implemented in the class describing the model.
   virtual void WriteReference()=0;
   //! Member returning the vector with the critical data
   const std::vector<float> GetVector() { return criticaldata; }
} ;

//! This is the <A HREF="../pdf/Liu_charge_state_studies_of_low_energy_ions_passing_through_hydrogen_and_helium_gas.pdf">Liu model (pdf)</A>.
/*!
  The documentation of the class is provided as an example. The choice to document this class is due
to the fact that this model is medium (target) dependent which make it slightly more complex.
*/
class Charge_Liu : public ChargeStateDist
{
 private :
   //! Gas name H or He
   std::string Gas;
   //! Check if the interacting medium (target) or gas is selected
   int GasSet;
   //! Parameter dependant or the gas (target)
   double par[3];
   //! Return the not normalized fraction of the z at energy E. Actual implementation of the GetUnNormalizefq function.
   double GetUnNormalizefq(double z,double E);
   //! Fill the vector critical Data with N,Z and a flag for the gas (1 for H 2 for He)
   const std::vector<float> GetCriticalData();
 public :
   //! Constructor with ion choice but no gas. If j=0 the reference are not displayed.
   Charge_Liu(Isotope i,int j=1):ChargeStateDist(i){if(j)WriteReference();GasSet=0;}
   //! Constructor. If j=0 the reference are not displayed.
   Charge_Liu(Isotope i, std::string s,int j=1);
   //! Set the gas choice.
   int SetGas(std::string s);
   //! Return 1 if the gas is selected
   int GetGasCheck(){return GasSet;}
   //! Return the gas name (H or He)
   std::string GetGas(){return Gas;}
   //! Return qbar. Actual implementation.
   double GetqBar(double E);
   //! Display the reference. Actual implementation.
   void WriteReference();
};

//! This is the <A HREF="../pdf/Sayer_Rev_Phys_App12_1977.pdf">Sayer model (pdf)</A>.
/*!
  See the Liu model for class documentations.
*/
class Charge_Sayer : public ChargeStateDist
{
 private :
   std::string Medium;
   int MediumSet;
   double par[10];
   double GetUnNormalizefq(double z,double E);
   //! Fill the vector critical Data with N,Z and a flag for solid (1) or gas (2)
   const std::vector<float> GetCriticalData();
 public :
   Charge_Sayer(Isotope i, int j=1):ChargeStateDist(i){if(j)WriteReference();MediumSet=0;}
   Charge_Sayer(Isotope i, std::string s,int j=1);
   int SetMedium(std::string s);
   int GetMediumCheck(){return MediumSet;}
   std::string GetMedium(){return Medium;}
   double GetqBar(double E);
   void WriteReference();
};

//! This is the <A HREF="../pdf/Nikolaev_Dmitriev_Phys_lett28A_1968.pdf">Nikolaev and Dmitriev 1968 model (pdf)</A>.
/*!
  See the Liu model for class documentations.
*/
class Charge_Zdist : public ChargeStateDist
{
 private :
   double GetUnNormalizefq(double z,double E);
   //! Fill the vector critical Data with N,Z
   const std::vector<float> GetCriticalData();
 public :
   Charge_Zdist(Isotope i,int j=1):ChargeStateDist(i){if(j)WriteReference();}
   double GetqBar(double E);
   void WriteReference();
};

//! This is the <A HREF="../pdf/Dmitriev_Nikolaev_JETP20_1965.pdf">Dmitriev and Nikolaev 1965 model (pdf)</A>.
/*!
  See the Liu model for class documentations.
*/
class Charge_Dmitriev : public ChargeStateDist
{
 private :
   double GetUnNormalizefq(double z,double E);
   //! Fill the vector critical Data with N,Z
   const std::vector<float> GetCriticalData();
 public :
   Charge_Dmitriev(Isotope i,int j=1):ChargeStateDist(i){if(j)WriteReference();}
   double GetqBar(double E);
   void WriteReference();
};

//! This is the <A HREF="../pdf/Baudinet-Robinet_NIM190_1981.pdf">Baudinet-Robinet 1981 model (pdf)</A>.
/*!
  See the Liu model for class documentations.
*/
class Charge_Baudinet : public ChargeStateDist
{
 private :
   double gamma(double nu);
   double GetUnNormalizefq(double z,double E);
   //! Fill the vector critical Data with N,Z
   const std::vector<float> GetCriticalData();
 public :
   Charge_Baudinet(Isotope i,int j=1):ChargeStateDist(i){if(j)WriteReference();}
   double GetqBar(double E);
   void WriteReference();
};

#endif
