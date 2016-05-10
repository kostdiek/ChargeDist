#ifndef INCtoih
#define INCtoih

#include <string>

/** \file toi.h
 * \brief File in which the Isotpe class is declared.
 */

using namespace std;

class Isotope;

//! This class describe an isotope from <A HREF="http://www.nndc.bnl.gov/masses/"> The Ame2003 atomic mass evaluation </A>
/*!
  by G.Audi, A.H.Wapstra and C.Thibault Nuclear Physics A729 p. 337-676, December 22, 2003.
  It basicly read the AME2003 table and fill the class variable.
*/
class Isotope
{ private :
    //! 1 if contain all the require information else 0.
    int isOK ;
    //! mass excess.
    double mex ;
    //! Number of nucleon.
    int A;
    //! Number of proton.
    int Z;
    //! Number of neutron.
    int N;
    //! Mass in MeV/c2.
    double mass;
    //! Symbol.
    string   symb;
    //! Find the isotpe fill the variable.
    void FindIsotope(int _N, int _Z, int _A, string _Symb);

  public :
    //! Default constructor.
    Isotope(){isOK=0;}
    //! Constructor using number of Neutron (N) and atomic number (Z).
    Isotope(int _N, int _Z);
    //! Constructor using Number of nucleon (A) and Symbol.
    Isotope(int _A, string _Symb);
    //! Constructor using full name.
    Isotope(string _Name);
    //! Copy constructor.
    Isotope(Isotope& i);
    //! Set isotope using full name.
    int SetIsotope(string _Name);
    //! Set isotope using number of Neutron (N) and atomic number (Z).
    int SetIsotope(int _N, int _Z);
    //! Return the mass excess.
    double GetMassExcess() { return isOK==0 ? -1: mex ; }
    //! return the mass in Mev/c2.
    double GetMass() { return isOK==0 ? -1:mass ; }
    //! return the full name (e.g. 12C).
    string GetName();
    //! return the isotope symbol (e.g. C).
    string GetSymbol() {if(isOK==0)return "-1";else return symb;}
    //! return the number of nucleon.
    int GetA() { return isOK==0 ? -1:A ; }
    //! return the atomic number.
    int GetZ() { return isOK==0 ? -1:Z ; }
    //! return the number of neutron.
    int GetN() { return isOK==0 ? -1:N ; }
    //! return 1 if the isotope is fully described 0 otherwise.
    int GetCheck() {return isOK;}
} ;
#endif
