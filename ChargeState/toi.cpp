#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "toi.h"
#include "constant.h"

/** \file toi.cpp
 * \brief File in which the Isotpe class is implemented.
 */

using namespace std;

/*! Nice macro which behave like printf but with it include the std type. example
 * tmp=FSTR(A<<symb) tmp and symb are strings, A is a double tmp finally is the concatenate of
 * A and symb.
*/
#define FSTR( value ) (((ostringstream&)(ostringstream().flush() << value)).str())

Isotope::Isotope(int _N, int _Z)
{
 SetIsotope(_N,_Z);
}

Isotope::Isotope(int _A, string _Symb)
{
 isOK=0; 
 if (!(_Symb.size()>2)) FindIsotope(0,0,_A,_Symb);
}

Isotope::Isotope(string _Name)
{
 SetIsotope(_Name);
}

Isotope::Isotope(Isotope& i)
{
 isOK=i.GetCheck();
 mex=i.GetMassExcess();
 A=i.GetA(); Z=i.GetZ(); N=i.GetN();
 mass=i.GetMass();
 symb=i.symb;
}

int Isotope::SetIsotope(int _N, int _Z)
{
 FindIsotope(_N,_Z,0,"");
 return 1;
}

int Isotope::SetIsotope(string _Name)
{
 int a;
 string s;
 string s_a, s_s;
 isOK=0;
 if (!(_Name.size()>5)){
 //extract digit and symbol
 for(unsigned int i=0;i!=_Name.size();i++){
  if (isdigit(_Name[i])) s_a=s_a+_Name[i];
  else s_s=s_s+_Name[i];}
 istringstream iss(s_a); iss>>a;
 if (!(s_s.size()>2)) FindIsotope(0,0,a,s_s);
 }
 return 1;
}

string Isotope::GetName()
{
 if(isOK==0) return "-1";
 else
 {
  string tmp(symb);
  tmp=FSTR(A<<symb);
  return tmp;}
}

void Isotope::FindIsotope(int _N, int _Z, int _A, string _Symb)
{
 isOK=0;
 bool found(false);
 bool ok(true);
 char dummy[256];
 int t_A,t_Z, t_N;
 string t_name, tt_name;
 string cc(_Symb);
 double t_mex;
 transform(cc.begin(), cc.end(), cc.begin(), (int(*)(int))std::tolower);

 ifstream *file;
 file = new ifstream("/afs/nd.edu/user30/mcouder/Public/ChargeState/mass.csv");
 file->getline(dummy,256);

 while((!(file->eof()))&&(!found)&&(ok))
 {
  file->getline(dummy,256);
  istringstream iss(dummy);
  if (!(iss >> t_N >> t_Z >> t_A >> t_name >> t_mex)) ok=false;
  else
   if ((_N!=0)&&(_Z!=0)) {if ((t_N==_N)&&(t_Z==_Z)) found=true;}
   else
    if((_A!=0)&&(_Symb.c_str()!="")){
      tt_name=t_name;
      transform(tt_name.begin(), tt_name.end(), tt_name.begin(), (int(*)(int))std::tolower);
     if ((t_A==_A)&&(tt_name==cc)) found=true;}
 }
 if (!found)cout << "The ion has not been found in the file " << endl;
 else
 if (!ok) cout << "There was a problem during the reading of the file."<< endl;
 else
 if(ok&&found){
  A=t_A; Z=t_Z; N=t_N; mex=t_mex/1000.; symb=t_name;
  mass=MeV_amu*A+mex;
  isOK=1;}
 else {cout << "An unknown problem occur during the parsing of the mass file" << endl;}
 file->close();
return;
}

