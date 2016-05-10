#include <iostream>
#include <math.h>

#include "chargemodel.h"
#include "toi.h"
#include "constant.h"

/** \file chargemodel.cpp
 * \brief File in which the charge distribution class are implemented.
 */

using namespace std;

double ChargeStateDist::Getfq(double z,double E)
{
 if (E==energy&&CompareCritical()){;}
 else
 {
  fqvector.clear();
  if (GetIon().GetCheck())
  {
   double sum;
   sum=0.;
   for (int i=0;i<=GetIon().GetZ();i++){
    fqvector.push_back(GetUnNormalizefq(double(i),E));
    if (fqvector[i]<0.) continue;
    sum = sum + fqvector[i];
   }
   for (int i=0;i<=GetIon().GetZ();i++)fqvector[i]=fqvector[i]/sum;
   energy=E;
   FillCritical();
  }
  else {cout << "The ion is not set -- Sorry. " << endl;return 0.;}
 }
 return fqvector[int(z)];
}

double ChargeStateDist::GetMaxq(double E)
{
 double fq,q,max;
 max=0; q=0.;
 for (int i=0;i<=GetIon().GetZ();i++){
   fq=Getfq(double(i),E);
   if (fq<0.) return -1.;
   if (fq>max) {q=double(i); max=fq;}
  }
 return q;
}

bool ChargeStateDist::FillCritical()
{
 criticaldata = GetCriticalData();
 return true;
}

bool ChargeStateDist::CompareCritical()
{
 std::vector<float> tmp = GetCriticalData();
 if (tmp.size()==criticaldata.size()){
  for(unsigned int i=0; i<tmp.size();i++) if (tmp[i]!=criticaldata[i]) return false;}
 else return false;
 return true;
}

//----------------------LIU-------------------------------------
Charge_Liu::Charge_Liu(Isotope i, string s,int j):ChargeStateDist(i)
{
 SetGas(s);
 if(j)WriteReference();
}

int Charge_Liu::SetGas(string s)
{
  string h("h"), H("H");
  string he("he"), He("He"), HE("HE");
  if ((s==H)||(s==h)){
   Gas="H";
   par[0]=1.4211;
   par[1]=0.44515;
   par[2]=0.4495;
   GasSet=1;
  }
  else if ((s==He)||(s==he)||s==HE){
   Gas="He";
   par[0]=1.1326;
   par[1]=0.44515;
   par[2]=0.3449;
   GasSet=1;
  }
  else{
   GasSet=0;
   cout << "The choosen gas is not valid -- Sorry. " << endl;
  }
  return GasSet;
}

double Charge_Liu::GetqBar(double E)
{
 if (GetIon().GetCheck()&&GetGasCheck())
 {
 double fitval;
 fitval=-1.*par[0]*sqrt(E/(GetIon().GetA()*0.067635))/pow(GetIon().GetZ(),par[1]);
 fitval=fitval+par[2];
 fitval=exp(fitval);
 fitval=1.-fitval;
 fitval=GetIon().GetZ()*fitval;

 return fitval;
 }
 else {cout << "The gas is not set -- Sorry. " << endl;return 0.;}
}

double Charge_Liu::GetUnNormalizefq(double z,double E)
{
 if (GetIon().GetCheck()&&GetGasCheck())
 {
 double qbar;
 double sigma;

 qbar=GetqBar(E);

 sigma = 0.23675*pow(GetIon().GetZ(),0.54772);
 if (sigma == 0) return 1.e30;
 double arg = (z-qbar)/sigma;
 double res = exp(-0.5*arg*arg);
 return res/(2.50662827463100024*sigma); //sqrt(2*Pi)=2.50662827463100024
 }
 else {cout << "The gas is not set -- Sorry. " << endl;return 0.;}
}

void Charge_Liu::WriteReference()
{
 cout << "---------------------------------------------------------------"<< endl;
 cout << "You are using: " << endl;
 cout << "Charge state studies of low energy heavy ions passing through hydrogen and helium gas."<<endl;
 cout << "W Liu et. al" << endl;
 cout << "Nuclear Instruments and Methods in Physics Research A496 (2003) 198-214" << endl;
 cout << "---------------------------------------------------------------"<< endl;
}

const std::vector<float> Charge_Liu::GetCriticalData()
{
 std::vector<float> vec(3);
 vec[0]=GetIon().GetN();
 vec[1]=GetIon().GetZ();
 if (GetGas()=="H") vec[2]=1.;
 else if (GetGas()=="He") vec[2]=2.;
 else vec[2]=-1.;
 return vec;
}

//----------------------SAYER-------------------------------------
Charge_Sayer::Charge_Sayer(Isotope i, string s, int j):ChargeStateDist(i)
{
 SetMedium(s);
 if(j)WriteReference();
}

int Charge_Sayer::SetMedium(string s)
{
  string d("solid"), D("Solid");
  string g("gas"), G("Gas");
  if ((s==D)||(s==d)){
   Medium="Solid";
   par[0]=47.3;
   par[1]=0.380;
   par[2]=0.860;
   par[3]=1.03;
   par[4]=0.48;
   par[5]=0.45;
   par[6]=0.26;
   par[7]=0.0007;
   par[8]=0.;
   par[9]=-0.7;
   MediumSet=1;
  }
  else if ((s==g)||(s==G)){
   Medium="Gas";
   par[0]=80.1;
   par[1]=0.506;
   par[2]=0.996;
   par[3]=1.08;
   par[4]=0.35;
   par[5]=0.55;
   par[6]=0.27;
   par[7]=0.17;
   par[8]=0.0012;
   par[9]=-3.3;
   MediumSet=1;
  }
  else{
   MediumSet=0;
   cout << "The choosen medium is not valid -- Sorry. " << endl;
  }
  return MediumSet;
}

double Charge_Sayer::GetqBar(double E)
{
 if (GetIon().GetCheck()&&GetMediumCheck())
 {
   double fitval = sqrt((E*E)+2*E*GetIon().GetMass())/(E+GetIon().GetMass());
   fitval=pow(fitval,par[2]);
   fitval=par[0]*pow(GetIon().GetZ(),-1.*par[1])*fitval;
   fitval=exp(-1*fitval);
   fitval=1.-(par[3]*fitval);
   fitval=GetIon().GetZ()*fitval;
   return fitval;
 }
 else {cout << "The medium and/or the ion are not set -- Sorry. " << endl;return 0.;}
}
double Charge_Sayer::GetUnNormalizefq(double z,double E)
{
 if (GetIon().GetCheck()&&GetMediumCheck())
 {
  double rho;
  double qz;
  double qbar;
  qbar=GetqBar(E);
  rho=0.;
  qz = qbar/GetIon().GetZ();
  if (qz<=0.) rho=0.;
  else{
  rho =pow(qz*(1.-qz),par[6]);
  rho=par[4]*rho*pow(GetIon().GetZ(),par[5]);}
  double fitval ;
  fitval=0.;
  if (rho==0.) return fitval;
  double t= (z-qbar)/rho;
  if (t==0.) return fitval;
  double eps = rho * (par[7]+par[8]*GetIon().GetZ()+par[9]*(sqrt((E*E)+2*E*GetIon().GetMass())/(E+GetIon().GetMass())));
  fitval =(-0.5*t*t)/(1+(eps*t));
  if ((eps*t)<-1.08)return 0.;
  if (fitval>50) return 0.;
  fitval=exp(fitval);
  return fitval;
 }
 else {cout << "The medium and/or the ion are not set -- Sorry. " << endl;return 0.;}
}

void Charge_Sayer::WriteReference()
{
 cout << "---------------------------------------------------------------"<< endl;
 cout << "You are using: " << endl;
 cout << "Semi-empirical formulas for heavy-ion stripping data."<<endl;
 cout << "R.O. Sayer" << endl;
 cout << "Revue de Physique Appliquee, 12(1977) 1543" << endl;
 cout << "---------------------------------------------------------------"<< endl;
}

const std::vector<float> Charge_Sayer::GetCriticalData()
{
 std::vector<float> vec(3);
 vec[0]=GetIon().GetN();
 vec[1]=GetIon().GetZ();
 if (GetMedium()=="Solid") vec[2]=1.;
 else if (GetMedium()=="Gas") vec[2]=2.;
 else vec[2]=-1.;
 return vec;
}

//----------------------Zdist-------------------------------------
double Charge_Zdist::GetqBar(double E)
{
 double vel;
 double qmean;

 vel=13.89E8*sqrt(E/(GetIon().GetMass()/MeV_amu));
 qmean=3.6E08 *pow(GetIon().GetZ(),0.45);
 qmean=pow(vel/qmean,1./0.6);
 qmean=(1./qmean)+1.;
 qmean=GetIon().GetZ()/pow(qmean,0.6);
 return qmean;
}

double Charge_Zdist::GetUnNormalizefq(double z,double E)
{
  double d,probz;
  d = 0.38*pow(GetIon().GetZ(),0.40);
  probz =-1.*pow((z-GetqBar(E)),2.)/(2.*d*d);
  probz = exp(probz)/(d*2.50662827463100024); //sqrt(2*Pi)=2.50662827463100024
  return probz;
}

void Charge_Zdist::WriteReference()
{
 cout << "--------------------------------------------------------------------------------------"<< endl;
 cout << "You are using: " << endl;
 cout << "On the equilibrium charge distribution in heavy elements beams." << endl; 
 cout << "Nikolaev and Dmitriev" << endl;
 cout << "Physics Letters 28A, 277 (1968)." << endl;
 cout << "--------------------------------------------------------------------------------------"<< endl;
}

const std::vector<float> Charge_Zdist::GetCriticalData()
{
 std::vector<float> vec(2);
 vec[0]=GetIon().GetN();
 vec[1]=GetIon().GetZ();
 return vec;
}

//----------------------Dmitriev-------------------------------------
double Charge_Dmitriev::GetqBar(double E)
{
 double l, alpha1, alpha2, rm, rn, emb, clu, vb, rnza, vzam, qmean;
 l=c*100.;
 alpha1=0.1;
 alpha2=0.6;
 rm=1.2;
 rn=5.0;

 emb=E/GetIon().GetA();
 clu=l*1.0e-8;
 vb=clu*sqrt(1.0-1.0/pow((emb/MeV_amu+1.0),2));
 vb=vb/1.6;
 rnza=rn* pow(GetIon().GetZ(),alpha2);
 vzam=vb*pow(GetIon().GetZ(),alpha1)/rm;
 qmean=GetIon().GetZ()*log10(vzam)/log10(rnza);
 return qmean;
}

double Charge_Dmitriev::GetUnNormalizefq(double z,double E)
{
 double sd, sd0, psi,probz,qmean;

 sd0=0.38;
 psi=0.4;
 sd=sd0*pow(GetIon().GetZ(),psi);
 qmean=GetqBar(E);
 probz=-1.*pow((z-qmean),2.)/(2.*pow(sd,2.));
 probz=exp(probz)/(sd*2.50662827463100024); //sqrt(2*Pi)=2.50662827463100024
 return probz;
}

void Charge_Dmitriev::WriteReference()
{
 cout << "---------------------------------------------------------------"<< endl;
 cout << "You are using: " << endl;
 cout << "Gaussian charge state distribution using the semi-empirical method." << endl; 
 cout << "Dmitriev and Nikolaev " << endl;
 cout << "Soviet Physics JETP 20, 409 (1965)." << endl;
 cout << "---------------------------------------------------------------"<< endl;
}

const std::vector<float> Charge_Dmitriev::GetCriticalData()
{
 std::vector<float> vec(2);
 vec[0]=GetIon().GetN();
 vec[1]=GetIon().GetZ();
 return vec;
}

//----------------------Baudinet-------------------------------------
double Charge_Baudinet::GetqBar(double E)
{
 double qmean;
 qmean=GetIon().GetZ()*(1.-exp(-3.8585*sqrt(E*MeV_amu/GetIon().GetMass())/pow(GetIon().GetZ(),0.447)));
 return qmean;
}
double Charge_Baudinet::gamma(double nu)
{
 double gamma;
 if (nu>0) {
 // Coefficients for the series expansion
 double c[7] = { 2.5066282746310005, 76.18009172947146, -86.50532032941677
                 ,24.01409824083091,  -1.231739572450155, 0.1208650973866179e-2
                 ,-0.5395239384953e-5};
 double x   = nu;
 double y   = x;
 double tmp = x+5.5;
 tmp = (x+0.5)*log(tmp)-tmp;
 double ser = 1.000000000190015;
 for (int i=1; i<7; i++) {
    y   += 1;
    ser += c[i]/y;
 }
 gamma = tmp+log(c[0]*ser/x);
 }
 else gamma=0.;
 return exp(gamma);
}

double Charge_Baudinet::GetUnNormalizefq(double z,double E)
{
 double gam,nu;
 double probz,cc,s,v,t;

 v = sqrt(2*E/GetIon().GetMass())*c*1E-6 ;
 if( v/(3.6*pow(GetIon().GetZ(),0.45))<1.)
 {
  cout << "The velocity of " <<  GetIon().GetName() <<   " is to low for this model -- sorry." << endl;
  return -1.;
 }
 s=pow(GetIon().GetZ(),0.4)*(0.426-0.0571/3.6*v*pow(GetIon().GetZ(),(-0.45)));
 cc=2.*(-1.*GetqBar(E)+GetIon().GetZ()+2.)/pow(s,2.);

 nu=2.0*(pow((-1.*GetqBar(E)+GetIon().GetZ()+2.0),2.0))/pow(s,2.);
 gam=gamma(nu/2.);

  t=cc*(GetIon().GetZ()-z+2.);
  probz=cc/(pow(2.0,(nu/2.0))*gam)*(pow(t,(nu/2.0-1.0)))*exp(-t/2.0);
 return probz;
}

void Charge_Baudinet::WriteReference()
{
 cout << "---------------------------------------------------------------"<< endl;
 cout << "You are using: " << endl;
 cout << "Equilibrium charge state distribution using the semi-empirical asymmetric method." << endl;
 cout << "Baudinet-Robinet"<<endl;
 cout << "Nucl. Inst. and Meth. 190, 197 (1981)." <<endl;
 cout << "---------------------------------------------------------------"<< endl;
}

const std::vector<float> Charge_Baudinet::GetCriticalData()
{
 std::vector<float> vec(2);
 vec[0]=GetIon().GetN();
 vec[1]=GetIon().GetZ();
 return vec;
}
