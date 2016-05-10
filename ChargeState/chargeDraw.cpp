#include <string>
#include <sstream>

#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>

#include "toi.h"
#include "chargemodel.h"

/** \file chargeDraw.cpp
 * \brief An example for drawing charge state.
 * 
 * This file depend on the ROOT framework. It display 2 canvas.
 * The first one is on plot of the fraction of each charge state
 * of 12C between E=1 MeV and E=15 MeV. The second one is a plot of
 * the maximum charge state in function of the energy.
 *
 * \section trick Trick
 * Despite the fact the ChargeStateDist provide a double returning the 
 * fraction of the charge at a given energy we have two define a function
 * fq which will use the global ChargeStateDist *t to return the wanted value.
 * This is because the TF1 class of ROOT can not directly take a class member.
 * The same apply to maxq.
 */

using namespace std;

//! Global pointer to a ChargeStateDist see section trick in the file description.
ChargeStateDist *t;

/*! Return the charge state fraction. This is required to plot the Getfq member using the TF1 class from ROOT, see section trick in the file description. */
double fq(double *x, double *par)
{
 //x[0] = kinetic energy
 //par[0] = selected charge state
 double fitval;
 fitval=t->Getfq(par[0],x[0]);
 return fitval;
}

/*! Return the maximum charge state. This is required to plot the GetMaxq member using the TF1 class from ROOT, see section trick in the file description. */
double maxq(double *x, double *par)
{
 //x[0] = kinetic energy
 double fitval;
 fitval=t->GetMaxq(x[0]);
 return fitval;
}

int GetModel()
{
 int choice;
 bool isOK(false);
 choice=0;
 cout << "You have the choice between those different model:" << endl;
 cout << "1: Baudinet-Robinet (qdist3 at nsl)" << endl;
 cout << "2: Dmitriev and Nikolaev for solid 1965 (qdist2 at nsl)" << endl;
 cout << "3: Nikolaev and Dmitriev 1968 (zdist at nsl)" << endl;
 cout << "4: Sayer for solid"<< endl;
 cout << "5: Sayer for dilute gas" << endl;
 cout << "6: Liu for H" << endl;
 cout << "7: Liu for He" << endl;
 while (((choice>0)||(choice<8))&&(!isOK)){
 cout << "Please enter your choice : " ;
 cin  >> choice;
 if ((choice<1)||(choice>7)) cout << "Not valid choice" <<endl;
 else isOK=true;
 }
 return choice;
}

int main(int argc, char **argv)
{
 string ion_name;
 Isotope ion;
 bool isOK(false);
 string s_min_energy;
 string s_max_energy;
 double min_energy;
 double max_energy;
 double max;
 int store_max;

 TApplication theApp("App", &argc, argv);

   if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
   }
 while (!isOK){
  cout << "Enter ion name (e.g. 12C) : ";
  cin >> ion_name;
  ion.SetIsotope(ion_name);
  if (ion.GetCheck()) isOK=true;
  else cout << "This is not an ion name" << endl;
 }isOK=false;
 while (!isOK){
  cout << "Enter the min energy (MeV) : " ;
  cin >> s_min_energy ;
  istringstream min_stream(s_min_energy);
  if (min_stream >> min_energy) isOK=true;
  else cout << "This is not a valid energy" << endl;
 }
 isOK=false;
 while (!isOK){
  cout << "Enter the max energy (MeV) : " ;
  cin >> s_max_energy;
  istringstream max_stream(s_max_energy);
  if (max_stream >> max_energy) isOK=true;
  else cout << "This is not a valid energy" << endl;
 }
 switch (GetModel())
 {
  case 1:
       t = new Charge_Baudinet(ion);
       break;
  case 2:
       t = new Charge_Dmitriev(ion);
       break;
  case 3:
       t = new Charge_Zdist(ion);
       break;
  case 4:
       t = new Charge_Sayer(ion,"solid");
       break;
  case 5:
       t = new Charge_Sayer(ion,"gas");
       break;
  case 6:
       t = new Charge_Liu(ion,"H");
       break;
  case 7:
       t = new Charge_Liu(ion,"He");
       break;
  default:
       cout << "It seems that you trigger a problem that was not forseen" << endl;
       exit(-1);
       break;
 }

 TCanvas *c1;
 c1 = new TCanvas("c1");

 max=0.; store_max=0;
 TF1 *arrayfq[ion.GetZ()+1];
 for (int i=0;i<ion.GetZ()+1;i++){
  arrayfq[i] = new TF1(Form("fq%d",i),fq,min_energy,max_energy,1);
  arrayfq[i]->SetParameter(0,double(i));
  arrayfq[i]->SetLineColor(i+1);
  if (arrayfq[i]->GetMaximum(min_energy,max_energy)>max) {store_max=i; max=arrayfq[i]->GetMaximum(min_energy,max_energy);}
 }

 arrayfq[store_max]->GetHistogram()->SetTitle(Form("^{%d}%s charge state distribution",ion.GetA(),(ion.GetSymbol()).c_str()));
 arrayfq[store_max]->GetHistogram()->GetXaxis()->SetTitle("Energy (MeV) ");
 arrayfq[store_max]->GetHistogram()->GetXaxis()->SetTitleFont(42);
 arrayfq[store_max]->Draw();
 for (int i=0;i<ion.GetZ()+1;i++) arrayfq[i]->Draw("SAME");

 TLegend *legend=new TLegend(0.6,0.65,0.88,0.95);
 legend->SetTextFont(72);
 legend->SetTextSize(0.04);
 for (int i=0;i<ion.GetZ()+1;i++)legend->AddEntry(arrayfq[i],Form("Charge %d",i),"l");
 legend->Draw();
 /*TCanvas *c2;
 c2 = new TCanvas("c2");
 TF1 *g = new TF1("maxq",maxq,1.,15.,0);
 g->Draw(); */
 theApp.Run();
 return 1;
}
