#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
//#include <cstdlib>


#include "toi.h"
#include "chargemodel.h"

using namespace std;

/** \file ChargeDist.cpp
 * \brief Program which calculate the charge state distribution.
 *
 * This program asks for the name of the ion you want to use. It asks
 * for the model you want to use. It displays some properties of your
 * calculation and the charge state fraction for all the charge which
 * represent more than 0.01%. It finally displays the average charge state.
 *
 */

using namespace std;

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

void printarray(double a[][2], int length);

int main()
{
cout<< "Single Stripper or Double Stripper (type 's' or 'd'): "<< endl;
string answer;
cin>>answer;

if (answer=="d"){

string ion_name;
string s_energy;
bool isOK(false);
Isotope ion;
double energy;
double average;
double fq;

while (!isOK){
cout << "Enter ion name (e.g. 12C) : ";
cin >> ion_name;
ion.SetIsotope(ion_name);
if (ion.GetCheck()) isOK=true;
else cout << "This is not an ion name" << endl;
}
isOK=false;
while (!isOK){
cout << "Enter the beam energy (MeV) : " ;
cin >> s_energy;
istringstream stream(s_energy);
if (stream >> energy) isOK=true;
else cout << "This is not a valid energy" << endl;
}
ChargeStateDist *t;
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
       //exit(-1);
       break;
}

cout << endl << "---------------------------------------------------------------------"<<endl;
cout << "Inital Energy of Beam (MeV) = " << energy << endl;
cout << "Atomic Number (Z) = " << ion.GetZ() << endl;
cout << "Mass of Ion (amu) =  " << ion.GetA() << endl;
cout << "Maximum of the charge state distribution = " << t->GetqBar(energy) << endl;


cout << endl << "Charge     Percentage (Not printed if <0.01%)" << endl;
cout << "---------------------------------------------------------------------"<<endl;


average=0.; fq=0.;
double array[50][2];
int a=0;
int b=0;
int length=0;
ofstream arrayData("array.txt");
arrayData<< "Beam Energy (MeV) - at First Stripper = " <<"	"<<energy<< endl;
for (int i=0; i<=ion.GetZ();i++)
{
fq=t->Getfq(double(i),energy);
if (fq<0.) break;
average = average+double(i)*fq;
if (fq>0.0001){
cout << "  " <<i<<"       "<<100.* fq << endl;
array[a][b]=i;
b++;
array[a][b]=100.*fq;
b=0;
a++;
length++;
arrayData<<i<< "	"<<100*fq<<endl; }

}
if (fq>=0.){
//cout << "---------------------------------------------------------------------"<<endl;
//cout << "Average charge state = " << average << endl;
//cout << "---------------------------------------------------------------------"<<endl;
}

//printarray(array, length);
double energy2=0;
//cout<< "New energies for allowed charge states at second stripper"<< endl;
for (int n=0; n<length; n++){
	energy2 = energy + (energy/2)*array[n][0];
	//cout<< array[n][0] << "	" << energy2<< endl;
	//cout<< "--------------------------------------------------"<< endl;
	//cout<< "Energy of Beam (MeV) = " << energy2 << endl;
	//cout<< "Maximum of the charge state distribution = " << t->GetqBar(energy2) << endl;
	double average2=0.; double fq2=0.;
	double array2[50][2];
	int a2=0;
	int b2=0;
	int length2 =0;
	
arrayData<<"--------------------" << endl;
arrayData<< "Beam Energy (MeV)at Second Stripper = "<< energy2 << ", Charge State = " << array[n][0]<< endl;
	for (int m=0; m<=ion.GetZ();m++)
{
fq2=t->Getfq(double(m),energy2);
if (fq2<0.) break;
average2 = average2+double(m)*fq2;
if (fq2>0.0001){
//cout << "  " <<m<<"       "<<100.* fq2 << endl;
array2[a2][b2]=m;
b2++;
array2[a2][b2]=100.*fq2;
b2=0;
a2++;
length2++;
double energyHE;
energyHE = energy2 + (energy/2)*m;
arrayData<<m<< "	"<<100*fq2<<"	"<<" Energy at HE cup: "<< energyHE << endl; }
}
}
if (fq>=0.){
//cout << "---------------------------------------------------------------------"<<endl;
//cout << "Average charge state = " << average << endl;
//cout << "---------------------------------------------------------------------"<<endl;
}
}else if (answer=="s"){

	string ion_name;
string s_energy;
bool isOK(false);
Isotope ion;
double energy;
double average;
double fq;

while (!isOK){
cout << "Enter ion name (e.g. 12C) : ";
cin >> ion_name;
ion.SetIsotope(ion_name);
if (ion.GetCheck()) isOK=true;
else cout << "This is not an ion name" << endl;
}
isOK=false;
while (!isOK){
cout << "Enter the beam energy (MeV) : " ;
cin >> s_energy;
istringstream stream(s_energy);
if (stream >> energy) isOK=true;
else cout << "This is not a valid energy" << endl;
}
ChargeStateDist *t;
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
       //exit(-1);
       break;
}
cout << endl << "---------------------------------------------------------------------"<<endl;
cout << "Inital Energy of Beam (MeV) = " << energy << endl;
cout << "Atomic Number (Z) = " << ion.GetZ() << endl;
cout << "Mass of Ion (amu) =  " << ion.GetA() << endl;
cout << "Maximum of the charge state distribution = " << t->GetqBar(energy) << endl;


cout << endl << "Charge     Percentage (Not printed if <0.01%)" << endl;
cout << "---------------------------------------------------------------------"<<endl;
average=0.; fq=0.;

double array[50][2];
int a=0;
int b=0;
int length=0;
ofstream arrayData("array.txt");
arrayData<< "Beam Energy (MeV) - at First Stripper = " <<"	"<<energy<< endl;
for (int i=0; i<=ion.GetZ();i++)
{
fq=t->Getfq(double(i),energy);
if (fq<0.) break;
average = average+double(i)*fq;
if (fq>0.0001){
cout << "  " <<i<<"       "<<100.* fq << endl;
array[a][b]=i;
b++;
array[a][b]=100.*fq;
b=0;
a++;
length++;
arrayData<<i<< "	"<<100*fq<<endl; }

}
if (fq>=0.){
//cout << "---------------------------------------------------------------------"<<endl;
//cout << "Average charge state = " << average << endl;
//cout << "---------------------------------------------------------------------"<<endl;
}
}






return 1;
}

void printarray(double a[][2], int length){
	for(int x=0; x<length; x++){
	cout<< a[x][0]<< " ";
	cout<< a[x][1]<< " ";
	cout<< "\n";
	}
}
