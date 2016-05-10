#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
//#include <cstdlib>
#include <stdlib.h>
#include <locale>
#include <algorithm>


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
 * Modified from Manoel Couder's program to include the option of a 
 * second stripper half way between the terminal and HE end of FN.
 */


//Choice to get difference models.------------------------------------------------------------------
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


//To print array to the screen.------------------------------------------------------------------
void printarray(double a[][2], int length){
	for(int x=0; x<length; x++){
	cout<< a[x][0]<< " ";
	cout<< a[x][1]<< " ";
	cout<< "\n";
	}
}





//Main program.-----------------------------------------------------------------------------------
int main()
{

string ion_name;
string s_energy;
bool isOK(false);
Isotope ion;
double energy;
double average;



while (!isOK){
cout << "Enter ion name (e.g. 12C) : ";
cin >> ion_name;
ion.SetIsotope(ion_name);
if (ion.GetCheck()) isOK=true;
else cout << "This is not an ion name" << endl;
}
isOK=false;
while (!isOK){
cout << "Enter the Terminal Voltage (MV) : " ;
cin >> s_energy;
istringstream stream(s_energy);
if (stream >> energy) isOK=true;
else cout << "This is not a valid terminal voltage" << endl;
}


std::locale loc;
  
  int mass;
  if (isdigit(ion_name[0],loc))
  {
    
    std::stringstream(ion_name) >> mass;
    //std::cout << "The Mass is " << mass << ".\n";
    ion_name.erase(std::remove_if(ion_name.begin(), ion_name.end(), ::isdigit),ion_name.end());
    
    //std::cout<< ion_name << mass<< "\n";
  }
  //std::cout<< mass <<"\n";
  //std::cout<< ion_name<< "\n";

  std::ostringstream flip;
  flip<< ion_name << mass;
  //std::cout << flip.str() <<"\n";


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
cout << "Maximum of the charge state distribution = "<< t->GetqBar(energy)<< endl;
cout << endl << "Charge     Percentage (Not printed if <0.01%)" << endl;
cout << "---------------------------------------------------------------------"<<endl;
cout << "Charge State"<< "	"<<"Probability"<< "	"<< "Energy at HE cup (MeV)"<< endl;

double array[50][2]; //charge states and fq's from first stripper
double array2[50][2];
//double array3[50][2];
ofstream arrayData("array.txt");
double fq=0;
double fq2=0;
//double fq3=0;

double energyHE1=0;//HE energy for 1 stripper
double energy2=0;//Energy at second stripper for one charge state
double energyHE2=0; //HE energy for one charge state with two strippers
double energy3=0;// Energy at second stripper for full charge state set
double energyHE3=0;// HE energy for full charge state set with 2 strippers

int a=0;
int b=0;
int a2=0;
int b2=0;
int length=0;
//int length2=0;
//int length3=0;
double finalprob=0;
double finalprob2=0;

double minj=0;
double mbreak=0;
cout << "Mass injected: " <<endl;
cin >>minj;
cout << "Mass after breakup: "<<endl;
cin >>mbreak;



arrayData<< "Beam Energy (MeV) - at First Stripper = "<< "	"<<energy<< endl;
arrayData<< "Charge State"<< "	"<<"Probability"<< "	"<< "Energy at HE cup (MeV)"<< endl;

for (int i=0; i<=ion.GetZ();i++)
{
fq=t->Getfq(double(i),energy);
if (fq<0.) break;
average = average+double(i)*fq;
if (fq>0.0001){
array[a][b]=i;
b++;
array[a][b]=100.*fq;
b=0;
a++;
length++;
energyHE1 = energy*((mbreak/minj)+i);
cout << "  " <<i<<"       "<<100.* fq << "	"<<energyHE1<<endl;
arrayData<<i<< "	"<<100*fq<<"	"<<energyHE1<< endl; }

}
if (fq>=0.){
//cout << "---------------------------------------------------------------------"<<endl;
//cout << "Average charge state = " << average << endl;
//cout << "---------------------------------------------------------------------"<<endl;
}
cout<< "Single Stripper or Double Stripper (type 's' or 'd'): "<< endl;
string answer;
cin>>answer;
if(answer=="s"){
//DONE
arrayData<< "__________________________________________________________"<<endl;
}else if(answer=="d"){
arrayData<< "__________________________________________________________"<<endl;
arrayData<< "Beam Energy (MeV) - with Second Stripper "<< endl;

cout<<"One charge state or Full Set at 2nd stripper? (type 'o' or 'f')? "<< endl;
string answer2;
cin>>answer2;

	if(answer2=="o"){
	arrayData<< "Charge State"<< "	"<<"Probability"<< "	"<< "Energy at HE cup (MeV)"<< "	"<<"Final Charge State"<<"	"<<"Final Probability"<<endl;
	cout<< "Which charge state would you like to use? " << endl;
	int charge=0;
	cin>>charge;
	cout<< "What charge state do you want after the second stripper? "<< endl;
	int charge2=0;
	cin>>charge2;	

	energy2 = energy*((mbreak/minj)+(charge*0.5));
	fq2=t->Getfq(charge2, energy2);
	energyHE2=energy*((mbreak/minj)+(charge*0.5)+(charge2*0.5));

	for (int j=0; j<=length; j++){
		if(charge == array[j][0]){
		finalprob = array[j][1]*fq2;		
		}
	}
	
	arrayData<<charge<< "	"<<fq2*100<< "	"<<energyHE2<< "	"<<charge2<<"	"<<finalprob <<endl;
	cout<< "	"<<charge2<<"	"<<fq2*100<<"	"<<energyHE2<<"	"<<finalprob<< endl;	

//ion_name string
//charge2 int
//energyHE2 double


string nmr;
stringstream sc2, seHE2;
sc2<<charge2;
seHE2<<energyHE2;
nmr = "echo $'s\\n" + flip.str() + "\\n" + sc2.str() + "\\n" + seHE2.str() + "' | nmrfreq";
cout<< nmr<< endl;
system(nmr.c_str());
	
	}else if (answer2=="f"){
		//printarray(array,6);
		
		for(int n=0; n<length; n++){
		energy3 = energy*(1+(0.5*array[n][0]));
		arrayData<<"Charge State from 1st stripper = "<<array[n][0]<<"	"<<"Energy at Second Stripper = "<<energy3<< endl;
		arrayData<< "Charge State"<< "	"<<"Probability"<< "	"<< "Energy at HE cup (MeV)"<<"	"<<"Final Probability"<<endl;
			for(int m=0; m<=ion.GetZ(); m++){
			fq2=t->Getfq(double(m), energy3);
			if(fq2<0.)break;
			if (fq2>0.0001){
			array2[a2][b2]=m;
			b2++;
			array2[a2][b2]=100*fq2;
			b2=0;
			a2++;
			energyHE3 = energy*((mbreak/minj)+(0.5*array[n][0])+(m*0.5));
			finalprob2 = array[n][1]*fq2; 
			arrayData<<m<< "	"<<100*fq2<< "	"<<energyHE3<<"	"<<finalprob2<<endl;
			finalprob=0;
			}//end of if
			}//end of for (m, GetZ)
		arrayData<< "__________________________________________________________"<<endl;
				
	}//end of for (n, length
	cout<< "Data has been printed to the file array.txt." << endl;
}//end of answer2 = f
}//end of answer = d







return 1;

}//end of main function



