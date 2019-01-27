#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

//this program relies on both spin states having the same energy axis
//the merged file should also be generated in this manner. See merge_DOS.cpp
int main(){
	ifstream infile("Co_total_pdos_l5.dat");//make more general?
	string line;
	double s1,s2;
	getline(infile, line);
	istringstream iss(line);
	iss >> s1;
	getline(infile, line);
	istringstream iss2(line);
	iss2 >> s2;

	double DOS = 0; //initialise
	double E = 0; //initialise
	double step = s2 - s1;
	double a, b;
	double group = 9; //Co is in group 9
	double no_atoms = 1; // There are 5 atoms in the unit cell
	double Fermi = group*no_atoms;
	//This to calculate Fermi level
	while (E < Fermi) 
	{
		getline(infile, line);
		istringstream iss(line);
		if (!(iss >> a >> b)) {break;}
		DOS+=b;
		E = DOS*step;
	}
	cout<<"Fermi level = "<<a<<"eV"<<endl;
	cout<<"Termination energy = "<<E<<"eV, against reference "<<Fermi<<"eV"<<endl;
	ifstream infile2("Co_up_pdos_l5.dat");//make more general?
	ifstream infile3("Co_dn_pdos_l5.dat");//make more general?
	double E1, E2, DOS1, DOS2;
	DOS1 = 0;
	DOS2 = 0;
	double fermi_level = a;
	a = -23;//this guarantees a is smaller than fermi_level when the loop starts
	//Now we intend to find what proportion of atoms come from which spin state
	while (a < fermi_level) 
	{
		getline(infile2, line);
		istringstream iss(line);
		iss >> a >> b;
		DOS1+=b;
		E1 = DOS1*step;
		getline(infile3, line);
		istringstream iss2(line);
		iss2 >> a >> b;
		DOS2+=b;
		E2 = DOS2*step;
	}
	// 2 in denominator for Wannier only, as it is multiplied by 2...
	cout<<"spin up has on average "<<E1/(2.*no_atoms)<<" electrons"<<endl;
	cout<<"spin dn has on average "<<E2/(2.*no_atoms)<<" electrons"<<endl;

	return 0;
}
