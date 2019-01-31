#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

//this program relies on both spin states having the same energy axis
//the merged file should also be generated in this manner. See merge_DOS.cpp
int main(){
	ifstream infile("Co_field_total_ldos.dat");//make more general?
	string line;
	double s1,s2;
	getline(infile, line);
	istringstream iss(line);
	iss >> s1;
	getline(infile, line);
	istringstream iss2(line);
	iss2 >> s2;

	double DOS = 0; //initialise
	double E = 0.; //initialise
	double step = s2 - s1;
	double a, b;
	double group = 9; //Co is in group 9
	double no_atoms = 5; // There are 5 atoms in the unit cell
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
	ifstream infile2("Co_field_up_ldos.dat");//make more general?
	ifstream infile3("Co_field_dn_ldos.dat");//make more general?
	ifstream infile4("Co_field_total_pdos_l1.dat");//make more general?
	/* ifstream infile5("Co_field_total_pdos_l2.dat");//make more general? */
	/* ifstream infile6("Co_field_total_pdos_l3.dat");//make more general? */
	/* ifstream infile7("Co_field_total_pdos_l4.dat");//make more general? */
	ifstream infile8("Co_field_total_pdos_l5.dat");//make more general?
	double E1, E2, E3, E4, E5, E6, E7, DOS1, DOS2, DOS3, DOS4, DOS5, DOS6, DOS7;
	DOS1 = 0; DOS2 = 0; DOS3 = 0; DOS4 = 0; DOS5 = 0; DOS6 = 0; DOS7 = 0;
	double fermi_level = a;
	double e,f,g;
	a = -23;//this guarantees a is smaller than fermi_level when the loop starts
	//Now we intend to find what proportion of atoms come from which spin state
	while (a <= fermi_level) 
	{
		getline(infile2, line);
		istringstream iss(line);
		iss >> a >> b;
		DOS1+=b;
		E1 = DOS1*step;
	}
	e = -23;
	while (e <= fermi_level) 
	{
		getline(infile3, line);
		istringstream iss(line);
		iss >> e >> b;
		DOS2+=b;
		E2 = DOS2*step;
	}
	f = -23;
	while (f <= fermi_level) 
	{
		getline(infile4, line);
		istringstream iss(line);
		iss >> f >> b;
		DOS3+=b;
		E3 = DOS3*step;
	}
	/* a = -23; */
	/* while (a < fermi_level) */ 
	/* { */
	/* 	getline(infile5, line); */
	/* 	istringstream iss(line); */
	/* 	iss >> a >> b; */
	/* 	DOS4+=b; */
	/* 	E4 = DOS4*step; */
	/* } */
	/* a = -23; */
	/* while (a < fermi_level) */ 
	/* { */
	/* 	getline(infile6, line); */
	/* 	istringstream iss(line); */
	/* 	iss >> a >> b; */
	/* 	DOS5+=b; */
	/* 	E5 = DOS5*step; */
	/* } */
	/* a = -23; */
	/* while (a < fermi_level) */ 
	/* { */
	/* 	getline(infile7, line); */
	/* 	istringstream iss(line); */
	/* 	iss >> a >> b; */
	/* 	DOS6+=b; */
	/* 	E6 = DOS6*step; */
	/* } */
	g = -23;
	while (g <= fermi_level) 
	{
		getline(infile8, line);
		istringstream iss(line);
		iss >> g >> b;
		DOS7+=b;
		E7 = DOS7*step;
	}
	// 2 in denominator for Wannier only, as it is multiplied by 2...
	cout<<"spin up has on average "<<E1/(2.*no_atoms)<<" electrons"<<endl;
	/* cout<<a<<endl; */
	cout<<"spin dn has on average "<<E2/(2.*no_atoms)<<" electrons"<<endl;
	/* cout<<e<<endl; */
	cout<<"layer one contains "<<E3<<" electrons"<<endl;
	/* cout<<f<<endl; */
	/* cout<<"layer two contains "<<E4<<" electrons"<<endl; */
	/* cout<<"layer three contains "<<E5<<" electrons"<<endl; */
	/* cout<<"layer four contains "<<E6<<" electrons"<<endl; */
	cout<<"layer five contains "<<E7<<" electrons"<<endl;
	/* cout<<g<<endl; */

	return 0;
}
