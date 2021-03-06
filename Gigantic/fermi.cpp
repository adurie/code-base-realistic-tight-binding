#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(){
	/* ifstream infile("Ag_spin_down.txt"); */
	/* ifstream infile("Cu_111_spin_down2.txt"); */
	/* ifstream infile("Cu_spin_down.txt"); */
	ifstream infile("wan-tdos.txt");
	string line;
	double DOS = 0;
	double E = 0;
	double step = 0.0026;
	double a, b;
	while (E < 8) 
	/* while (E < 5.5) */ 
	{
		getline(infile, line);
		istringstream iss(line);
		if (!(iss >> a >> b)) {break;}
		DOS+=b;
		E = DOS*step;
	}
	cout<<"Fermi level = "<<a<<"Ry"<<endl;
	cout<<"Termination energy = "<<E<<"Ry"<<endl;
	return 0;
}
