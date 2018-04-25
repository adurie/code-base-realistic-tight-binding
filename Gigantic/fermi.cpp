#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(){
	ifstream infile("Ag_spin_down.txt");
	/* ifstream infile("Fe_tdos.txt"); */
	string line;
	double DOS = 0;
	double E = 0;
	double step = 0.0026;
	double a, b;
	/* while (E < 8) */ 
	while (E < 5.5) 
	{
		getline(infile, line);
		istringstream iss(line);
		if (!(iss >> a >> b)) {break;}
		DOS+=b;
		E = DOS*step;
	}
	cout<<a<<endl;
	cout<<E<<endl;
	return 0;
}
