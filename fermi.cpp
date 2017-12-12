#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(){
	/* ifstream infile("tdosrem2.txt"); */
	/* ifstream infile("rem.txt"); */
	ifstream infile("diamond.txt");
	/* ifstream infile("mn.txt"); */
	/* ifstream infile("gsl.txt"); */
	/* ifstream infile("mnm.txt"); */
	/* ifstream infile("2mmm.txt"); */
	string line;
	double DOS = 0;
	double E = 0;
	/* double step = 0.001; */
	double step = 0.0026;
	/* double step = 0.0001; */
	double a, b;
	while (E < 22) 
	/* while (E < 11) */ 
	/* while (E < 11.308) */ 
	/* while (E < 4.8762) */ 
	/* while (E < 5.80353) */ 
	/* while (E < 4.8745) */ 
	/* while (E < 2) */ 
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
