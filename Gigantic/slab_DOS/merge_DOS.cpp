#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[]){
	string file1 = argv[1];
	string file2 = argv[2];

	ifstream infile1(file1);
	ifstream infile2(file2);
	string line;
	double a, b, c, d;
	string Mydata = "merge_DOS.dat";
	ofstream Myfile;	
	Myfile.open( Mydata.c_str(),ios::trunc );
	while (!infile1.eof()){
		getline(infile1, line);
		istringstream iss(line);
		iss >> a >> b;
		getline(infile2, line);
		istringstream iss2(line);
		iss2 >> c >> d;
		Myfile<<(a+c)/2.<<" "<<(b+d)/2.<<endl; //divide by two on the second here when using Wannier
	}
	Myfile.close();
	return 0;
}
