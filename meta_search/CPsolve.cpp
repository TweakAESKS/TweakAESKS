#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

#include "CPsolve.hpp"

using namespace std;

string exec(const char* cmd){
//Exec cmd and grab the stdout
    array<char, 128> buffer;
    string result;
    shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

void createPermFile(vector<uint8_t> const & Pk, vector<uint8_t> const & Ps, uint const Nround){
	ofstream file;
	file.open("dataPerm_MCeq.dzn");
	file << "Pk = array1d(0..15, [";
	for(auto x: Pk){
		file << uint(x) << ",";
	}
	file << "]);" << endl;
  file << "Ps = array1d(0..15, [";
	for(auto x: Ps){
		file << uint(x) << ",";
	}
	file << "]);" << endl;
	file << "n = " << Nround << ";" << endl;
	file.close();
}

uint shortestPath(std::vector<uint8_t> const & Pk, std::vector<uint8_t> const & Ps, uint const Nround){
	createPermFile(Pk,Ps,Nround);
	cout << "Searching a path with Pk = ";
	for(auto const & x: Pk)	cout << int(x) << " ";
	cout << endl;
	auto s = exec("mzn-gecode modelPerm_MCeq.mzn dataPerm_MCeq.dzn");

	string delimiter = "\n";
	string token;
	//The first line is the number of active sboxes
	size_t pos = s.find(delimiter);
	token = s.substr(0,pos);
	return stoi(token);

	// if( s == "=====UNSATISFIABLE=====\n"){
	// 	cout << "UNSATISFIABLE problem with permutation ";
	// 	for(auto const & x: Pk) cout << int(x) << " ";
	// 	cout << endl;
	// 	return 0;
	// }
	// else{
	// 	string delimiter = "\n";
	// 	string token;

	// 	//The first line is the number of active sboxes
	// 	size_t pos = s.find(delimiter);
	// 	token = s.substr(0,pos);
	// 	return stoi(token);
	// }
}