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

void createPermFile(vector<uint8_t> const & Pk, vector<uint8_t> const & Ps, unsigned const Nround, unsigned const bound){
	ofstream file;
	file.open("dataPerm.dzn");
	file << "Pk = array1d(0..15, [";
	for(auto x: Pk){
		file << unsigned(x) << ",";
	}
	file << "]);" << endl;
  file << "Ps = array1d(0..15, [";
	for(auto x: Ps){
		file << unsigned(x) << ",";
	}
	file << "]);" << endl;
	file << "objective_bound = " << bound-1 << ";" << endl;
	file << "n = " << Nround << ";" << endl;
	file.close();
}


bool shortestPathAndState(vector<uint8_t> const & Pk, vector<uint8_t> const & Ps, unsigned const Nround, unsigned const bound){

	createPermFile(Pk,Ps,Nround,bound);

	string str_cmd = "mzn-gecode modelPerm.mzn dataPerm.dzn";
	const char* cmd = str_cmd.c_str();
	auto s = exec(cmd);

	if( s == "=====UNSATISFIABLE=====\n") return true;
  else return false;
}
