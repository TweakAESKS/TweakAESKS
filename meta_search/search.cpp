#include <vector>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>
#include "CPsolve.hpp"
#include "allPathsNoEquations.hpp"

using namespace std;



std::random_device rd; 
std::mt19937 rng(rd());
std::uniform_int_distribution<int> randint(0,15);
std::uniform_real_distribution<double> randouble(0,1);

vector<uint8_t> epsilon(vector<uint8_t> const & Pk){
	//Generate a neighbour of Pk
	//This is done by switching two elements of Pk randomly
	vector<uint8_t> Pk2(Pk);
	auto i1 = randint(rng);
	auto i2 = randint(rng);
	while(i2 == i1) i2 = randint(rng);
	swap(Pk2[i1], Pk2[i2]);
	return Pk2;
}

double alpha(double const T, double const beta){
	//Cooling schedule
	return T/(1+beta*T);
}

int main(){

	clock_t start = clock();
	double duration;

	DicStateKey dic;
	vector<uint8_t> Pk({5,2,3,8,9,6,7,12,13,10,11,0,1,14,15,4}); //Initial permutation
	PermTrans p = vecToPermTrans(Pk);
	vector<uint8_t> Ps({0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3}); //ShiftRows
	double T = 2.0; //Initial temperature
	double beta = 0.001; //beta is used in the cooling schedule
	uint bound_active = 21;
	uint activeK = 7;
	uint n_rounds = 6;

	
	//Precomputation for the quicksearch

	auto const & state_conf = dic.getStateConfs();
	auto const n_conf = state_conf.size();
	// cout << "start search" << endl;
	/* Generate transitions */
	vector<vector<vector<Subset>>> transitions1r (activeK+1);
	vector<vector<vector<Subset>>> transitions1r_ji (activeK+1);
	//vector<vector<vector<Subset>>> transitions1r_ij (activeK+1);
	for (unsigned a = 1; a <= activeK; ++a) {
		auto transitions1r_a = generateTransitions(bound_active, a, dic, Ps);
		// cout << "\r - generate all transitions (" <<  a << "): 100%" << endl;
		//vector<vector<Subset>> transitions1r_ka_ij (dic.nbKeys(a), vector<Subset> (n_conf));
		vector<vector<Subset>> transitions1r_ka_ji (dic.nbKeys(a), vector<Subset> (n_conf));
		for (unsigned i = 0; i < n_conf; ++i) {
		  for (unsigned j = 0; j < n_conf; ++j) {
		    transitions1r_a[i][j].apply([&transitions1r_ka_ji,i,j](unsigned k){
		      transitions1r_ka_ji[k][j] += i;
		      //transitions1r_ka_ij[k][i] += j;
		    });
		  }
		}
		transitions1r_ji[a] = move(transitions1r_ka_ji);
		//transitions1r_ij[a] = move(transitions1r_ka_ij);
		transitions1r[a] = move(transitions1r_a);
	}

	vector<vector<vector<Subset>>> transitions1r_noeq (activeK+1);
	vector<vector<vector<Subset>>> transitions1r_ji_noeq (activeK+1);
	//vector<vector<vector<Subset>>> transitions1r_ij (activeK+1);
	for (unsigned a = 1; a <= activeK; ++a) {
		auto transitions1r_a = generateTransitionsNoEquation(bound_active, a, dic, Ps);
		// cout << "\r - generate all transitions (" <<  a << "): 100%" << endl;
		//vector<vector<Subset>> transitions1r_ka_ij (dic.nbKeys(a), vector<Subset> (n_conf));
		vector<vector<Subset>> transitions1r_ka_ji (dic.nbKeys(a), vector<Subset> (n_conf));
		for (unsigned i = 0; i < n_conf; ++i) {
		  for (unsigned j = 0; j < n_conf; ++j) {
		    transitions1r_a[i][j].apply([&transitions1r_ka_ji,i,j](unsigned k){
		      transitions1r_ka_ji[k][j] += i;
		      //transitions1r_ka_ij[k][i] += j;
		    });
		  }
		}
		transitions1r_ji_noeq[a] = move(transitions1r_ka_ji);
		//transitions1r_ij[a] = move(transitions1r_ka_ij);
		transitions1r_noeq[a] = move(transitions1r_a);
	}

	ofstream file;
	file.open("log_timing.txt", ios::out | ios::app);
	file << "trans init ok" << endl;

	// uint min = NoPathsNoEq(p, bound_active, n_rounds, transitions1r_ji, transitions1r_ji_noeq, dic);
	// file << min << endl;

	uint64_t nbiter = 0;
	uint64_t bestiter = 0;
	vector<uint8_t> bestPk(16); uint bestfPk = 0;
	uint fPk = NoPathsNoEq(p, bound_active, n_rounds, transitions1r_ji, transitions1r_ji_noeq, dic); //Initial short path length
	while(fPk < bound_active){
		nbiter++;
		vector<uint8_t> Pk2 = epsilon(Pk); //Get a neighbour of Pk
		PermTrans p2 = vecToPermTrans(Pk2);
		file << "Pk = ";
		for(auto const & x: Pk2) file << int(x) << " ";
		file << endl;
		uint fPk2 = NoPathsNoEq(p2, bound_active, n_rounds, transitions1r_ji, transitions1r_ji_noeq, dic); //Get a short path from Pk2
		duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		file << "Length " << fPk2 << "(" << nbiter << "-th, " << duration << "s)" << endl;

		if(fPk2 >= bound_active){ 
			//If the shortest path seems to be of length bound_active, compute the actual shortest path
			//No matter what the actual shortest path length is, choose this permutation for the next iteration
			Pk = move(Pk2);
			fPk = shortestPath(Pk,Ps,n_rounds);
			file << "Real shortest path : " << fPk << endl;
			if(fPk > bestfPk){ //Save the best solution encountered
				bestfPk = fPk;
				bestPk = Pk;
				bestiter = nbiter;
			}
		}
		else if(fPk2 > fPk){
			Pk = move(Pk2);
			fPk = fPk2;
			if(fPk > bestfPk){ //Save the best solution encountered
				bestfPk = fPk;
				bestPk = Pk;
				bestiter = nbiter;
			}
		}
		else{
			double r = randouble(rng);
			if(r < exp((fPk2-fPk)/T)){
				Pk = move(Pk2);
				fPk = fPk2;
			}
		}
		T = alpha(T,beta);
	}

	file << "In " << nbiter << "iterations" << endl;
	file << "Best permutation found : ";
	for(auto const & x : bestPk) file << int(x) << " ";
	file << "(" << bestiter << "-th iteration)" << endl;
	file << "With " << bestfPk << " active Sboxes" << endl;
	
}