#ifndef H_CPSOLVE
#define H_CPSOLVE

#include <vector>
#include <string>
#include <cstdint>

typedef unsigned int uint;
std::string exec(const char* cmd);
//Exec cmd and grab the stdout

void createPermFile(std::vector<uint8_t> const & Pk, std::vector<uint8_t> const & Ps, uint const nround);

uint shortestPath(std::vector<uint8_t> const & Pk, std::vector<uint8_t> const & Ps, uint const Nround);

#endif
