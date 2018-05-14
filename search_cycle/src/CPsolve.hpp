#ifndef H_CPSOLVE
#define H_CPSOLVE

#include <vector>
#include <string>
#include <cstdint>


std::string exec(const char* cmd);
//Exec cmd and grab the stdout

void createPermFile(std::vector<uint8_t> const & Pk, std::vector<uint8_t> const & Ps, unsigned const nround, unsigned const bound);


bool shortestPathAndState(std::vector<uint8_t> const & Pk, std::vector<uint8_t> const & Ps, unsigned const Nround, unsigned const bound);

#endif
