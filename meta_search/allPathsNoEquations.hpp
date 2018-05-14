#ifndef H_ALLPATH
#define H_ALLPATH
#include <vector>
#include <list>
#include <map>
#include <set>
#include <iostream>
#include <memory>

#include "Subset.hpp"
#include "DicStateKey.hpp"
#include "ShortestPath.hpp"
#include "PermTrans.hpp"

std::vector<std::vector<Subset>> generateTransitionsNoEquation(unsigned const bound_active, unsigned const activeK, DicStateKey dic, std::vector<uint8_t> const & Ps);
std::vector<std::vector<Subset>> generateTransitions(unsigned const bound_active, unsigned const activeK, DicStateKey dic, std::vector<uint8_t> const & Ps);
std::map<unsigned,Subset> updateRight(std::map<unsigned,Subset> const & bound, std::vector<Subset> const & transitions, unsigned const bound_active, DicStateKey const & dic);
unsigned minAfterUpdateRight(std::map<unsigned,Subset> const & bound, std::vector<Subset> const & transitions, unsigned const bound_active, DicStateKey const & dic);
std::map<unsigned,Subset> updateLeft(std::map<unsigned,Subset> const & bound, std::vector<Subset> const & transitions, unsigned const bound_active, DicStateKey const & dic);
bool NoPaths(PermTrans const & p, unsigned const bound_active, unsigned const n_rounds, std::vector<std::vector<std::vector<Subset>>> const & transitions1r_k, DicStateKey const & dic);
unsigned NoPathsNoEq(PermTrans const & p, unsigned const bound_active, unsigned const n_rounds, std::vector<std::vector<std::vector<Subset>>> const & transitions1r_k, std::vector<std::vector<std::vector<Subset>>> const & transitions1r_k_noeq, DicStateKey const & dic);
bool compatibleMaps(std::map<unsigned, Subset> const & map_min, std::map<unsigned, Subset> const & map_max);

template <typename T1, typename T2>
void printMap(std::map<T1, T2> const & my_map) {
  for (auto const & my_pair : my_map) {
    std::cout << my_pair.first << ": " << my_pair.second << std::endl;
  }
  getchar();
};

void printTransitions(std::vector<uint16_t> const & vec);
void searchPermRandomMax(DicStateKey const & dic, unsigned const activeK, unsigned const bound_active, std::vector<uint8_t> const & Ps, unsigned const n_rounds);
void searchIntuition(DicStateKey const & dic, unsigned const activeK, unsigned const bound_active, std::vector<uint8_t> const & Ps, unsigned const n_rounds);
PermTrans vecToPermTrans(std::vector<uint8_t> const & P);

#endif