#ifndef HPP_SHORTESTPATH
#define HPP_SHORTESTPATH

#include <cstdint>

#include "DicStateKey.hpp"

class ShortestPath {
public:
  ShortestPath (unsigned n_rounds = 0, unsigned n_conf = 0) : dec0 (n_conf), dec1 (n_conf*(n_rounds-1)) {
    shortest_path = new uint8_t [dec0*dec1];
  };

  ShortestPath (std::vector<std::vector<Subset>> const & transitions1r, unsigned bound_active, unsigned n_rounds, DicStateKey const & dic) {
    auto const & state_conf = dic.getStateConfs();
    dec0 = state_conf.size();
    dec1 = dec0*(n_rounds-1);
    shortest_path = new uint8_t [dec0*dec1];
    for (unsigned i = 0; i < dec0*dec1; ++i) shortest_path[i] = bound_active;
    for (unsigned i0 = 0; i0 < dec0; ++i0) {
      for (unsigned i1 = 0; i1 < dec0; ++i1) {
        if (!transitions1r[i0][i1].empty()) (*this)(0,i0,i1) = fast_popcnt16(state_conf[i0]) + fast_popcnt16(state_conf[i1]);
      }
    }
    for (unsigned r = 1; r < n_rounds - 1; ++r) {
      for (unsigned i0 = 0; i0 < dec0; ++i0) {
        for (unsigned i1 = 0; i1 < dec0; ++i1) {
          auto const pathi0i1 = (*this)(r-1,i0,i1);
          if (pathi0i1 < bound_active) {
            auto const pop_i1 = fast_popcnt16(state_conf[i1]);
            for (unsigned i2 = 0; i2 < dec0; ++i2) {
              uint8_t const tmp = pathi0i1 - pop_i1 + (*this)(0,i1,i2);
              (*this)(r,i0,i2) = std::min((*this)(r,i0,i2), tmp);
            }
          }
        }
      }
    }
  }

  ShortestPath (ShortestPath const & s) : dec0(s.dec0), dec1(s.dec1) {
    auto const bound = dec0*dec1;
    shortest_path = new uint8_t [bound];
    for (unsigned i = 0; i < bound; ++i) shortest_path[i] = s.shortest_path[i];
  }
  ShortestPath (ShortestPath && s) : dec0(s.dec0), dec1(s.dec1), shortest_path(s.shortest_path) {
    s.shortest_path = nullptr;
    s.dec0 = 0;
    s.dec1 = 0;
  }
  ShortestPath & operator=(ShortestPath const & s) {
    dec0 = s.dec0;
    dec1 = s.dec1;
    delete[] shortest_path;
    auto const bound = dec0*dec1;
    shortest_path = new uint8_t [bound];
    for (unsigned i = 0; i < bound; ++i) shortest_path[i] = s.shortest_path[i];
    return *this;
  }
  ShortestPath & operator=(ShortestPath && s) {
    dec0 = s.dec0;
    dec1 = s.dec1;
    delete[] shortest_path;
    shortest_path = s.shortest_path;
    s.shortest_path = nullptr;
    s.dec0 = 0;
    s.dec1 = 0;
    return *this;
  }

  ~ShortestPath() {delete[] shortest_path;};

  unsigned size() const {return dec0;};
  unsigned deep() const {return dec1/dec0;};

  uint8_t operator()(unsigned r, unsigned i0, unsigned i1) const {return shortest_path[dec1*i1 + dec0*r + i0];};
  uint8_t & operator()(unsigned r, unsigned i0, unsigned i1) {return shortest_path[dec1*i1 + dec0*r + i0];};

private:
  unsigned dec0;
  unsigned dec1;
  uint8_t * shortest_path;
};

uint64_t countPaths(ShortestPath const & s, unsigned i0, unsigned ir, unsigned bound_active, unsigned n_rounds, DicStateKey const & dic);

std::vector<std::vector<Subset>> transitionsFromPaths(ShortestPath const & s, unsigned i0, unsigned ir, unsigned bound_active, unsigned n_rounds, DicStateKey const & dic, std::vector<std::vector<Subset>> const & transitions1r);

bool transitions1rIsValid(unsigned const i0, unsigned const i1, ShortestPath const & s, unsigned const bound_active, unsigned const n_rounds, DicStateKey const & dic);

bool transitions2rIsValid(unsigned const i0, unsigned const i1, unsigned const i2, ShortestPath const & s, unsigned const bound_active, unsigned const n_rounds, DicStateKey const & dic);



#endif
