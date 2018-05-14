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

uint64_t countPaths(ShortestPath const & s, unsigned i0, unsigned ir, unsigned bound_active, unsigned n_rounds, DicStateKey const & dic) {
  if (n_rounds == 2) return 1;
  auto const & state_conf = dic.getStateConfs();
  uint64_t cpt = 0;
  for (unsigned i = 0; i < s.size(); ++i) {
    auto const pop_real_i = fast_popcnt16(state_conf[i]);
    if (s(0,i0,i) + s(n_rounds-3,i,ir) < bound_active + pop_real_i) cpt += countPaths(s, i, ir, bound_active + pop_real_i - s(0,i0,i), n_rounds - 1, dic);
  }
  return cpt;
}

std::vector<std::vector<Subset>> transitionsFromPaths(ShortestPath const & s, unsigned i0, unsigned ir, unsigned bound_active, unsigned n_rounds, DicStateKey const & dic, std::vector<std::vector<Subset>> const & transitions1r) {
  if (n_rounds < 2) return std::vector<std::vector<Subset>> ();
  if (n_rounds == 2) return std::vector<std::vector<Subset>> (1, std::vector<Subset> (1, transitions1r[i0][ir]));
  std::vector<std::vector<Subset>> res;
  auto const & state_conf = dic.getStateConfs();
  for (unsigned i = 0; i < s.size(); ++i) {
    auto const pop_real_i = fast_popcnt16(state_conf[i]);
    if (s(0,i,ir) + s(n_rounds-3,i0,i) < bound_active + pop_real_i) {
      auto vec = transitionsFromPaths(s, i0, i, bound_active + pop_real_i - s(0,i,ir), n_rounds-1, dic, transitions1r);
      for (auto & v : vec) {
        v.emplace_back(transitions1r[i][ir]);
        res.emplace_back(move(v));
      }
    }
  }
  return res;
}

bool transitions1rIsValid(unsigned const i0, unsigned const i1, ShortestPath const & s, unsigned const bound_active, unsigned const n_rounds, DicStateKey const & dic) {
  if (s(0,i0,i1) >= bound_active) return false;
  auto const & state_conf = dic.getStateConfs();
  auto const n_conf = state_conf.size();
  auto const pop_i0 = fast_popcnt16(state_conf[i0]);
  auto const pop_i1 = fast_popcnt16(state_conf[i1]);
  for (unsigned i = 0; i < n_conf; ++i) if (pop_i0 + s(n_rounds-3,i1,i) < bound_active) return true;
  for (unsigned i = 0; i < n_conf; ++i) if (pop_i1 + s(n_rounds-3,i,i0) < bound_active) return true;
  for (unsigned r = 0; r <= n_rounds-4; ++r) {
    auto min_i = s(r,0,i0);
    for (unsigned i = 1; i < n_conf; ++i) min_i = std::min(min_i, s(r,i,i0));
    for (unsigned j = 0; j < n_conf; ++j) if (min_i + s(n_rounds-r-4,i1,j) < bound_active) return true;
  }
  return false;
}

bool transitions2rIsValid(unsigned const i0, unsigned const i1, unsigned const i2, ShortestPath const & s, unsigned const bound_active, unsigned const n_rounds, DicStateKey const & dic) {
  auto const & state_conf = dic.getStateConfs();
  auto const n_conf = state_conf.size();
  auto const a_2r = s(0,i0,i1) + s(0,i1,i2) - fast_popcnt16(state_conf[i1]);
  if (a_2r >= bound_active) return false;
  auto const pop_i2 = fast_popcnt16(state_conf[i2]);
  for (unsigned i = 0; i < n_conf; ++i) if (a_2r - pop_i2 + s(n_rounds-4,i2,i) < bound_active) return true;
  auto const pop_i0 = fast_popcnt16(state_conf[i0]);
  for (unsigned i = 0; i < n_conf; ++i) if (a_2r - pop_i0 + s(n_rounds-4,i,i0) < bound_active) return true;
  for (unsigned r = 0; n_rounds >= 5 + r; ++r) {
    for (unsigned i = 0; i < n_conf; ++i) {
      auto const a = a_2r - pop_i0 + s(r,i,i0);
      if (a < bound_active) {
        for (unsigned j = 0; j < n_conf; ++j) if (a - pop_i2 + s(n_rounds-(5+r),i2,j) < bound_active) return true;
      }
    }
  }
  return false;
}



#endif
