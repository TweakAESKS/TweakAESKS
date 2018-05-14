#include "ShortestPath.hpp"

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