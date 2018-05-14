#ifndef HPP_PermTrans
#define HPP_PermTrans

#include <vector>
#include <cstdint>
#include <random>

#include "DicStateKey.hpp"

class PermTrans {
public:
  PermTrans () : transitions (1, std::make_pair(0xFFFF, 0xFFFF)) {};
  PermTrans (uint16_t i0, uint16_t i1) : transitions ({std::make_pair(i0, i1),std::make_pair(~i0, ~i1)}) {};
  PermTrans (std::pair<uint16_t, uint16_t> const & p) : transitions ({p, std::make_pair(~p.first, ~p.second)}) {};

  bool addTransition(uint16_t const i0, uint16_t const i1) {
    auto const n_transitions = transitions.size();
    for (unsigned i = 0; i != n_transitions; ++i) {
      auto const tmp0 = transitions[i].first & i0;
      auto const tmp1 = transitions[i].second & i1;
      if (fast_popcnt16(tmp0) != fast_popcnt16(tmp1)) return false;
      if (tmp0 != 0 && tmp0 != transitions[i].first) {
        transitions[i].first &= ~i0;
        transitions[i].second &= ~i1;
        transitions.emplace_back(tmp0, tmp1);
      }
    }
    return true;
  }
  bool addTransition(std::pair<uint16_t, uint16_t> const & p) {return addTransition(p.first, p.second);};

  bool transitionIsCompatible(uint16_t const i0, uint16_t const i1) {
    for (auto const & my_pair : transitions) {
      if (fast_popcnt16(my_pair.first & i0) != fast_popcnt16(my_pair.second & i1)) return false;
    }
    return true;
  }

  bool isSet() const {return (transitions.size() == 16);};

  uint16_t getPreImage(uint16_t const x) const {
    uint16_t res = 0;
    for (auto const & my_pair : transitions) {
      if ((my_pair.second & x) != 0) res |= my_pair.first;
    }
    return res;
  }

  uint16_t getImage(uint16_t const x) const {
    uint16_t res = 0;
    for (auto const & my_pair : transitions) {
      if ((my_pair.first & x) != 0) res |= my_pair.second;
    }
    return res;
  }

  void composeWithTransition(int i, int j) {
    std::swap(transitions[i].first, transitions[j].first);
  }

  std::vector<std::pair<uint16_t, uint16_t>> const & getTransitions() const {return transitions;};

  std::vector<uint8_t> convertToPermutation() const {
    std::vector<uint8_t> res (16);
    std::map<uint16_t, unsigned> my_map;
    for (unsigned i = 0; i < 16; ++i) my_map[1 << i] = i;
    for (auto const & my_pair : transitions) {
      res[my_map[my_pair.first]] = my_map[my_pair.second];
    }
    return res;
  }

  friend std::ostream & operator<<( std::ostream &flux, PermTrans const & p) {
    for (auto const & my_pair : p.transitions) {
      auto x = my_pair.first;
      flux << "{";
      unsigned i = 0;
      while (x != 0) {
        if (x & 1) flux << i << ",";
        x >>= 1;
        ++i;
      }
      flux << "} -> {";
      x = my_pair.second;
      i = 0;
      while (x != 0) {
        if (x & 1) flux << i << ",";
        x >>= 1;
        ++i;
      }
      flux << "}" << std::endl;
    }
    return flux;
  }

  friend bool operator<(PermTrans const & p1, PermTrans const & p2) {return p1.transitions < p2.transitions;};

private:
  std::vector<std::pair<uint16_t, uint16_t>> transitions;
};

std::pair<bool, PermTrans> merge(PermTrans const & p1, PermTrans const & p2) {
  PermTrans p (p1);
  auto const & transitions = p2.getTransitions();
  for (auto const & my_pair : transitions) {
    if (!p.addTransition(my_pair)) return std::make_pair(false, std::move(p));
  }
  return std::make_pair(true, std::move(p));
}

#endif
