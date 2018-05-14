#ifndef HPP_DICSTATEKEY
#define HPP_DICSTATEKEY

#include <vector>
#include <map>

#include "Subset.hpp"

std::vector<unsigned> init_fast_popcnt256();

unsigned fast_popcnt16(uint16_t x);


template <typename itT>
bool next_ord_binary (itT begin, itT end) { // similar to next_permutation but generates all binomial combinations
  //first combination: 0.....01...1
  //last combination: 1...10.....0
  int nb_zero_last = 0;
  while (begin != end && *(--end) == 0) ++nb_zero_last;
  if (begin == end && *end == 0) return false;
  int nb_one_last = 1;
  while (begin != end && *(--end) == 1) ++nb_one_last;
  bool non_last = !(begin == end && *end == 1);
  if (non_last) {
    *end = 1;
    ++end;
    ++nb_zero_last;
    --nb_one_last;
  }
  while (nb_zero_last > 0) {
    *end = 0;
    ++end;
    --nb_zero_last;
  }
  while (nb_one_last > 0) {
    *end = 1;
    ++end;
    --nb_one_last;
  }
  return non_last;
};

class DicStateKey {
public:
  DicStateKey () :  key_to_subset (17), subset_to_key (17), unsigned_to_key (17), state_conf (1,0), all_keys (17), key_to_unsigned (1 << 16) {
    for (unsigned activeK = 0; activeK <= 16; ++activeK) {
      unsigned cpt = 0;
      std::vector<int> K (16,1);
      for (unsigned i = 0; i < 16-activeK; ++i) K[i] = 0;
      do { // loop on k, activeK active variables
        uint16_t k = K[0];
        for (unsigned i = 1; i < 16; ++i) k |= K[i] << i;
        Subset k_subset (cpt);
        key_to_subset[activeK].emplace(k, k_subset);
        subset_to_key[activeK].emplace(k_subset, k);
        unsigned_to_key[activeK].emplace_back(k);
        key_to_unsigned[k] = cpt;
        ++cpt;
      } while(next_ord_binary(K.begin(), K.end()));
    }
    {
      unsigned n_columns[5] = {0,1,3,7,15}; // 0000, 0001, 0011, 0111, 1111
      for (unsigned c = 0; c < 4; ++c) {
        std::vector<uint16_t> tmp;
        for (auto conf : n_columns) {
          auto new_conf = conf << 4*c;
          for (auto actual_conf : state_conf) tmp.emplace_back(actual_conf | new_conf);
        }
        state_conf = move(tmp);
      }
      auto const size_state_conf = state_conf.size();
      sort(state_conf.begin(), state_conf.end(), [](uint16_t i, uint16_t j){return fast_popcnt16(i) < fast_popcnt16(j);});
      for (unsigned i = 0; i < size_state_conf; ++i) reverse_state_conf.emplace(state_conf[i],i);
      for (unsigned i = 0; i < size_state_conf; ++i) pop_state_conf.emplace_back(fast_popcnt16(state_conf[i]));
    }
    {
      unsigned n_columns[5] = {0,1,3,7,15}; // 0000, 0001, 0011, 0111, 1111
      uint16_t x = 0;
      do { // loop on x1 (2^16 configurations)
        uint16_t x_col = 0;
        for (unsigned c = 0; c < 4; ++c) {
          uint16_t mask = 0xF << 4*c;
          x_col |= n_columns[fast_popcnt16(x & mask)] << 4*c;
        }
        state_to_state_conf.emplace_back(x_col);
        ++x;
      } while (x != 0);
    }
    {
      for (unsigned activeK = 0; activeK <= 16; ++activeK) {
        Subset s;
        for (auto const & my_pair : subset_to_key[activeK]) s += my_pair.first;
        all_keys[activeK] = std::move(s);
      }
    }
  };
  uint16_t getKey(Subset const & s, unsigned const activeK) const {return subset_to_key[activeK].at(s);};
  Subset const & getSubset(uint16_t const k) const {return key_to_subset[fast_popcnt16(k)].at(k);};
  std::map<uint16_t, Subset> const & getMapKeyToSubset(unsigned const activeK) const {return key_to_subset[activeK];};
  unsigned nbKeys(unsigned const activeK) const {return unsigned_to_key[activeK].size();};
  std::vector<std::map<uint16_t, Subset>> const & getMapKeyToSubset() const {return key_to_subset;};
  std::map<Subset, uint16_t> const & getMapSubsetToKey(unsigned const activeK) const {return subset_to_key[activeK];};
  std::vector<uint16_t> convertSubsetToKeys(Subset const & s, unsigned const activeK) const {
    std::vector<uint16_t> res;
    s.apply([&res,this,activeK](unsigned x){res.emplace_back(unsigned_to_key[activeK][x]);});
    return res;
  };
  uint16_t convertUnsignedToKey(unsigned const k, unsigned const activeK) const {return unsigned_to_key[activeK][k];};
  unsigned convertKeyToUnsigned(uint16_t const k) const {return key_to_unsigned[k];};
  std::vector<uint16_t> const & getStateConfs() const {return state_conf;};
  unsigned getIndexStateConf(uint16_t const state) const {return reverse_state_conf.at(state);};
  uint16_t convertStateToStateConf(uint16_t const state) const {return state_to_state_conf[state];};
  Subset const & getAllKeys(unsigned const activeK) const {return all_keys[activeK];};
  unsigned popStateConf(unsigned const state) const {return pop_state_conf[state];};
  unsigned nbConf() const {return state_conf.size();};

private:
  std::vector<std::map<uint16_t, Subset>> key_to_subset;
  std::vector<std::map<Subset, uint16_t>> subset_to_key;
  std::vector<std::vector<uint16_t>> unsigned_to_key;
  std::vector<uint16_t> state_conf;
  std::map<uint16_t, unsigned> reverse_state_conf;
  std::vector<uint16_t> state_to_state_conf;
  std::vector<Subset> all_keys;
  std::vector<unsigned> key_to_unsigned;
  std::vector<unsigned> pop_state_conf;
};




#endif
