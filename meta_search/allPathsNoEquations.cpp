#include <vector>
#include <list>
#include <map>
#include <set>
#include <iostream>
#include <memory>

#include "allPathsNoEquations.hpp"
#include "Subset.hpp"
#include "DicStateKey.hpp"
#include "ShortestPath.hpp"
#include "PermTrans.hpp"


using namespace std;



vector<vector<Subset>> generateTransitionsNoEquation(unsigned const bound_active, unsigned const activeK, DicStateKey dic, vector<uint8_t> const & Ps) {
  unsigned n_columns[5] = {0,1,3,7,15};

  auto const & state_conf = dic.getStateConfs();
  auto const n_conf = state_conf.size();
  auto const & map_key_to_subset = dic.getMapKeyToSubset(activeK);


  uint16_t * tab = new uint16_t [5*5*5*5];

  vector<vector<Subset>> res (n_conf, vector<Subset> (n_conf, Subset()));

  uint16_t const mask[4] = {0x000F, 0x00F0, 0x0F00, 0xF000};

  uint16_t x = 0;
  do { // loop on x1 (2^16 configurations)
    unsigned y0c[4]; // minimal number of active variables in columns of y0 to avoid equations
    for (unsigned c = 0; c < 4; ++c) {
      y0c[c] = 4 - fast_popcnt16(x & mask[c]);
    }
    unsigned y1 = 0;
    for (unsigned i = 0; i < 16; ++i) y1 |= ((x >> i) & 1) << Ps[i];
    uint16_t y1_col = dic.convertStateToStateConf(y1);
    unsigned i_y1_col = dic.getIndexStateConf(y1_col);
    //cout << "\r - generate all transitions (" <<  activeK << "): " << ((100*x)/0x10000) << "%" << " ";
    if (fast_popcnt16(y1_col) < bound_active) {
      unsigned const new_bound = bound_active - fast_popcnt16(y1_col);
      for (auto const & my_pair : map_key_to_subset) { // loop on k, activeK active variables
        uint16_t z = x | my_pair.first;
        unsigned n_tab = 1;
        tab[0] = 0;
        for (unsigned c = 0; c < 4; ++c) {
          auto const i_c = max(y0c[c], 5 - fast_popcnt16(z & mask[c]));
          if ((my_pair.first & mask[c]) == (x & mask[c]))  { // col c of k and x1 are equals
            if (i_c != 5) {
              uint16_t const tmp = n_columns[i_c] << (4*c);
              unsigned const n_tab_tmp = n_tab;
              for (unsigned i = 0; i < n_tab_tmp; ++i) {
                tab[n_tab] = tab[i] | tmp;
                if (fast_popcnt16(tab[n_tab]) < new_bound) ++n_tab;
              }
            }
          }
          else {
            if (i_c != 5) {
              uint16_t const tmp = n_columns[i_c] << (4*c);
              unsigned i = 0;
              while (i < n_tab) {
                tab[i] |= tmp;
                if (fast_popcnt16(tab[i]) < new_bound) ++i;
                else tab[i] = tab[--n_tab];
              }
            }
            else goto lab_break;
          }
        }
        for (unsigned i = 0; i < n_tab; ++i) {
          res[dic.getIndexStateConf(tab[i])][i_y1_col] += my_pair.second;
        }
        lab_break:
        {}
      }
    }
    ++x;
  } while(x != 0);

  delete[] tab;

  for (unsigned i1 = 0; i1 < n_conf; ++i1) {
    for (unsigned i2 = 0; i2 < i1; ++i2) {
      unsigned c = 0;
      while (c < 4 && (((state_conf[i1] & mask[c]) != 0 ) == ((state_conf[i2] & mask[c]) != 0 ))) ++c;
      if (c != 4) continue;
      ////cout << (state_conf[i1] & mask[0]) << " " << (state_conf[i1] & mask[1]) << " " << (state_conf[i1] & mask[2]) << " " << (state_conf[i1] & mask[3]) << endl;
      ////cout << (state_conf[i2] & mask[0]) << " " << (state_conf[i2] & mask[1]) << " " << (state_conf[i2] & mask[2]) << " " << (state_conf[i2] & mask[3]) << endl;
      c = 0;
      while (c < 4 && (state_conf[i1] & mask[c]) <= (state_conf[i2] & mask[c])) ++c;
      if (c == 4) {
        for (unsigned i3 = 0; i3 < n_conf; ++i3) res[i2][i3] += res[i1][i3];
        continue;
      }
      c = 0;
      while (c < 4 && (state_conf[i1] & mask[c]) >= (state_conf[i2] & mask[c])) ++c;
      if (c == 4) {
        for (unsigned i3 = 0; i3 < n_conf; ++i3) res[i1][i3] += res[i2][i3];
      }
    }
  }

  return res;
}

vector<vector<Subset>> generateTransitions(unsigned const bound_active, unsigned const activeK, DicStateKey dic, vector<uint8_t> const & Ps) {
  unsigned n_columns[5] = {0,1,3,7,15};

  auto const & state_conf = dic.getStateConfs();
  auto const n_conf = state_conf.size();
  auto const & map_key_to_subset = dic.getMapKeyToSubset(activeK);


  uint16_t * tab = new uint16_t [5*5*5*5];

  vector<vector<Subset>> res (n_conf, vector<Subset> (n_conf, Subset()));

  uint16_t const mask[4] = {0x000F, 0x00F0, 0x0F00, 0xF000};

  uint16_t x = 0;
  do { // loop on x1 (2^16 configurations)
    unsigned y1 = 0;
    for (unsigned i = 0; i < 16; ++i) y1 |= ((x >> i) & 1) << Ps[i];
    uint16_t y1_col = dic.convertStateToStateConf(y1);
    unsigned i_y1_col = dic.getIndexStateConf(y1_col);
    //cout << "\r - generate all transitions (" <<  activeK << "): " << ((100*x)/0x10000) << "%" << " ";
    if (fast_popcnt16(y1_col) < bound_active) {
      unsigned const new_bound = bound_active - fast_popcnt16(y1_col);
      for (auto const & my_pair : map_key_to_subset) { // loop on k, activeK active variables
        uint16_t z = x | my_pair.first;
        unsigned n_tab = 1;
        tab[0] = 0;
        for (unsigned c = 0; c < 4; ++c) {
          auto const i_c = 5 - fast_popcnt16(z & mask[c]);
          if ((my_pair.first & mask[c]) == (x & mask[c]))  { // col c of k and x1 are equals
            if (i_c != 5) {
              uint16_t const tmp = n_columns[i_c] << (4*c);
              unsigned const n_tab_tmp = n_tab;
              for (unsigned i = 0; i < n_tab_tmp; ++i) tab[n_tab++] = tab[i] | tmp;
            }
          }
          else {
            if (i_c != 5) {
              uint16_t const tmp = n_columns[i_c] << (4*c);
              for (unsigned i = 0; i < n_tab; ++i) tab[i] |= tmp;
            }
            else goto lab_break;
          }
        }
        for (unsigned i = 0; i < n_tab; ++i) {
          if (fast_popcnt16(tab[i]) < new_bound) res[dic.getIndexStateConf(tab[i])][i_y1_col] += my_pair.second;
        }
        lab_break:
        {}
      }
    }
    ++x;
  } while(x != 0);

  delete[] tab;

  for (unsigned i1 = 0; i1 < n_conf; ++i1) {
    for (unsigned i2 = 0; i2 < i1; ++i2) {
      unsigned c = 0;
      while (c < 4 && (((state_conf[i1] & mask[c]) != 0 ) == ((state_conf[i2] & mask[c]) != 0 ))) ++c;
      if (c != 4) continue;
      ////cout << (state_conf[i1] & mask[0]) << " " << (state_conf[i1] & mask[1]) << " " << (state_conf[i1] & mask[2]) << " " << (state_conf[i1] & mask[3]) << endl;
      ////cout << (state_conf[i2] & mask[0]) << " " << (state_conf[i2] & mask[1]) << " " << (state_conf[i2] & mask[2]) << " " << (state_conf[i2] & mask[3]) << endl;
      c = 0;
      while (c < 4 && (state_conf[i1] & mask[c]) <= (state_conf[i2] & mask[c])) ++c;
      if (c == 4) {
        for (unsigned i3 = 0; i3 < n_conf; ++i3) res[i2][i3] += res[i1][i3];
        continue;
      }
      c = 0;
      while (c < 4 && (state_conf[i1] & mask[c]) >= (state_conf[i2] & mask[c])) ++c;
      if (c == 4) {
        for (unsigned i3 = 0; i3 < n_conf; ++i3) res[i1][i3] += res[i2][i3];
      }
    }
  }

  return res;
}




map<unsigned,Subset> updateRight(map<unsigned,Subset> const & bound, vector<Subset> const & transitions, unsigned const bound_active, DicStateKey const & dic) {
  if (bound.empty()) return map<unsigned,Subset> ();
  unsigned pop = 0;
  auto end = bound.end();
  map<unsigned,Subset> new_bound;
  auto const n_conf = dic.nbConf();
  for (unsigned j = 0; j < n_conf; ++j) {
    if (dic.popStateConf(j) != pop) {
      ++pop;
      end = bound.lower_bound(bound_active - pop);
      if (bound.begin() == end) return new_bound;
    }
    auto const & s = transitions[j];
    auto it = bound.begin();
    do {
      if (s.shareElements(it->second)) {
        new_bound[it->first + pop] += j;
        break;
      }
      ++it;
    } while(it != end);
  }
  return new_bound;
}

unsigned minAfterUpdateRight(map<unsigned,Subset> const & bound, vector<Subset> const & transitions, unsigned const bound_active, DicStateKey const & dic) {
  if (bound.empty()) return bound_active;
  unsigned const bound_min = bound.begin()->first;
  unsigned my_min = bound_active;
  unsigned bound_pop = min(17u, my_min - bound_min);
  auto const n_conf = dic.nbConf();
  auto const end = bound.end();
  unsigned j = 0;
  for (unsigned pop = 0; pop < bound_pop; ++pop) {
    Subset s = transitions[j];
    while (++j < n_conf && dic.popStateConf(j) == pop) s += transitions[j];
    auto it = bound.begin();
    while (it != end && it->first + pop < my_min) {
      if (s.shareElements(it->second)) {
        my_min = it->first + pop;
        bound_pop = min(bound_pop, my_min - bound_min);
        break;
      }
      else ++it;
    }
  }
  return my_min;
}

map<unsigned,Subset> updateLeft(map<unsigned,Subset> const & bound, vector<Subset> const & transitions, unsigned const bound_active, DicStateKey const & dic) {
  auto const n_conf = dic.nbConf();
  map<unsigned, Subset> bound2;
  for (auto const & my_pair : bound) {
    my_pair.second.apply([&bound2,&my_pair,&dic](unsigned j){
      unsigned const pop_j =  dic.popStateConf(j);
      if (my_pair.first > pop_j) bound2[my_pair.first - pop_j] += j;
    });
  }
  auto const end = bound2.rend();
  map<unsigned,Subset> new_bound;
  for (unsigned i = 0; i < n_conf; ++i) {
    auto const & s = transitions[i];
    auto it = bound2.rbegin();
    while (it != end) {
      if (s.shareElements(it->second)) {
        new_bound[it->first] += i;
        break;
      }
      ++it;
    }
  }
  return new_bound;
}

bool NoPaths(PermTrans const & p, unsigned const bound_active, unsigned const n_rounds, vector<vector<vector<Subset>>> const & transitions1r_k, DicStateKey const & dic) {
  auto const n_conf = dic.nbConf();

  map<unsigned,Subset> first;
  for (unsigned i = 0; i < n_conf; ++i) {
    if (bound_active > dic.popStateConf(i)) first[dic.popStateConf(i)] += i;
  }

  unsigned const nb_activeK = transitions1r_k.size();
  for (unsigned a = 1; a < nb_activeK; ++a) {
    auto const nb_keys = dic.nbKeys(a);
    for (unsigned k = 0; k < nb_keys; ++k) {
      auto map_min = updateRight(first, transitions1r_k[a][k], bound_active, dic);
      uint16_t real_ki = dic.convertUnsignedToKey(k,a);
      for (unsigned r = 1; r < n_rounds-2; ++r) {
        real_ki = p.getImage(real_ki);
        map_min = updateRight(map_min, transitions1r_k[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
      }
      real_ki = p.getImage(real_ki);
      if (minAfterUpdateRight(map_min, transitions1r_k[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic) < bound_active) return false;
    }
  }

  return true;
}

unsigned NoPathsNoEq(PermTrans const & p, unsigned const bound_active, unsigned const n_rounds, vector<vector<vector<Subset>>> const & transitions1r_k, vector<vector<vector<Subset>>> const & transitions1r_k_noeq, DicStateKey const & dic) {
  auto const n_conf = dic.nbConf();

  map<unsigned,Subset> first;
  for (unsigned i = 0; i < n_conf; ++i) {
    if (bound_active > dic.popStateConf(i)) first[dic.popStateConf(i)] += i;
  }

  unsigned const nb_activeK = transitions1r_k.size();
  for (unsigned a = 1; a < nb_activeK; ++a) {
    auto const nb_keys = dic.nbKeys(a);
    for (unsigned k = 0; k < nb_keys; ++k) {
      {
        auto map_min = updateRight(first, transitions1r_k[a][k], bound_active, dic);
        uint16_t real_ki = dic.convertUnsignedToKey(k,a);
        for (unsigned r = 1; r < n_rounds-2; ++r) {
          real_ki = p.getImage(real_ki);
          map_min = updateRight(map_min, transitions1r_k_noeq[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
        }
        real_ki = p.getImage(real_ki);
        unsigned min = minAfterUpdateRight(map_min, transitions1r_k_noeq[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
        if(min < bound_active) return min;
        // if (minAfterUpdateRight(map_min, transitions1r_k_noeq[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic) < bound_active) return false;
      }
      auto map_min_backup = updateRight(first, transitions1r_k[a][k], bound_active, dic);
      uint16_t real_ki_backup = dic.convertUnsignedToKey(k,a);
      for (unsigned r_eq = 1; r_eq < n_rounds-1; ++r_eq) {
        auto map_min = map_min_backup;
        auto real_ki = real_ki_backup;
        for (unsigned r = 1; r < n_rounds-1; ++r) {
          real_ki = p.getImage(real_ki);
          if (r < n_rounds-2) {
            if (r == r_eq) map_min = updateRight(map_min, transitions1r_k[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
            else map_min = updateRight(map_min, transitions1r_k_noeq[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
          }
          else {
            if (r == r_eq) {
              unsigned min = minAfterUpdateRight(map_min, transitions1r_k[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
              if(min < bound_active) return min;
              // if (minAfterUpdateRight(map_min, transitions1r_k[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic) < bound_active) return false;
            }
            else {
              unsigned min = minAfterUpdateRight(map_min, transitions1r_k_noeq[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic);
              if(min < bound_active) return min;
              // if (minAfterUpdateRight(map_min, transitions1r_k_noeq[a][dic.convertKeyToUnsigned(real_ki)], bound_active, dic) < bound_active) return false;
            }
          }
        }
      }
    }
  }

  return bound_active;
}


bool compatibleMaps(map<unsigned, Subset> const & map_min, map<unsigned, Subset> const & map_max) {
    auto it_min = map_min.begin();
    auto const end_min = map_min.end();
    auto it_max = map_max.begin();
    auto const end_max = map_max.end();
    Subset s;
    while (it_max != end_max) {
      while (it_min != end_min && it_min->first < it_max->first) {
        s += it_min->second;
        ++it_min;
      }
      if (s.shareElements(it_max->second)) return false;
      ++it_max;
    }
    return true;
}

// template <typename T1, typename T2>
// void printMap(map<T1, T2> const & my_map) {
//   for (auto const & my_pair : my_map) {
//     //cout << my_pair.first << ": " << my_pair.second << endl;
//   }
//   getchar();
// }


void printTransitions(vector<uint16_t> const & vec) {
  for (unsigned l = 0; l < 4; ++l) {
    for (auto k : vec) {
      for (unsigned c = 0; c < 4; ++c) {
        //cout << ((k >> (4*l + c)) & 1) << " ";
      }
      //cout << "  -->  ";
    }
    //cout << endl;
  }
  //cout << endl;
}


// #include "CPsolve.hpp"

void searchPermRandomMax(DicStateKey const & dic, unsigned const activeK, unsigned const bound_active, vector<uint8_t> const & Ps, unsigned const n_rounds) {
  auto const & state_conf = dic.getStateConfs();
  auto const n_conf = state_conf.size();
  //cout << "start search" << endl;
  /* Generate transitions */
  vector<vector<vector<Subset>>> transitions1r (activeK+1);
  vector<vector<vector<Subset>>> transitions1r_ji (activeK+1);
  //vector<vector<vector<Subset>>> transitions1r_ij (activeK+1);
  for (unsigned a = 1; a <= activeK; ++a) {
    auto transitions1r_a = generateTransitions(bound_active, a, dic, Ps);
    //cout << "\r - generate all transitions (" <<  a << "): 100%" << endl;
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
    //cout << "\r - generate all transitions (" <<  a << "): 100%" << endl;
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

  map<unsigned, Subset> map_x0;
  for (unsigned i = 0; i < n_conf; ++i) {
    map_x0[dic.popStateConf(i)] += i;
  }


  {
    unsigned const a = 4;
    unsigned const nb_keys = dic.nbKeys(a);
    vector<unsigned> all_keys (nb_keys);
    vector<vector<unsigned>> path_min (nb_keys, vector<unsigned> (nb_keys));
    for (unsigned k0 = 0; k0 < nb_keys; ++k0) {
      //cout << "\r" << k0 << flush;
      auto map_x1 = updateRight(map_x0, transitions1r_ji[a][k0], bound_active, dic);
      auto map_x1_noeq = updateRight(map_x0, transitions1r_ji_noeq[a][k0], bound_active, dic);
      for (unsigned k1 = 0; k1 < nb_keys; ++k1) {
        path_min[k0][k1] = minAfterUpdateRight(map_x1, transitions1r_ji_noeq[a][k1], bound_active, dic);
        path_min[k0][k1] = min(path_min[k0][k1], minAfterUpdateRight(map_x1_noeq, transitions1r_ji[a][k1], bound_active, dic));
      }
    }
    //cout << endl;
    std::default_random_engine generator;
    for (unsigned i = 0; i < nb_keys; ++i) all_keys[i] = i;
    unsigned cpt = 0;
    for (;;) {
      //cout << "\r" << ++cpt << flush;
      PermTrans p;
      while (!p.isSet()) {
        std::uniform_int_distribution<int> distribution(0,nb_keys-1);
        unsigned k0 = distribution(generator);
        while (fast_popcnt16(p.getImage(dic.convertUnsignedToKey(k0,a))) == a) k0 = (k0 + 1) % nb_keys;
        uint16_t real_k0 = dic.convertUnsignedToKey(k0, a);
        vector<unsigned> tmp;
        unsigned my_min = 0;
        for (unsigned k1 = 0; k1 < nb_keys; ++k1) {
          if (!p.transitionIsCompatible(real_k0, dic.convertUnsignedToKey(k1, a))) continue;
          if (path_min[k0][k1] < my_min) continue;
          if (path_min[k0][k1] > my_min) {
            my_min = path_min[k0][k1];
            tmp.clear();
          }
          tmp.emplace_back(k1);
        }
        distribution = std::uniform_int_distribution<int> (0, tmp.size()-1);
        p.addTransition(real_k0, dic.convertUnsignedToKey(tmp[distribution(generator)], a));
      }
      unsigned best_bound = 20;
      if (NoPathsNoEq(p, best_bound, n_rounds, transitions1r_ji, transitions1r_ji_noeq, dic) >= best_bound) {
        // if (shortestPathAndState(p.convertToPermutation(), Ps, 6, best_bound)) {
          ++best_bound;
          //cout << p << endl << endl;
        // }
      }
    }
  }
}

void searchIntuition(DicStateKey const & dic, unsigned const activeK, unsigned const bound_active, vector<uint8_t> const & Ps, unsigned const n_rounds) {
  auto const & state_conf = dic.getStateConfs();
  auto const n_conf = state_conf.size();
  //cout << "start search" << endl;
  /* Generate transitions */
  vector<vector<vector<Subset>>> transitions1r (activeK+1);
  vector<vector<vector<Subset>>> transitions1r_ji (activeK+1);
  //vector<vector<vector<Subset>>> transitions1r_ij (activeK+1);
  for (unsigned a = 1; a <= activeK; ++a) {
    auto transitions1r_a = generateTransitions(bound_active, a, dic, Ps);
    //cout << "\r - generate all transitions (" <<  a << "): 100%" << endl;
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
    //cout << "\r - generate all transitions (" <<  a << "): 100%" << endl;
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

  {
    uint16_t columns[4] = {0xF, 0xF0, 0xF00, 0xF000};
    uint16_t sr[4];
    for (unsigned i = 0; i < 4; ++i) {
      uint16_t tmp = 0;
      for (unsigned j = 0; j < 16; ++j) tmp |= ((columns[i] >> Ps[j]) & 1) << j;
      sr[i] = tmp;
    }
    vector<uint16_t> best_sr;
    {
      uint16_t k = 0;
      do {
        if (fast_popcnt16(k) == 4) {
          if (fast_popcnt16(k & sr[0]) == 1 && fast_popcnt16(k & sr[1]) == 1 && fast_popcnt16(k & sr[2]) == 1 && fast_popcnt16(k & sr[3]) == 1) {
              best_sr.emplace_back(k);
          }
        }
        ++k;
      } while(k != 0);
    }
    vector<uint16_t> best_col;
    {
      uint16_t k = 0;
      do {
        if (fast_popcnt16(k) == 4) {
          if (fast_popcnt16(k & columns[0]) == 1 && fast_popcnt16(k & columns[1]) == 1 && fast_popcnt16(k & columns[2]) == 1 && fast_popcnt16(k & columns[3]) == 1) {
            best_col.emplace_back(k);
          }
        }
        ++k;
      } while(k != 0);
    }

    //cout << "best: " << best_sr.size() << " - " << best_col.size() << endl;
    vector<PermTrans> all_perm (1, PermTrans());
    for (unsigned i = 0; i < 4; ++i) {
      vector<PermTrans> all_perm_tmp;
      for (auto const & p : all_perm) {
        for (auto k : best_sr) {
          auto p1 (p);
          if (p1.addTransition(columns[i], k)) all_perm_tmp.emplace_back(move(p1));
        }
      }
      all_perm = move(all_perm_tmp);
    }
    for (unsigned i = 0; i < 4; ++i) {
      vector<PermTrans> all_perm_tmp;
      for (auto const & p : all_perm) {
        for (auto k : best_col) {
          auto p1 (p);
          if (p1.addTransition(k,sr[i])) all_perm_tmp.emplace_back(move(p1));
        }
      }
      all_perm = move(all_perm_tmp);
    }
    //cout << "all: " << all_perm.size() << endl;
    unsigned best_bound = 5;
    uint64_t cpt = 0;
    for (auto const & p : all_perm) {
      //cout << "\r" << ++cpt << flush;
      if (NoPathsNoEq(p, best_bound, n_rounds, transitions1r_ji, transitions1r_ji_noeq, dic) >= best_bound) {
        do {
          ++best_bound;
        } while(NoPathsNoEq(p, best_bound, n_rounds, transitions1r_ji, transitions1r_ji_noeq, dic) >= best_bound);
        //cout << endl << "minimal path: " << best_bound - 1 << endl << p << endl;
      }
    }
    getchar();
  }
}

PermTrans vecToPermTrans(std::vector<uint8_t> const & P){
  PermTrans p;
  for(uint8_t x = 0; x < 16; x++){
    uint16_t i0 = (1 << x);
    uint16_t i1 = (1 << P[x]);
    p.addTransition(i0,i1);
  }
  return p;
}

// int main(int argc, char const *argv[]) {
//   vector<uint8_t> Ps = {0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3}; // ShiftRows
//   DicStateKey dic;
//   searchPermRandomMax(dic, 7, 22, Ps, 6);
//   return 0;
// }
