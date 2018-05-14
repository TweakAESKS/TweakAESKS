#include "DicStateKey.hpp"

std::vector<unsigned> init_fast_popcnt256() {
  std::vector<unsigned> res (256);
  res[0] = 0;
  for (int i = 1; i < 256; ++i) res[i] = (i & 1) + res[i / 2];
  return res;
};

unsigned fast_popcnt16(uint16_t x) {
  static std::vector<unsigned> wordbits = init_fast_popcnt256();
  uint8_t * p = (uint8_t *) &x;
  return wordbits[p[0]] + wordbits[p[1]];
};