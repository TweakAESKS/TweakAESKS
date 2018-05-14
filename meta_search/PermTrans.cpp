#include "PermTrans.hpp"

std::pair<bool, PermTrans> merge(PermTrans const & p1, PermTrans const & p2) {
  PermTrans p (p1);
  auto const & transitions = p2.getTransitions();
  for (auto const & my_pair : transitions) {
    if (!p.addTransition(my_pair)) return std::make_pair(false, std::move(p));
  }
  return std::make_pair(true, std::move(p));
}