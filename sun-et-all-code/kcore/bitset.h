#ifndef _BITSET_H
#define _BITSET_H

#include <vector>

#define BITS_PER_CHAR 8
#define BITS_PER_INT (sizeof(int) * BITS_PER_CHAR)

class Bitset {
public:
  Bitset() : vec(4){};
  inline void set(int index) {
    set(index / BITS_PER_INT, index & (BITS_PER_INT - 1));
  }
  inline void reset(int index) {
    reset(index / BITS_PER_INT, index & (BITS_PER_INT - 1));
  }
  inline bool get(int index) const {
    if (index / BITS_PER_INT >= (int)vec.size()) {
      return false;
    }
    return ((vec[index / BITS_PER_INT] >> (index & (BITS_PER_INT - 1))) & 1) ==
           1;
  }
  int nextAvailable() {
    int i;
    for (i = 0; i < (int)vec.size() * BITS_PER_INT; i++) {
      if (!get(i)) {
        return i;
      }
    }
    vec.push_back(0);
    return i;
  }

  void mix(const Bitset &other) {
    if (other.vec.size() > vec.size()) {
      vec.resize(other.vec.size());
    }
    for (int i = 0; i < vec.size(); i++) {
      vec[i] |= other.vec[i];
    }
  }

private:
  inline void set(int index, int delta) {
    if (index >= (int)vec.size()) {
      vec.resize(index + 1);
    }
    vec[index] |= (1 << delta);
  }
  inline void reset(int index, int delta) {
    if (index >= (int)vec.size()) {
      vec.resize(index);
    }
    vec[index] &= (~(1 << delta));
  }
  std::vector<int> vec;
};

#endif