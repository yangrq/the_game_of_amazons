#ifndef BOARD_H
#define BOARD_H
#include <cstdint>
#include <exception>
#include "bitmap.hpp"

namespace yrq {
  class board {
    bitmap queen;
    bitmap arrow;
  public:
    board() noexcept :queen(), arrow() {};
    board(bitmap _queen, bitmap _arrow) noexcept :queen(_queen), arrow(_arrow) {};
    ~board() noexcept {};
    bitmap& get_queen_map() { return queen; }
    bitmap& get_arrow_map() { return arrow; }
    bitmap accessible(int y, int x) {
      bitmap r1, r2, r3, r4, obstacle(queen | arrow);
      for (int i = x + 1; i < 8; ++i)
        if (!obstacle[i][y]) r1[i][y] = 1;
        else break;
      for (int i = x - 1; i >= 0; --i)
        if (!obstacle[i][y]) r2[i][y] = 1;
        else break;
      for (int i = y + 1; i < 8; ++i)
        if (!obstacle[x][i]) r3[x][i] = 1;
        else break;
      for (int i = y - 1; i >= 0; --i)
        if (!obstacle[x][i]) r4[x][i] = 1;
        else break;
      for (int i = x + 1, j = y + 1; i < 8 && j < 8; ++i, ++j)
        if (!obstacle[i][j]) r1[i][j] = 1;
        else break;
      for (int i = x + 1, j = y - 1; i < 8 && j >= 0; ++i, --j)
        if (!obstacle[i][j]) r2[i][j] = 1;
        else break;
      for (int i = x - 1, j = y + 1; i >= 0 && j < 8; --i, ++j)
        if (!obstacle[i][j]) r3[i][j] = 1;
        else break;
      for (int i = x - 1, j = y - 1; i >= 0 && j >= 0; --i, --j)
        if (!obstacle[i][j]) r4[i][j] = 1;
        else break;
      return r1 | r2 | r3 | r4;
    }
    bitmap accessible_raw(int offset) {
      int x = offset % 8; int y = offset / 8;
      return accessible(y, x);
    }
  private:
  };
}
#endif