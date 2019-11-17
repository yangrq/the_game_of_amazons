#ifndef BOARD_H
#define BOARD_H
#include <cstdint>
#include "bitmap.hpp"

namespace yrq {
  class board {
    bitmap amazon;
    bitmap arrow;
  public:
    class teil {
      friend teil;
      uint8_t idx;
    public:
      teil() noexcept :idx(0) {}
      teil(uint8_t x, uint8_t y, bool is_obstacle = false) noexcept :idx(x << 5 | y << 2 | (uint8_t)is_obstacle) {}
      teil(uint8_t _idx) noexcept : idx(_idx) {}
      uint8_t x() { return idx >> 5; }
      uint8_t y() { return (idx & 0x1C) >> 2; }
      bool is_obstacle() { return idx & 1; }
      bool operator==(const teil& v) {
        return v.idx == idx;
      }
      uint8_t eigen_value() {
        return idx >> 2;
      }
      static uint8_t eigen_value(uint8_t x, uint8_t y) {
        return x << 3 | y;
      }
    };
  public:
    board() noexcept :amazon(), arrow() {};
    board(bitmap amazon, bitmap _arrow = 0) noexcept :amazon(amazon), arrow(_arrow) {};
    ~board() noexcept {};
    bitmap& get_queen_map() { return amazon; }
    bitmap& get_arrow_map() { return arrow; }
    bitmap accessible(int y, int x) {
      bitmap r1, r2, r3, r4, obstacle(amazon | arrow);
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
    bool is_obstacle(teil t) { return (amazon | arrow)[t.y()][t.x()]; }
    bool is_obstacle(uint8_t x, uint8_t y) { return (amazon | arrow)[y][x]; }
  private:
  };
}
#endif