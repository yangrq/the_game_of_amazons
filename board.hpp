#ifndef BOARD_H
#define BOARD_H
#ifndef ALL_IN_ONE
#include <cstdint>
#endif
#include "bitmap.hpp"

namespace yrq {
  class board {
    bitmap amazon;
    bitmap arrow;
  public:
    class piece {
      uint8_t idx;
    public:
      piece() noexcept :idx(0) {}
      piece(uint8_t x, uint8_t y, bool is_obstacle = false) noexcept :idx(x << 5 | y << 2 | (uint8_t)is_obstacle) {}
      piece(uint8_t _idx) noexcept : idx(_idx) {}
      uint8_t x() const { return idx >> 5; }
      uint8_t y() const { return (idx & 0x1C) >> 2; }
      bool is_obstacle() { return idx & 1; }
      bool operator==(const piece& v) const {
        return v.idx == idx;
      }
      uint8_t eigen_value() const {
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
    const bitmap& get_queen_map() const { return amazon; }
    bitmap& get_arrow_map() { return arrow; }
    const bitmap& get_arrow_map() const { return arrow; }
    bool is_obstacle(piece t) { return (amazon | arrow)[t.y()][t.x()]; }
    bool is_obstacle(uint8_t x, uint8_t y) { return (amazon | arrow)[y][x]; }
  private:
  };
}
#endif