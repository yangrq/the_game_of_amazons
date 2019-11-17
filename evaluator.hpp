#ifndef EVALUATOR_H
#define EVALUATOR_H
#include <vector>
#include <cstdint>
#include <bitset>
#include <algorithm>
#include <iomanip>
#include "bitmap.hpp"
#include "board.hpp"

namespace yrq {

  class evaluator {
    board bd;
    uint8_t distance[8][8];
  public:
    using teil = board::teil;
    evaluator(board _bd) noexcept :bd(_bd) { memset(distance, -1, 64); }
    evaluator() noexcept :bd() { memset(distance, -1, 64); };
    ~evaluator() noexcept {};
    void distance_output() noexcept {
      for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
          std::cout << std::setw(4) << (unsigned)distance[j][i];
        std::cout << std::endl;
      }
    }
    void queen_min_moves(teil from) {
      std::vector<teil> open;
      std::bitset<64> closed;
      open.push_back(from);
      distance[from.x()][from.y()] = 0;

      while (!open.empty()) {
        teil tmp = open.back();
        uint8_t x = tmp.x();
        uint8_t y = tmp.y();
        uint8_t w = distance[x][y];

        open.pop_back();
        closed[tmp.eigen_value()] = 1;

        for (int i = x + 1; i < 8; ++i) {
          if (bd.is_obstacle(i, y)) break;
          if (!closed[teil::eigen_value(i, y)]) {
            open.emplace_back(i, y);
            distance[i][y] = std::min(distance[i][y], (uint8_t)(w + 1));
          }
        }
        for (int i = x - 1; i >= 0; --i) {
          if (bd.is_obstacle(i, y)) break;
          if (!closed[teil::eigen_value(i, y)]) {
            open.emplace_back(i, y);
            distance[i][y] = std::min(distance[i][y], (uint8_t)(w + 1));
          }
        }
        for (int i = y + 1; i < 8; ++i) {
          if (bd.is_obstacle(x, i)) break;
          if (!closed[teil::eigen_value(x, i)]) {
            open.emplace_back(x, i);
            distance[x][i] = std::min(distance[x][i], (uint8_t)(w + 1));
          }
        }
        for (int i = y - 1; i >= 0; --i) {
          if (bd.is_obstacle(x, i)) break;
          if (!closed[teil::eigen_value(x, i)]) {
            open.emplace_back(x, i);
            distance[x][i] = std::min(distance[x][i], (uint8_t)(w + 1));
          }
        }
        for (int i = x + 1, j = y + 1; i < 8 && j < 8; ++i, ++j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[teil::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        for (int i = x - 1, j = y + 1; i >= 0 && j < 8; --i, ++j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[teil::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        for (int i = x + 1, j = y - 1; i < 8 && j >= 0; ++i, --j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[teil::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        for (int i = x - 1, j = y - 1; i >= 0 && j >= 0; --i, --j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[teil::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        distance[from.x()][from.y()] = -1;
      }

    }
    void king_min_moves(teil from) {
      std::vector<teil> open;
      std::bitset<64> closed;
      open.push_back(from);
      distance[from.x()][from.y()] = 0;

      while (!open.empty()) {
        teil tmp = open.back();
        uint8_t x = tmp.x();
        uint8_t y = tmp.y();
        uint8_t w = distance[x][y];

        open.pop_back();
        closed[tmp.eigen_value()] = 1;

        if (x + 1 < 8) if (!bd.is_obstacle(x + 1, y) && !closed[teil::eigen_value(x + 1, y)]) {
          open.emplace_back(x + 1, y);
          distance[x + 1][y] = std::min(distance[x + 1][y], (uint8_t)(w + 1));
        }
        if (x - 1 >= 0) if (!bd.is_obstacle(x - 1, y) && !closed[teil::eigen_value(x - 1, y)]) {
          open.emplace_back(x - 1, y);
          distance[x - 1][y] = std::min(distance[x - 1][y], (uint8_t)(w + 1));
        }
        if (y + 1 < 8) if (!bd.is_obstacle(x, y + 1) && !closed[teil::eigen_value(x, y + 1)]) {
          open.emplace_back(x, y + 1);
          distance[x][y + 1] = std::min(distance[x][y + 1], (uint8_t)(w + 1));
        }
        if (y - 1 >= 0) if (!bd.is_obstacle(x, y - 1) && !closed[teil::eigen_value(x, y - 1)]) {
          open.emplace_back(x, y - 1);
          distance[x][y - 1] = std::min(distance[x][y - 1], (uint8_t)(w + 1));
        }
        if (x + 1 < 8 && y + 1 < 8) if (!bd.is_obstacle(x + 1, y + 1) && !closed[teil::eigen_value(x + 1, y + 1)]) {
          open.emplace_back(x + 1, y + 1);
          distance[x + 1][y + 1] = std::min(distance[x + 1][y + 1], (uint8_t)(w + 1));
        }
        if (x - 1 >= 0 && y + 1 < 8) if (!bd.is_obstacle(x - 1, y + 1) && !closed[teil::eigen_value(x - 1, y + 1)]) {
          open.emplace_back(x - 1, y + 1);
          distance[x - 1][y + 1] = std::min(distance[x - 1][y + 1], (uint8_t)(w + 1));
        }
        if (x + 1 < 8 && y - 1 >= 0) if (!bd.is_obstacle(x + 1, y - 1) && !closed[teil::eigen_value(x + 1, y - 1)]) {
          open.emplace_back(x + 1, y - 1);
          distance[x + 1][y - 1] = std::min(distance[x + 1][y - 1], (uint8_t)(w + 1));
        }
        if (x - 1 >= 0 && y - 1 >= 0) if (!bd.is_obstacle(x - 1, y - 1) && !closed[teil::eigen_value(x - 1, y - 1)]) {
          open.emplace_back(x - 1, y - 1);
          distance[x - 1][y - 1] = std::min(distance[x - 1][y - 1], (uint8_t)(w + 1));
        }
        distance[from.x()][from.y()] = -1;
      }

    }
  };
}
#endif