#ifndef EVALUATOR_H
#define EVALUATOR_H
#include <vector>
#include <cstdint>
#include <cmath>
#include <bitset>
#include <array>
#include <algorithm>
#include <utility>
#include <initializer_list>
#include <iomanip>
#include "bitmap.hpp"
#include "board.hpp"

namespace yrq {

  class evaluator {
    board bd;
    using distance_matrix = uint8_t[8][8];
  public:
    using teil = board::teil;
    evaluator(board _bd) noexcept :bd(_bd) {}
    evaluator() noexcept :bd() {};
    ~evaluator() noexcept {};
  private:
    double _territory_determine_delta(uint8_t m, uint8_t n) {
      if (m == 255 && n == 255) return 0;
      if (m == n) return 0.125;
      if (n < m) return 1;
      return -1;
    }
    double _territory_ingredient(const std::array<teil, 4>& p1, const std::array<teil, 4>& p2) {
      auto f_w_t1 = [](double w) {
        return 0.1 + 0.9 * std::pow(2, -w);
      };
      auto f_w_t2 = [](double w) { return std::pow(w, 0.2); };
      auto f_w_c1 = [](double w) { return 0.3 * std::pow(std::log(w + 1) / (w + 5), 0.5); };
      auto f_w_t4 = [](double w) { return 0.08 * std::atan(w); };
      auto [t1, c1, w] = _t1_c1_w(p1, p2);
      auto [t2, c2] = _t2_c2(p1, p2);
    }
    std::tuple<double, double, double> _t1_c1_w(const std::array<teil, 4>& p1, const std::array<teil, 4>& p2) {
      distance_matrix dm1, dm2;
      memset(dm1, -1, sizeof(dm1));
      memset(dm2, -1, sizeof(dm2));
      for (auto& m : p1)
        _single_queen_min_moves(m, dm1);
      for (auto& m : p2)
        _single_queen_min_moves(m, dm2);
      double t1 = 0, c1 = 0, w = 0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          t1 += _territory_determine_delta(dm1[i][j], dm2[i][j]);
          c1 += std::pow(2, -dm1[i][j]) - std::pow(2, -dm2[i][j]);
          w += std::pow(2, -std::abs(dm1[i][j] - dm2[i][j]));
        }
      c1 *= 2;
      return std::make_tuple(t1, c1, w);
    }
    std::pair<double, double>  _t2_c2(const std::array<teil, 4>& p1, const std::array<teil, 4>& p2) {
      distance_matrix dm1, dm2;
      memset(dm1, -1, sizeof(dm1));
      memset(dm2, -1, sizeof(dm2));
      for (auto& m : p1)
        _single_king_min_moves(m, dm1);
      for (auto& m : p2)
        _single_king_min_moves(m, dm2);
      double t2 = 0, c2 = 0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          t2 += _territory_determine_delta(dm1[i][j], dm2[i][j]);
          c2 += std::min(1.0, std::max(-1.0, (double)(dm2[i][j] - dm1[i][j]) / 6.0));
        }
      return std::make_pair(t2, c2);
    }
    void _single_queen_min_moves(teil from, distance_matrix& distance) {
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
    void _single_king_min_moves(teil from, distance_matrix& distance) {
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