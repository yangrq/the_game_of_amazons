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
    board bd; //位棋盘
    using distance_matrix = uint8_t[8][8];  //距离矩阵
    using distance_matrix_group = std::array<distance_matrix, 4>; //距离矩阵组
    distance_matrix_group dm1_1{ (uint8_t)-1 }, dm2_1{ (uint8_t)-1 }; //玩家1,2的queen距离矩阵组
    distance_matrix_group dm1_2{ (uint8_t)-1 }, dm2_2{ (uint8_t)-1 }; //玩家1,2的king距离矩阵组
    distance_matrix merged_dm1_1{ (uint8_t)-1 }, merged_dm2_1{ (uint8_t)-1 }; //玩家1,2的queen距离矩阵组合并后的最小queen距离矩阵
    distance_matrix merged_dm1_2{ (uint8_t)-1 }, merged_dm2_2{ (uint8_t)-1 }; //玩家1,2的king距离矩阵组合并后的最小king距离矩阵
  public:
    using player = std::array<board::piece, 4>;  //玩家：包含四个棋子
    using piece = board::piece;
    player player1, player2;
    evaluator(board _bd) noexcept :bd(_bd) {}
    evaluator() noexcept :bd() {};
    ~evaluator() noexcept {};
  private:
    //根据双方最小距离取值的delta函数
    double _territory_determine_delta(uint8_t m, uint8_t n) {
      if (m == 255 && n == 255) return 0;
      if (m == n) return 0.125;
      if (n < m) return 1;
      return -1;
    }
    //生成之后所需的距离矩阵
    void _generate_distance_matrix(const player& p1, const player& p2) {
      int idx = 0;
      for (auto& m : p1)
        _single_queen_min_moves(m, dm1_1[idx++]);
      idx = 0;
      for (auto& m : p2)
        _single_queen_min_moves(m, dm2_1[idx++]);
      idx = 0;
      for (auto& m : p1)
        _single_king_min_moves(m, dm1_2[idx++]);
      idx = 0;
      for (auto& m : p2)
        _single_king_min_moves(m, dm2_2[idx++]);

      _merge_distance_matrix(merged_dm1_1, dm1_1);
      _merge_distance_matrix(merged_dm2_1, dm2_1);
      _merge_distance_matrix(merged_dm1_2, dm1_2);
      _merge_distance_matrix(merged_dm2_2, dm2_2);
    }
    //生成领地（territory）参量
    double _territory_ingredient() {
      auto f_w_t1 = [](double w) { return w; };
      auto f_w_t2 = [](double w) { return w; };
      auto f_w_c1 = [](double w) { return w; };
      auto f_w_c2 = [](double w) { return w; };
      auto [t1, c1, w] = _t1_c1_w();
      auto [t2, c2] = _t2_c2();
      return f_w_t1(w) * t1 + f_w_c1(w) * c1 + f_w_t2(w) * t2 + f_w_c2(w) * c2;
    }
    //将距离矩阵组合并为最小距离矩阵
    void _merge_distance_matrix(distance_matrix& out, const distance_matrix_group& in) {
      for (const auto& mat : in)
        for (int i = 0; i < 8; ++i)
          for (int j = 0; j < 8; ++j)
            out[i][j] = std::min(mat[i][j], out[i][j]);
    }
    //两两计算所有amazon之间几何距离，取得棋子整体分布特征
    double _amazons_distribution(const player& p) {
      double sum = 0;
      for (const auto& amazon1 : p)
        for (const auto& amazon2 : p)
          sum += std::sqrt(std::pow(std::abs((double)amazon1.x() - amazon2.x()), 2) + std::pow(std::abs((double)amazon1.y() - amazon2.y()), 2));
      return sum;
    }
    //计算t1,c1和局势进度特征值w
    std::tuple<double, double, double> _t1_c1_w() {
      double t1 = 0, c1 = 0, w = 0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          t1 += _territory_determine_delta(merged_dm1_1[i][j], merged_dm2_1[i][j]);
          c1 += std::pow(2, -merged_dm1_1[i][j]) - std::pow(2, -merged_dm2_1[i][j]);
          w += std::pow(2, -std::abs(merged_dm1_1[i][j] - merged_dm2_1[i][j]));
        }
      c1 *= 2;
      return std::make_tuple(t1, c1, w);
    }
    //计算t2,c2
    std::tuple<double, double>  _t2_c2() {
      double t2 = 0, c2 = 0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          t2 += _territory_determine_delta(merged_dm1_2[i][j], merged_dm2_2[i][j]);
          c2 += std::min(1.0, std::max(-1.0, (double)(merged_dm2_2[i][j] - merged_dm1_2[i][j]) / 6.0));
        }
      return std::make_tuple(t2, c2);
    }
    //计算相邻的空闲方格数量
    size_t _empty_neighbor_num(uint8_t x, uint8_t y) {
      size_t sum = 0;
      if (x + 1 < 8 && !bd.is_obstacle(x + 1, y)) ++sum;
      if (x - 1 >= 0 && !bd.is_obstacle(x - 1, y)) ++sum;
      if (y + 1 < 8 && !bd.is_obstacle(x, y + 1)) ++sum;
      if (y - 1 >= 0 && !bd.is_obstacle(x, y - 1)) ++sum;
      if (x + 1 < 8 && y + 1 < 8 && !bd.is_obstacle(x + 1, y + 1)) ++sum;
      if (x - 1 >= 0 && y + 1 < 8 && !bd.is_obstacle(x - 1, y + 1)) ++sum;
      if (x + 1 < 8 && y - 1 >= 0 && !bd.is_obstacle(x + 1, y - 1)) ++sum;
      if (x - 1 >= 0 && y - 1 >= 0 && !bd.is_obstacle(x - 1, y - 1)) ++sum;
      return sum;
    }
    //计算单个amazon的queen距离矩阵
    void _single_queen_min_moves(piece from, distance_matrix& distance) {
      std::vector<piece> open;
      std::bitset<64> closed;
      open.push_back(from);
      distance[from.x()][from.y()] = 0;

      while (!open.empty()) {
        piece tmp = open.back();
        uint8_t x = tmp.x();
        uint8_t y = tmp.y();
        uint8_t w = distance[x][y];

        open.pop_back();
        closed[tmp.eigen_value()] = 1;

        for (int i = x + 1; i < 8; ++i) {
          if (bd.is_obstacle(i, y)) break;
          if (!closed[piece::eigen_value(i, y)]) {
            open.emplace_back(i, y);
            distance[i][y] = std::min(distance[i][y], (uint8_t)(w + 1));
          }
        }
        for (int i = x - 1; i >= 0; --i) {
          if (bd.is_obstacle(i, y)) break;
          if (!closed[piece::eigen_value(i, y)]) {
            open.emplace_back(i, y);
            distance[i][y] = std::min(distance[i][y], (uint8_t)(w + 1));
          }
        }
        for (int i = y + 1; i < 8; ++i) {
          if (bd.is_obstacle(x, i)) break;
          if (!closed[piece::eigen_value(x, i)]) {
            open.emplace_back(x, i);
            distance[x][i] = std::min(distance[x][i], (uint8_t)(w + 1));
          }
        }
        for (int i = y - 1; i >= 0; --i) {
          if (bd.is_obstacle(x, i)) break;
          if (!closed[piece::eigen_value(x, i)]) {
            open.emplace_back(x, i);
            distance[x][i] = std::min(distance[x][i], (uint8_t)(w + 1));
          }
        }
        for (int i = x + 1, j = y + 1; i < 8 && j < 8; ++i, ++j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[piece::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        for (int i = x - 1, j = y + 1; i >= 0 && j < 8; --i, ++j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[piece::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        for (int i = x + 1, j = y - 1; i < 8 && j >= 0; ++i, --j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[piece::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        for (int i = x - 1, j = y - 1; i >= 0 && j >= 0; --i, --j) {
          if (bd.is_obstacle(i, j)) break;
          if (!closed[piece::eigen_value(i, j)]) {
            open.emplace_back(i, j);
            distance[i][j] = std::min(distance[i][j], (uint8_t)(w + 1));
          }
        }
        distance[from.x()][from.y()] = -1;
      }
    }
    //计算单个amazon的king距离矩阵
    void _single_king_min_moves(piece from, distance_matrix& distance) {
      std::vector<piece> open;
      std::bitset<64> closed;
      open.push_back(from);
      distance[from.x()][from.y()] = 0;

      while (!open.empty()) {
        piece tmp = open.back();
        uint8_t x = tmp.x();
        uint8_t y = tmp.y();
        uint8_t w = distance[x][y];

        open.pop_back();
        closed[tmp.eigen_value()] = 1;

        if (x + 1 < 8) if (!bd.is_obstacle(x + 1, y) && !closed[piece::eigen_value(x + 1, y)]) {
          open.emplace_back(x + 1, y);
          distance[x + 1][y] = std::min(distance[x + 1][y], (uint8_t)(w + 1));
        }
        if (x - 1 >= 0) if (!bd.is_obstacle(x - 1, y) && !closed[piece::eigen_value(x - 1, y)]) {
          open.emplace_back(x - 1, y);
          distance[x - 1][y] = std::min(distance[x - 1][y], (uint8_t)(w + 1));
        }
        if (y + 1 < 8) if (!bd.is_obstacle(x, y + 1) && !closed[piece::eigen_value(x, y + 1)]) {
          open.emplace_back(x, y + 1);
          distance[x][y + 1] = std::min(distance[x][y + 1], (uint8_t)(w + 1));
        }
        if (y - 1 >= 0) if (!bd.is_obstacle(x, y - 1) && !closed[piece::eigen_value(x, y - 1)]) {
          open.emplace_back(x, y - 1);
          distance[x][y - 1] = std::min(distance[x][y - 1], (uint8_t)(w + 1));
        }
        if (x + 1 < 8 && y + 1 < 8) if (!bd.is_obstacle(x + 1, y + 1) && !closed[piece::eigen_value(x + 1, y + 1)]) {
          open.emplace_back(x + 1, y + 1);
          distance[x + 1][y + 1] = std::min(distance[x + 1][y + 1], (uint8_t)(w + 1));
        }
        if (x - 1 >= 0 && y + 1 < 8) if (!bd.is_obstacle(x - 1, y + 1) && !closed[piece::eigen_value(x - 1, y + 1)]) {
          open.emplace_back(x - 1, y + 1);
          distance[x - 1][y + 1] = std::min(distance[x - 1][y + 1], (uint8_t)(w + 1));
        }
        if (x + 1 < 8 && y - 1 >= 0) if (!bd.is_obstacle(x + 1, y - 1) && !closed[piece::eigen_value(x + 1, y - 1)]) {
          open.emplace_back(x + 1, y - 1);
          distance[x + 1][y - 1] = std::min(distance[x + 1][y - 1], (uint8_t)(w + 1));
        }
        if (x - 1 >= 0 && y - 1 >= 0) if (!bd.is_obstacle(x - 1, y - 1) && !closed[piece::eigen_value(x - 1, y - 1)]) {
          open.emplace_back(x - 1, y - 1);
          distance[x - 1][y - 1] = std::min(distance[x - 1][y - 1], (uint8_t)(w + 1));
        }
        distance[from.x()][from.y()] = -1;
      }
    }
    //计算特定棋子的移动力
    double _amazon_mobility(bool is_player1, size_t amazon_idx) {
      double a = 0.0;
      auto opponent_distance = is_player1 ? merged_dm2_1 : merged_dm1_1;
      auto self_amazon_d1 = is_player1 ? dm1_1 : dm2_1;
      auto self_amazon_d2 = is_player1 ? dm1_2 : dm2_2;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
          if (opponent_distance[i][j] != 255 && self_amazon_d1[amazon_idx][i][j] <= 1)
            a += std::pow(2, -self_amazon_d2[amazon_idx][i][j]) * _empty_neighbor_num(i, j);
      return a;
    }
    //生成移动力（mobility）参量
    double _mobility_ingredient(double w) {
      double m1 = 0, m2 = 0;
      auto f_w_m1 = [w](double m) {return m; };
      auto f_w_m2 = [w](double m) {return m; };
      for (int i = 0; i < 4; ++i)
        m1 += _amazon_mobility(true, (size_t)i);
      for (int i = 0; i < 4; ++i)
        m2 += _amazon_mobility(false, (size_t)i);
      return f_w_m1(m1) - f_w_m2(m2);
    }
  };
}
#endif