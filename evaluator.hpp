#ifndef EVALUATOR_H
#define EVALUATOR_H
#include <vector>
#include <cstdint>
#include <cmath>
#include <bitset>
#include <array>
#include <algorithm>
#include <functional>
#include <utility>
#include <initializer_list>
#include <iomanip>
#include "bitmap.hpp"
#include "board.hpp"
#include "utility.hpp"

namespace yrq {
  class evaluator {
    board bd; //位棋盘
    using distance_matrix = uint8_t[8][8];  //距离矩阵
    using distance_matrix_group = std::array<distance_matrix, 4>; //距离矩阵组
    using player = std::array<board::piece, 4>;  //玩家：包含四个棋子
    using piece = board::piece;
    std::array<distance_matrix_group, 2> dm_1;    //玩家1,2的queen距离矩阵组
    std::array<distance_matrix_group, 2> dm_2;    //玩家1,2的king距离矩阵组
    std::array<distance_matrix, 2> merged_dm_1;   //玩家1,2的queen距离矩阵组合并后的最小queen距离矩阵
    std::array<distance_matrix, 2> merged_dm_2;   //玩家1,2的king距离矩阵组合并后的最小king距离矩阵
  public:
    std::array<player, 2> players;
    double w = 0;
    evaluator(board _bd) noexcept :bd(_bd) {
      for (auto& dmg : dm_1)
        for (auto& dm : dmg)
          memset(dm, (uint8_t)(-1), 64);
      for (auto& dmg : dm_2)
        for (auto& dm : dmg)
          memset(dm, (uint8_t)(-1), 64);

      memset(merged_dm_1[0], (uint8_t)(-1), 64);
      memset(merged_dm_1[1], (uint8_t)(-1), 64);
      memset(merged_dm_2[0], (uint8_t)(-1), 64);
      memset(merged_dm_2[1], (uint8_t)(-1), 64);
    }
    evaluator() noexcept {
      for (auto& dmg : dm_1)
        for (auto& dm : dmg)
          memset(dm, (uint8_t)(-1), 64);
      for (auto& dmg : dm_2)
        for (auto& dm : dmg)
          memset(dm, (uint8_t)(-1), 64);

      memset(merged_dm_1[0], (uint8_t)(-1), 64);
      memset(merged_dm_1[1], (uint8_t)(-1), 64);
      memset(merged_dm_2[0], (uint8_t)(-1), 64);
      memset(merged_dm_2[1], (uint8_t)(-1), 64);
    };
    ~evaluator() noexcept {};
  public:
    //输出距离矩阵
    void _debug_printf_distance_matrix(distance_matrix dm) {
      for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
          std::cout << std::setw(4) << (unsigned)dm[j][i];
        std::cout << std::endl;
      }
    }
    //根据双方最小距离取值的delta函数
    double _territory_determine_delta(uint8_t m, uint8_t n) {
      if (m == 255 && n == 255) return 0;
      if (m == n) return 0.125;
      if (m < n) return 1;
      return -1;
    }
    //生成之后所需的距离矩阵
    void _generate_distance_matrix() {
      int idx = 0;
      for (auto& m : players[0])
        _single_queen_min_moves(m, dm_1[0][idx++]);
      idx = 0;
      for (auto& m : players[1])
        _single_queen_min_moves(m, dm_1[1][idx++]);
      idx = 0;
      for (auto& m : players[0])
        _single_king_min_moves(m, dm_2[0][idx++]);
      idx = 0;
      for (auto& m : players[1])
        _single_king_min_moves(m, dm_2[1][idx++]);

      _merge_distance_matrix(merged_dm_1[0], dm_1[0]);
      _merge_distance_matrix(merged_dm_1[1], dm_1[1]);
      _merge_distance_matrix(merged_dm_2[0], dm_2[0]);
      _merge_distance_matrix(merged_dm_2[1], dm_2[1]);
#ifdef _DEBUG
      std::cout << "merged_dm_1 player0:" << std::endl;
      _debug_printf_distance_matrix(merged_dm_1[0]);
      std::cout << "merged_dm_1 player1:" << std::endl;
      _debug_printf_distance_matrix(merged_dm_1[1]);
      std::cout << "merged_dm_2 player0:" << std::endl;
      _debug_printf_distance_matrix(merged_dm_2[0]);
      std::cout << "merged_dm_2 player1:" << std::endl;
      _debug_printf_distance_matrix(merged_dm_2[1]);
#endif
    }
    //生成领地（territory）参量
    double _territory_ingredient() {
      auto [t1, c1, w] = _t1_c1_w();
      auto [t2, c2] = _t2_c2();
      // f(t1,w) = [ 0.75 * 1.1 ^ (-w) + 0.25 ] * t1
      auto f_w_t1 = [=](double v) { return (0.75 * std::pow(1.1, -w) + 0.25) * v; };
      // f(t2,w) = [ 0.08 * sqrt( max { w-1 , 0 } ) ] * t2
      auto f_w_t2 = [=](double v) { return (0.08 * std::sqrt(w - 1 > 0 ? w - 1 : 0))* v; };
      // f(c1,w) = [ 1 - f_t1(w) - f_t2(w) ] * [ 0.6 * 1.1 ^ (-w) + 0.4 ] * c1
      auto f_w_c1 = [=](double v) { return (1 - 0.75 * std::pow(1.1, -w) - 0.25 - 0.08 * std::sqrt(w - 1 > 0 ? w - 1 : 0))* (0.4 + 0.6 * std::pow(1.1, -w))* v; };
      // f(c2,w) = [ 1 - f_t1(w) - f_t2(w) ] * [ 0.6 - 0.6 * 1.1 ^ (-w) ] * c2
      auto f_w_c2 = [=](double v) { return (1 - 0.75 * std::pow(1.1, -w) - 0.25 - 0.08 * std::sqrt(w - 1 > 0 ? w - 1 : 0))* (1 - (0.4 + 0.6 * std::pow(1.1, -w)))* v; };
#ifdef _DEBUG
      std::printf("t1:%lf c1:%lf t2:%lf c2:%lf w:%lf\n", t1, c1, t2, c2, w);
      std::printf("f(t1,w):%lf f(c1,w):%lf f(t2,w):%lf f(c2,w):%lf\n", f_w_t1(t1), f_w_c1(c1), f_w_t2(t2), f_w_c2(c2));
      std::printf("t:%lf\n", f_w_t1(t1) + f_w_c1(c1) + f_w_t2(t2) + f_w_c2(c2));
      emit_key_value("t1", t1);
      emit_key_value("c1", c1);
      emit_key_value("t2", t2);
      emit_key_value("c2", c2);
      emit_key_value("w", w, true);
      emit_key_value("f(t1,w)", f_w_t1(t1));
      emit_key_value("f(c1,w)", f_w_c1(c1));
      emit_key_value("f(t2,w)", f_w_t2(t2));
      emit_key_value("f(c2,w)", f_w_c2(c2));
      emit_key_value("t", f_w_t1(t1) + f_w_c1(c1) + f_w_t2(t2) + f_w_c2(c2), true);
#endif
      return f_w_t1(t1) + f_w_c1(c1) + f_w_t2(t2) + f_w_c2(c2);
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
#ifdef _DEBUG
      std::printf("pieces distribution [(%d,%d) (%d,%d) (%d,%d) (%d,%d)]:%lf\n",
        p[0].x(), p[0].y(),
        p[1].x(), p[1].y(),
        p[2].x(), p[2].y(),
        p[3].x(), p[3].y(),
        std::sqrt(sum / 10.0) - 1.5);
      emit_key_value("d", std::sqrt(sum / 10.0) - 1.5, true);
#endif
      return std::sqrt(sum / 10.0) - 3.0;
    }
    //计算t1,c1和局势进度特征值w
    std::tuple<double, double, double> _t1_c1_w() {
      double t1 = 0, c1 = 0, w = 0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          if (bd.is_obstacle(i, j)) continue;
          t1 += _territory_determine_delta(merged_dm_1[0][i][j], merged_dm_1[1][i][j]);
          c1 += std::pow(2.0, -merged_dm_1[0][i][j]) - std::pow(2.0, -merged_dm_1[1][i][j]);
          w += std::pow(2.0, -std::abs(merged_dm_1[0][i][j] - merged_dm_1[1][i][j]));
        }
      c1 *= 2;
      this->w = w;
      return std::make_tuple(t1, c1, w);
    }
    //计算t2,c2
    std::tuple<double, double>  _t2_c2() {
      double t2 = 0, c2 = 0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          if (bd.is_obstacle(i, j)) continue;
          t2 += _territory_determine_delta(merged_dm_2[0][i][j], merged_dm_2[1][i][j]);
          c2 += std::min(1.0, std::max(-1.0, (double)(merged_dm_2[1][i][j] - merged_dm_2[0][i][j]) / 6.0));
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
    double _amazon_mobility(size_t player_idx, size_t amazon_idx) {
      double a = 0.0;
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
          if (merged_dm_1[1 - player_idx][i][j] != 255 && dm_1[player_idx][amazon_idx][i][j] <= 1)
            a += std::pow(2.0, -dm_2[player_idx][amazon_idx][i][j]) * _empty_neighbor_num(i, j);
#ifdef _DEBUG
      std::printf("the mobility of (%d,%d):%lf\n", players[player_idx][amazon_idx].x(), players[player_idx][amazon_idx].y(), a);
#endif
      return a;
    }
    //生成移动力（mobility）参量
    double _mobility_ingredient() {
      double m1 = 0, m2 = 0;
      auto f_w_m1 = [this](double m) {return 2 * (w < 10 ? 10 : w) * std::pow(1.2, -m); };
      auto f_w_m2 = [this, &f_w_m1](double m) {return f_w_m1(m); };
      for (int i = 0; i < 4; ++i)
        m1 += f_w_m1(_amazon_mobility(0, (size_t)i));
      for (int i = 0; i < 4; ++i)
        m2 += f_w_m2(_amazon_mobility(1, (size_t)i));
#ifdef _DEBUG
      std::printf("m1:%lf m2:%lf f(m1,w):%lf f(m2,w):%lf\n", m1, m2, f_w_m1(m1), f_w_m2(m2));
      std::printf("m:%lf\n", f_w_m1(m1) - f_w_m2(m2));
      emit_key_value("m1", m1);
      emit_key_value("m2", m2);
      emit_key_value("f(m1,w)", f_w_m1(m1));
      emit_key_value("f(m2,w)", f_w_m2(m2));
      emit_key_value("m", f_w_m2(m2) - f_w_m1(m1));
#endif
      return f_w_m2(m2) - f_w_m1(m1);
    }
    //扁平化二维距离数组
    class _distance_flat_wrapper {
      std::array<distance_matrix_group, 2>* p;
    public:
      _distance_flat_wrapper(std::array<distance_matrix_group, 2>& v) :p(&v) {};
      ~_distance_flat_wrapper() = default;
      distance_matrix& operator[](size_t idx) {
        return (*p)[idx / 4][idx % 4];
      }
    };
    _distance_flat_wrapper _flat_dm_1{ dm_1 };
    //只有一项满足
    std::tuple<bool, size_t> _only_one_satisfy(std::function<bool(size_t)> condition) {
      for (int i = 0; i < 8; ++i)
        if (condition(i)) {
          bool satisfy = true;
          for (int j = 0; j < 8; ++j) {
            if (i == j) continue;
            if (condition(j)) goto end;
          }
          return std::make_tuple(true, i);
        }
    end:
      return std::make_tuple(false, 0);
    }
    //计算amazon的独占区域集合
    std::array<std::array<size_t, 4>, 2> _amazon_exclusive_access_num() {
      auto is_reachable_closure = [this](uint8_t x, uint8_t y) {
        return [x, y, this](size_t piece_idx) {
          return _flat_dm_1[piece_idx][x][y] != 255;
        };
      };
      std::array<std::array<size_t, 4>, 2> exclusive_access_num = { 0 };
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
          auto [found, result] = _only_one_satisfy(is_reachable_closure(i, j));
          if (found) ++exclusive_access_num[result / 4][result % 4];
        }
#ifdef _DEBUG
      std::printf("player0: the number of exclusively accessible squares\n");
      for (int i = 0; i < 4; ++i)
        std::printf("%llu ", exclusive_access_num[0][i]);
      printf("\n");
      std::printf("player1: the number of exclusively accessible squares\n");
      for (int i = 0; i < 4; ++i)
        std::printf("%llu ", exclusive_access_num[1][i]);
      printf("\n");
#endif
      return exclusive_access_num;
    }
    //生成守卫（guard）参量
    double _guard_ingredient() {
      long long sum = 0;
      auto res = _amazon_exclusive_access_num();
      for (auto v : res[0]) sum += v;
      for (auto v : res[1]) sum -= v;
#ifdef _DEBUG
      std::printf("n:%llu g:%lf\n", sum, 0.2 * (sum > 0 ? std::pow(1.1, sum) : -std::pow(1.1, -sum)));
      emit_key_value("g(n)", 0.2 * (sum > 0 ? std::pow(1.1, sum) : -std::pow(1.1, -sum)), true);
#endif
      return 0.2 * (sum > 0 ? std::pow(1.1, sum) : -std::pow(1.1, -sum));
    }
  };
}
#endif 