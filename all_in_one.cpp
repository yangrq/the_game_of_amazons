# 1 "main.cpp"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "main.cpp"
# 1 "board.hpp" 1





# 1 "bitmap.hpp" 1
# 9 "bitmap.hpp"
namespace yrq {
  class bitmap {
    class strip {
      friend bitmap;
      bitmap* pbitmap;
      size_t idx;
      class reference {
        friend strip;
        strip* pstrip;
        size_t idx;
      public:
        reference() noexcept :idx(0), pstrip(nullptr) {}
        reference(strip& _strip, size_t _idx) :idx(_idx), pstrip(&_strip) {}
        ~reference() noexcept {}
        reference& operator=(bool x) {
          pstrip->set(idx, x);
          return *this;
        }
        reference& operator=(const reference& _bitref) noexcept {
          pstrip->set(idx, (bool)_bitref);
          return *this;
        }
        reference& flip(size_t _pos) noexcept {
          pstrip->set(_pos, !pstrip->get(_pos));
          return *this;
        }
        operator bool() const noexcept {
          return pstrip->get(idx);
        }
        bool operator~() const noexcept {
          return !pstrip->get(idx);
        }
      };
    private:
      void _set(size_t _pos, bool _val) {
        pbitmap->_set(idx * 8 + _pos, _val);
      }
      bool _get(size_t _pos) {
        return pbitmap->_get(idx * 8 + _pos);
      }
    public:
      strip() noexcept :idx(0), pbitmap(nullptr) {}
      strip(bitmap& _bitmap, size_t _idx) : idx(_idx), pbitmap(&_bitmap) {}
      ~strip() noexcept {}
      reference operator[](size_t _pos) {
        pbitmap->is_valid(_pos, 0, 7);
        return reference(*this, _pos);
      }
      void set(size_t _pos, bool _val) noexcept {
        pbitmap->is_valid(_pos, 0, 7);
        _set(_pos, _val);
      }
      bool get(size_t _pos) noexcept {
        pbitmap->is_valid(_pos, 0, 7);
        return _get(_pos);
      }
    };
  private:
    void _set(size_t _pos, bool _val) noexcept {
      if (_val) raw |= 1ULL << _pos;
      else raw &= ~(1ULL << _pos);
    }
    bool _get(size_t _pos) noexcept {
      return (raw & 1ULL << _pos) ? true : false;
    }
    void is_valid(size_t x, size_t lower_bound, size_t upper_bound) {



    }
    std::uint64_t raw;
  public:
    bitmap() noexcept :raw(0ULL) {};
    bitmap(const std::uint64_t v) :raw(v) {};
    ~bitmap() noexcept {}
    strip operator[](size_t _pos) {
      is_valid(_pos, 0, 7);
      return strip(*this, _pos);
    }
    operator uint64_t() {
      return raw;
    }
    void output() noexcept {
      for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
          std::cout << ((*this)[i][j] ? "1 " : "0 ");
        std::cout << std::endl;
      }
    }
    void output_xy() noexcept {
      for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
          std::cout << ((*this)[j][i] ? "1 " : "0 ");
        std::cout << std::endl;
      }
    }
    void set(size_t _pos, bool _val) noexcept {
      is_valid(_pos, 0, 63);
      return _set(_pos, _val);
    }
    bool get(size_t _pos) noexcept {
      is_valid(_pos, 0, 63);
      return get(_pos);
    }
    bitmap operator|(const bitmap& v) {
      return bitmap(raw | v.raw);
    }
  };
}
# 7 "board.hpp" 2

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
# 2 "main.cpp" 2
# 1 "evaluator.hpp" 1
# 17 "evaluator.hpp"
# 1 "utility.hpp" 1
# 18 "evaluator.hpp" 2

namespace yrq {
  class evaluator {
  public:
    using distance_matrix = uint8_t[8][8]; //距离矩阵
    using distance_matrix_group = std::array<distance_matrix, 4>; //距离矩阵组
    using player = std::array<board::piece, 4>; //玩家：包含四个棋子
    using piece = board::piece;
  private:
    board bd; //位棋盘
    std::array<distance_matrix_group, 2> dm_1; //玩家1,2的queen距离矩阵组
    std::array<distance_matrix_group, 2> dm_2; //玩家1,2的king距离矩阵组
    std::array<distance_matrix, 2> merged_dm_1; //玩家1,2的queen距离矩阵组合并后的最小queen距离矩阵
    std::array<distance_matrix, 2> merged_dm_2; //玩家1,2的king距离矩阵组合并后的最小king距离矩阵
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
    double evaluate() {
      double r = 0;
      _generate_distance_matrix();
      r += _territory_ingredient();
      r += _mobility_ingredient();
      r += _guard_ingredient();
      r += _distribution_ingredient();



      return r;
    }
  private:
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
# 134 "evaluator.hpp"
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
      return std::sqrt(sum / 10.0) - 1.5;
    }
    //生成分布（distribution）参量
    double _distribution_ingredient() {
      double d0 = _amazons_distribution(players[0]);
      double d1 = _amazons_distribution(players[1]);



      return w / 20.0 * (d0 - d1);
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
    std::tuple<double, double> _t2_c2() {
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
      return exclusive_access_num;
    }
    //生成守卫（guard）参量
    double _guard_ingredient() {
      long long sum = 0;
      auto res = _amazon_exclusive_access_num();
      for (auto v : res[0]) sum += v;
      for (auto v : res[1]) sum -= v;



      return 0.2 * (sum > 0 ? std::pow(1.1, sum) : -std::pow(1.1, -sum));
    }
  };
}
# 3 "main.cpp" 2

# 1 "searcher.hpp" 1





namespace yrq {
  struct move_action {
    board bd;
    board::piece from;
    board::piece to;
    board::piece arrow;
    move_action() {};
    move_action(board _bd,
      board::piece _from,
      board::piece _to,
      board::piece _arrow) {
      bd = _bd;
      from = _from;
      to = _to;
      arrow = _arrow;
    }
  };
  move_action make_move(board bd, board::piece from, board::piece to, board::piece arrow) {
    bd.get_arrow_map()[from.x()][from.y()] = 0;
    bd.get_arrow_map()[to.x()][to.y()] = 1;
    bd.get_arrow_map()[arrow.x()][arrow.y()] = 1;
    return { bd,from,to,arrow };
  }
  class greedy_searcher {
    std::array<evaluator::player, 2> players;
    board bd;
  public:
    greedy_searcher() {
      bd.get_arrow_map() = bitmap(0x2400810000810024);
      players[0] = {
        board::piece{2, 0},
        board::piece{0, 2},
        board::piece{5, 0},
        board::piece{7, 2} };
      players[1] = {
        board::piece{0, 5},
        board::piece{2, 7},
        board::piece{5, 7},
        board::piece{7, 5} };
    }
    ~greedy_searcher() {};
    void set_board(board bd) {
      this->bd = bd;
    }
    void set_players(std::array<evaluator::player, 2> players) {
      this->players = players;
    }
    move_action search_and_select(const std::vector<move_action>& mvs = {}) {
      _apply_all_move_actions(mvs);
      move_action mva;
      double current_evaluation = -1000.0;
      for (int i = 0; i < 4; ++i) {
        auto res = _possible_moves(players[0][i]);
        for (int j = 0; j < res.size(); ++j) {
          evaluator ev(res[j].bd);
          ev.players[0] = players[0];
          ev.players[1] = players[1];
          double evaluation = ev.evaluate();
          if (evaluation > current_evaluation) {
            current_evaluation = evaluation;
            mva = res[j];
          }
        }
      }
      mva.bd.get_arrow_map().output_xy();
      return mva;
    }
  private:
    void _apply_move_action(const move_action& mv) {
      for (int i = 0; i < 4; ++i)
        if (mv.from == players[0][i]) {
          players[0][i] = mv.to;
          bd = mv.bd;
        }
    }
    void _apply_all_move_actions(const std::vector<move_action>& mvs) {
      for (const auto& mv : mvs)
        _apply_move_action(mv);
    }
    std::vector<move_action> _possible_arrow_placements(std::vector<move_action>& mvs) {
      std::vector<move_action> res;
      for (auto& mv : mvs) {
        uint8_t x = mv.to.x();
        uint8_t y = mv.to.y();
        board nbd = mv.bd;
        for (int i = x + 1; i < 8; ++i) {
          if (mv.bd.is_obstacle((uint8_t)i, y)) break;
          nbd.get_arrow_map()[(uint8_t)i][y] = 1;
          res.push_back({ nbd, mv.from, mv.to, board::piece{(uint8_t)i,y} });
          nbd.get_arrow_map()[i][y] = 0;
        }
        for (int i = x - 1; i >= 0; --i) {
          if (mv.bd.is_obstacle((uint8_t)i, y)) break;
          nbd.get_arrow_map()[(uint8_t)i][y] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{(uint8_t)i,y} });
          nbd.get_arrow_map()[i][y] = 0;
        }
        for (int i = y + 1; i < 8; ++i) {
          if (bd.is_obstacle(x, (uint8_t)i)) break;
          nbd.get_arrow_map()[x][(uint8_t)i] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{x,(uint8_t)i} });
          nbd.get_arrow_map()[x][(uint8_t)i] = 0;
        }
        for (int i = y - 1; i >= 0; --i) {
          if (bd.is_obstacle(x, i)) break;
          nbd.get_arrow_map()[x][i] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{x,(uint8_t)i} });
          nbd.get_arrow_map()[x][i] = 0;
        }
        for (int i = x + 1, j = y + 1; i < 8 && j < 8; ++i, ++j) {
          if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{(uint8_t)i,(uint8_t)j} });
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
        }
        for (int i = x - 1, j = y + 1; i >= 0 && j < 8; --i, ++j) {
          if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{(uint8_t)i,(uint8_t)j} });
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
        }
        for (int i = x + 1, j = y - 1; i < 8 && j >= 0; ++i, --j) {
          if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{(uint8_t)i,(uint8_t)j} });
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
        }
        for (int i = x - 1, j = y - 1; i >= 0 && j >= 0; --i, --j) {
          if (bd.is_obstacle(i, j)) break;
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
          res.push_back({ nbd, mv.from,mv.to,board::piece{(uint8_t)i,(uint8_t)j} });
          nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
        }
      }
      return res;
    }
    std::vector<move_action> _possible_moves(board::piece amazon) {
      std::vector<move_action> res;
      uint8_t x = amazon.x();
      uint8_t y = amazon.y();
      board nbd = bd;
      nbd.get_arrow_map()[x][y] = 0;
      for (int i = x + 1; i < 8; ++i) {
        if (bd.is_obstacle((uint8_t)i, y)) break;
        nbd.get_arrow_map()[(uint8_t)i][y] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{(uint8_t)i,y},board::piece{} });
        nbd.get_arrow_map()[(uint8_t)i][y] = 0;
      }
      for (int i = x - 1; i >= 0; --i) {
        if (bd.is_obstacle((uint8_t)i, y)) break;
        nbd.get_arrow_map()[(uint8_t)i][y] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{(uint8_t)i,y},board::piece{} });
        nbd.get_arrow_map()[(uint8_t)i][y] = 0;
      }
      for (int i = y + 1; i < 8; ++i) {
        if (bd.is_obstacle(x, (uint8_t)i)) break;
        nbd.get_arrow_map()[x][(uint8_t)i] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{x,(uint8_t)i},board::piece{} });
        nbd.get_arrow_map()[x][(uint8_t)i] = 0;
      }
      for (int i = y - 1; i >= 0; --i) {
        if (bd.is_obstacle(x, (uint8_t)i)) break;
        nbd.get_arrow_map()[x][(uint8_t)i] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{x,(uint8_t)i},board::piece{} });
        nbd.get_arrow_map()[x][(uint8_t)i] = 0;
      }
      for (int i = x + 1, j = y + 1; i < 8 && j < 8; ++i, ++j) {
        if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{(uint8_t)i,(uint8_t)j},board::piece{} });
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
      }
      for (int i = x - 1, j = y + 1; i >= 0 && j < 8; --i, ++j) {
        if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{(uint8_t)i,(uint8_t)j},board::piece{} });
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
      }
      for (int i = x + 1, j = y - 1; i < 8 && j >= 0; ++i, --j) {
        if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{(uint8_t)i,(uint8_t)j},board::piece{} });
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
      }
      for (int i = x - 1, j = y - 1; i >= 0 && j >= 0; --i, --j) {
        if (bd.is_obstacle((uint8_t)i, (uint8_t)j)) break;
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 1;
        res.push_back({ nbd,board::piece{x,y},board::piece{(uint8_t)i,(uint8_t)j},board::piece{} });
        nbd.get_arrow_map()[(uint8_t)i][(uint8_t)j] = 0;
      }
      return _possible_arrow_placements(res);
    }
  };
}
# 5 "main.cpp" 2
# 1 "interactor.hpp" 1
# 9 "interactor.hpp"
namespace yrq {
  class interactor {
  public:
    interactor() {};
    ~interactor() {};
    void generate_output(move_action mv) {
      std::cout << (int)mv.from.x() << ' ' << (int)mv.from.y() << ' ' << (int)mv.to.x() << ' ' << (int)mv.to.y() << ' ' << (int)mv.arrow.x() << ' ' << (int)mv.arrow.y() << std::endl;
    }
    std::vector<move_action> parse_input() {
      std::vector<move_action> mvs;
      int turn, from_x, from_y, to_x, to_y, arrow_x, arrow_y;
      bool i_am_black = false;
      board bd;
      bd.get_arrow_map() = bitmap(0x2400810000810024);
      std::cin >> turn;
      for (int i = 0; i < turn; ++i) {
        std::cin >> from_x >> from_y >> to_x >> to_y >> arrow_x >> arrow_y;
        if (from_x == -1);
        else {
          auto mv = make_move(bd, board::piece(from_x, from_y), board::piece(to_x, to_y), board::piece(arrow_x, arrow_y));
          mvs.push_back(mv);
        }
        if (i < turn - 1) {
          std::cin >> from_x >> from_y >> to_x >> to_y >> arrow_x >> arrow_y;
          if (from_x >= 0) {
            auto mv = make_move(bd, board::piece(from_x, from_y), board::piece(to_x, to_y), board::piece(arrow_x, arrow_y));
            mvs.push_back(mv);
          }
        }
      }
      return mvs;
    }
  };
}
# 6 "main.cpp" 2




using namespace yrq;
using namespace std;
int main(int argc, char** argv) {
  interactor ita;
  greedy_searcher scher;
  auto r = ita.parse_input();
  auto mv = scher.search_and_select(r);
  ita.generate_output(mv);
  /*

  if (argc < 2) cerr << "please select a dataset directory" << endl;

  filesystem::directory_iterator fs_dir_it(argv[1]);

  vector<filesystem::directory_entry> files;

  for (const auto& file : fs_dir_it)

    files.push_back(file);

  sort(files.begin(), files.end());

  for (const auto& file : files) {

    if (!file.is_regular_file()) continue;

    cout << "+--------------------------------------------------------------------------------------+" << std::endl;

    cout << "FILE " << file.path() << endl;

    emit_key_value("{", ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", true);

    emit_key_value("file", file.path(), true);

    string info;

    bit_matrix_file result;

    try {

      info = bit_matrix_filename_parse(std::filesystem::path(file.path()).stem().string());

      cout << info;

      result = bit_matrix_file_load(file.path().string());

    }

    catch (const std::exception&) {

      return 1;

    }

    board bd(result.bd);

    cout << "BOARD" << std::endl;

    bd.get_queen_map().output();

    evaluator ev(bd);

    for (int j = 0; j < 4; ++j)

      ev.players[0][j] = board::piece(result.xy[j][0], result.xy[j][1]);

    for (int j = 0; j < 4; ++j)

      ev.players[1][j] = board::piece(result.xy[j + 4][0], result.xy[j + 4][1]);

    ev.evaluate();

    emit_key_value("}", ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", true);

    cout << "+--------------------------------------------------------------------------------------+" << std::endl;

  }

  */
# 54 "main.cpp"
}
