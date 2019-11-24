#ifndef INTERACTOR_H
#define INTERACTOR_H
#include "board.hpp"
#include "evaluator.hpp"
#include "utility.hpp"
#include "bitmap.hpp"
#include "searcher.hpp"

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
#endif