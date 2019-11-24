#ifndef SEARCHER_H
#define SEARCHER_H
#include "board.hpp"
#include "evaluator.hpp"
#include "utility.hpp"
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

#endif