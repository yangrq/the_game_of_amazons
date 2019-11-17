#include "board.hpp"
#include "evaluator.hpp"
#include <chrono>

using namespace yrq;
int main() {
  /*
  bitmap amazon(0xFF88888888888888);
  bitmap arrow;
  uint64_t no_use = 0;
  board bd(amazon, arrow), bdcp = bd;
  std::cout << "ORIGINAL BIT_TABLE:" << std::endl;
  amazon.output();
  std::cout << "ALL POSSIBLE MOVES AT (" << 3 << "," << 3 << "):" << std::endl;
  bd.accessible(3, 3).output();
  std::chrono::high_resolution_clock clk;

  auto t1 = clk.now();
  for (int j = 0; j < 1e7; ++j) {
    int tmp = j & 0x3F;
    bdcp.get_queen_map().set(tmp, false);
    bdcp.get_queen_map().set(tmp % 8, true);
    no_use += bdcp.accessible_raw(tmp);
  }
  auto t2 = clk.now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms" << std::endl << no_use << std::endl;
  */
  bitmap dis(0x201220223908557f);
  std::cout << "BOARD OUTPUT:" << std::endl;
  dis.output();
}