#include <chrono>
#include "board.hpp"
#include <Windows.h>

using namespace yrq;
int main() {
  bitmap queen(0xFF88888888888888);
  bitmap arrow;
  uint64_t no_use = 0;
  board bd(queen, arrow), bdcp = bd;
  std::cout << "ORIGINAL BIT_TABLE:" << std::endl;
  queen.output();
  std::cout << "ALL POSSIBLE MOVES AT (" << 3 << "," << 3 << "):" << std::endl;
  bd.accessible(3, 3).output();

  LARGE_INTEGER t1, t2, freq;
  QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&t1);
  for (int j = 0; j < 1e7; ++j) {
    int tmp = j & 0x3F;
    bdcp.get_queen_map().set(tmp, false);
    bdcp.get_queen_map().set(tmp % 8, true);
    no_use += bdcp.accessible_raw(tmp);
  }
  QueryPerformanceCounter(&t2);
  std::cout << (double)(t2.QuadPart - t1.QuadPart) / freq.QuadPart << "s" << std::endl << no_use;
}