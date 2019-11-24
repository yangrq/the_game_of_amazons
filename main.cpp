#include "board.hpp"
#include "evaluator.hpp"
#include "utility.hpp"
#include "searcher.hpp"
#include "interactor.hpp"
#ifndef ALL_IN_ONE
#include <chrono>
#endif

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
}