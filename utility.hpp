#ifndef UTILITY_H
#define UTILITY_H
#ifndef ALL_IN_ONE
#include <cstdint>
#include <bitset>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <regex>
#include <filesystem>
#endif

#if defined(UTILITY) || defined(_DEBUG)

std::ofstream ofs("records.log", std::ios::trunc);

//反转64bit整数
inline uint64_t bit_reverse(const uint64_t raw) {
  std::bitset<64> v = raw;
  for (uint64_t i = 0; i < 32; ++i)
    v[i] = v[63ull - i];
  return v.to_ullong();
}

struct bit_matrix_file {
  std::string id;
  int xy[8][2] = { 0 };
  uint64_t bd = 0;
};
template <typename Ty>
void emit_key_value(std::string key, Ty value, bool newline = false) {
  ofs.flush();
  ofs << key << ":" << value << "  ";
  ofs.flush();
  if (newline) ofs << std::endl;
  ofs.flush();
}

//读取bit_matrix(*.bm)文件名信息
inline std::string bit_matrix_filename_parse(std::string filename) {
  std::regex pattern(
    R"(([a-zA-z_]\w*)\-(\d+)\-(\d+)\-(0|1)\-wins)",
    std::regex::ECMAScript
    | std::regex::icase
    | std::regex::optimize
  );
  std::smatch match;
  std::string info;
  if (std::regex_match(filename, match, pattern)) {
    info += "GAME " + std::string(match[1]) + "\n";
    info += "TURN " + std::string(match[2]) + "/" + std::string(match[3]) + "\n";
    info += "PLAYER" + std::string(match[4]) + " WINS\n";
  }
  else {
    std::cout << "file name error!" << std::endl;
    throw std::exception();
  }
  return info;
}

//读入bit_matrix(*.bm)文件
inline bit_matrix_file bit_matrix_file_load(std::string filename) {
  std::ifstream ifs(filename);
  std::string all;
  std::regex pattern(
    R"(([a-zA-z_]\w*)\{\s{0,2}([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})([01][01][01][01][01][01][01][01]\s{0,2})\s{0,2}\}\s{0,2}pieces\s{0,2}(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2})(\{\d+,\d+\}\s{0,2}))",
    std::regex::ECMAScript
    | std::regex::icase
    | std::regex::optimize
  );
  std::smatch match;
  bit_matrix_file bm;
  if (!ifs.is_open())
    std::cout << "failed to open " << filename << std::endl;
  else
    for (;;) {
      std::string buf;
      if (!std::getline(ifs, buf)) break;
      all += buf;
    }
  if (std::regex_match(all, match, pattern)) {
    bm.id = match[1];
    std::bitset<64> bs;
    for (uint64_t i = 0; i < 8; ++i) {
      std::string s = match[i + 2];
      for (int j = 0; j < 8; ++j)
        bs[i * 8 + j] = s[j] - '0';
    }
    bm.bd = bs.to_ullong();
    for (uint64_t i = 0; i < 8; ++i) {
      std::string s = match[i + 10];
      int x, y;
      sscanf_s(s.c_str(), "{%d,%d}", &x, &y);
      bm.xy[i][0] = x;
      bm.xy[i][1] = y;
    }
  }
  else {
    std::cout << "file format error!" << std::endl;
    throw std::exception();
  }
  return bm;
}
#endif
#endif