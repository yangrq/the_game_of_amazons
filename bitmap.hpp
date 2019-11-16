#ifndef BITMAP_H
#define BITMAP_H
#include <cstdint>
#include <bitset>
#include <iostream>

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
#ifdef _DEBUG
      if (x < lower_bound || x > upper_bound) throw std::exception("index overflow");
#endif
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
    void set(size_t _pos, bool _val) noexcept {
      is_valid(_pos, 0, 63);
      return _set(_pos, _val);
    }
    bool get(size_t _pos) noexcept {
      is_valid(_pos, 0, 63);
      return get(_pos);
    }
  };
}
#endif