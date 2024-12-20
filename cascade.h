#pragma once

#include <climits>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>
#include <random>
#include <set>
#include <memory>
#include <stdio.h>
#include <atomic>
#include <cstdint>

#include "xxhash.h"
#include "bloom_filter.hpp"

/*Below Code Taken (with modification) from https://github.com/FastFilter/fastfilter_cpp*/
#include "random.h"
#include "timing.h"
#include "ribbon_impl.h"


#include "hashutil.h"

#include "sequence.h"

#include "primitives.h"


using namespace std;
using namespace hashing;
using namespace ribbon;

// #define CONTAIN_ATTRIBUTES  __attribute__ ((noinline))

uint64_t generate_random_number() {
    thread_local std::mt19937_64 generator(std::random_device{}());
    thread_local std::uniform_int_distribution<uint64_t> distribution;

    return distribution(generator);
}

struct XXH3Hash {
    static uint64_t hash(const unsigned char* data, std::size_t length, uint64_t seed) {
        // XXH3_64bits_withSeed returns a 64-bit hash value
        auto h = XXH3_64bits_withSeed(data, length, seed);
        h ^= generate_random_number();
        return h;
    }
};

struct RibbonTsHomog {
  static constexpr bool kIsFilter = true;
  static constexpr bool kHomogeneous = true;
  static constexpr bool kFirstCoeffAlwaysOne = true;
  static constexpr bool kUseSmash = false;
  using CoeffRow = uint64_t;
  using Hash = uint64_t;
  using Key = std::vector<uint8_t>;
  using Seed = XXH64_hash_t;
  using Index = size_t;
  using ResultRow = uint32_t;
  static constexpr bool kAllowZeroStarts = false;

  static Hash HashFn(const Key& input, Seed raw_seed) {
    auto h = XXH3_64bits_withSeed(input.data(), input.size(), raw_seed);
    h ^= generate_random_number();
    return h;
  }
};


class HomogRibbonFilter {
  using TS = RibbonTsHomog;
  using Key = std::vector<uint8_t>;
  IMPORT_RIBBON_IMPL_TYPES(TS);

  uint32_t kNumColumns;
  size_t num_slots;
  size_t bytes;
  unique_ptr<char[]> ptr;
  InterleavedSoln soln;
  Hasher hasher;
public:
  HomogRibbonFilter(size_t add_count, uint32_t k)
      : kNumColumns(k),
        num_slots(InterleavedSoln::RoundUpNumSlots((size_t)(((4.0 + k * 0.25) / (8.0 * sizeof(uint64_t)) + 1) * add_count))),
        bytes(static_cast<size_t>((num_slots * k + 7) / 8)),
        ptr(new char[bytes]),
        soln(ptr.get(), bytes, k) {}

  void AddAll(const parlay::sequence<Key> keys, const size_t start, const size_t end) {
    Banding b(num_slots);
    (void)b.AddRange(keys.begin() + start, keys.begin() + end);
    soln.BackSubstFrom(b);
  }
  bool Contain(Key key) const {
    return soln.FilterQuery(key, hasher);
  }
  size_t SizeInBytes() const {
    return bytes;
  }
};

template<typename Table>
struct Table_API {};

template<>
struct Table_API<HomogRibbonFilter> {
  using Key = std::vector<uint8_t>;
  using Table = HomogRibbonFilter;
  static Table ConstructFromAddCount(size_t add_count, double k) { return Table(add_count, std::ceil(k)); }
  static void Add(Key key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const parlay::sequence<Key> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(Key key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  static bool Contain(Key key, const Table * table) __attribute__ ((noinline)) {
    return table->Contain(key);
  }

  static size_t SizeInBytes(Table* table) {
    return table->SizeInBytes();
  }
};

/*End of copied code*/

class BloomFilter {
    using Key = std::vector<uint8_t>;

public:
    bloom_parameters parameters;
    bloom_filter<XXH3Hash> filter;

    BloomFilter(bloom_parameters p) : parameters(p), filter(bloom_filter<XXH3Hash>(p)) {}

    void AddAll(const parlay::sequence<Key>& keys) {
        for (const auto& key : keys) {
            if (!key.empty()) {
                filter.insert(reinterpret_cast<const unsigned char*>(key.data()), key.size());
            }
        }
    }

    bool Contain(const Key& key) const {
        if (key.empty()) {
            return false;
        }
        return filter.contains(reinterpret_cast<const unsigned char*>(key.data()), key.size());
    }

    size_t SizeInBytes() const {
        return (filter.size() + 7) / 8;
    }
};


template<>
struct Table_API<BloomFilter> {
  using Key = std::vector<uint8_t>;
  using Table = BloomFilter;

  static Table ConstructFromAddCount(size_t add_count, double k) {
    bloom_parameters p;
    p.projected_element_count = add_count;
    p.false_positive_probability = std::pow(2.0, -k);
    p.random_seed = generate_random_number();
    p.compute_optimal_parameters();
    return BloomFilter(p);
  }

  static void Add(Key key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const parlay::sequence<Key> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys);
  }

  static void Remove(Key key, Table * table) {
    throw std::runtime_error("Unsupported");
  }

  static bool Contain(Key key, const Table * table) __attribute__ ((noinline)) {
    return table->Contain(key);
  }

  static size_t SizeInBytes(Table* table) {
    return table->SizeInBytes();
  }
};


template<typename Table>
void build_cascade(int level, std::vector<Table>* cascade, const parlay::sequence<std::vector<uint8_t>> R, const parlay::sequence<std::vector<uint8_t>> S) {
    if (R.empty()) return;
    double k = (level == 0 && R.size() < S.size()) ? (-std::log2(double(R.size())/(std::sqrt(2)*double(S.size())))) : 1; 
    Table filter = Table_API<Table>::ConstructFromAddCount(R.size(),k);
    Table_API<Table>::AddAll(R, 0, R.size(), &filter);
    auto R_next = parlay::filter(S, [&] (auto&& s) {
      return Table_API<Table>::Contain(s, &filter);
    });
    cascade->push_back(std::move(filter));
    build_cascade<Table>(level + 1, cascade, R_next, R);
}

template <typename Table>
size_t cascade_size_in_bytes(std::vector<Table>* cascade) {
    size_t size = 0;
    for (int i = 0; i < cascade->size(); i++) {
        size += cascade->at(i).SizeInBytes();
    }

    return size;
}
