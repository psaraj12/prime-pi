// ============================================================================
// pi_query_wasm.cpp
// WASM module for π(x) queries against a precomputed prefix-sum table.
//
// Compile with Emscripten:
//   emcc -O3 -s WASM=1 -s EXPORTED_FUNCTIONS='["_init_table","_query_pi","_get_table_end","_get_checkpoint_interval","_get_num_checkpoints","_clear_query_cache","_get_cache_hits","_get_cache_misses","_malloc","_free"]' \
//        -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap","getValue"]' \
//        -s ALLOW_MEMORY_GROWTH=1 -s INITIAL_MEMORY=16777216 \
//        -o pi_query.js pi_query_wasm.cpp
// ============================================================================

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

using u64 = uint64_t;
using u32 = uint32_t;

// ============================================================================
// Wheel-30 tables
// ============================================================================

static constexpr int W30_RES[8] = {1, 7, 11, 13, 17, 19, 23, 29};
static constexpr int W30_IDX[30] = {
    -1, 0,-1,-1,-1,-1,-1, 1,-1,-1,
    -1, 2,-1, 3,-1,-1,-1, 4,-1, 5,
    -1,-1,-1, 6,-1,-1,-1,-1,-1, 7
};

static int g_k_target[30][8];

static void build_mod_tables() {
    memset(g_k_target, -1, sizeof(g_k_target));
    for (int p_mod = 0; p_mod < 30; p_mod++) {
        if (W30_IDX[p_mod] < 0 && p_mod != 1) continue;
        for (int ri = 0; ri < 8; ri++) {
            int res = W30_RES[ri];
            for (int m = 0; m < 30; m++) {
                if ((p_mod * m) % 30 == res) {
                    g_k_target[p_mod][ri] = m;
                    break;
                }
            }
        }
    }
}

// ============================================================================
// Base primes (small sieve up to sqrt of max query)
// ============================================================================

static std::vector<u32> g_base_primes;
static u64 g_base_primes_limit = 0;

static void build_base_primes(u64 limit) {
    if (limit <= g_base_primes_limit) return;  // already built
    u64 sieve_size = limit + 1;
    std::vector<bool> is_composite(sieve_size, false);
    for (u64 i = 2; i * i <= limit; i++)
        if (!is_composite[i])
            for (u64 j = i * i; j <= limit; j += i)
                is_composite[j] = true;
    g_base_primes.clear();
    for (u64 i = 2; i <= limit; i++)
        if (!is_composite[i])
            g_base_primes.push_back((u32)i);
    g_base_primes_limit = limit;
}

// ============================================================================
// Wheel-30 sieve for partial block
// ============================================================================

struct W30Sieve {
    u64* bits;
    u64 base, lo, hi, num_groups, num_words;
    bool owns_memory;

    W30Sieve() : bits(nullptr), owns_memory(false) {}
    ~W30Sieve() { if (owns_memory && bits) free(bits); }

    inline void mark_composite(u64 n) {
        u64 offset = n - base;
        u64 group = offset / 30;
        int bit = W30_IDX[offset % 30];
        if (bit < 0) return;
        u64 flat_bit = group * 8 + bit;
        bits[flat_bit / 64] &= ~(1ULL << (flat_bit % 64));
    }

    void init(u64 seg_lo, u64 seg_hi_inclusive) {
        lo = seg_lo;
        hi = seg_hi_inclusive + 1;
        base = (lo / 30) * 30;
        num_groups = (hi - base + 29) / 30;
        u64 total_bits = num_groups * 8;
        num_words = (total_bits + 63) / 64;

        if (owns_memory && bits) free(bits);
        bits = (u64*)calloc(num_words, sizeof(u64));
        owns_memory = true;

        // Set all bits to 1 (all prime)
        memset(bits, 0xFF, num_words * sizeof(u64));

        // Mask trailing bits
        u64 tail = total_bits % 64;
        if (tail != 0)
            bits[num_words - 1] &= (1ULL << tail) - 1;

        // Clear bits outside [lo, hi)
        u64 first_group_start = base;
        if (first_group_start < lo) {
            for (int b = 0; b < 8; b++) {
                u64 n = first_group_start + W30_RES[b];
                if (n < lo) {
                    u64 flat_bit = (u64)b;
                    bits[flat_bit / 64] &= ~(1ULL << (flat_bit % 64));
                }
            }
        }
        u64 last_group_start = base + (num_groups - 1) * 30;
        for (int b = 0; b < 8; b++) {
            u64 n = last_group_start + W30_RES[b];
            if (n > seg_hi_inclusive) {
                u64 flat_bit = (num_groups - 1) * 8 + b;
                bits[flat_bit / 64] &= ~(1ULL << (flat_bit % 64));
            }
        }
        if (base == 0) bits[0] &= ~1ULL;  // 1 is not prime
    }

    u64 count_primes() const {
        u64 count = 0;
        for (u64 w = 0; w < num_words; w++)
            count += __builtin_popcountll(bits[w]);
        return count;
    }

    u64 count_primes_to(u64 limit_inclusive) const {
        if (!bits || limit_inclusive < lo) return 0;
        if (limit_inclusive >= hi - 1) return count_primes();

        u64 rel = limit_inclusive - base;
        u64 groups_before = rel / 30;
        int rem = (int)(rel % 30);

        u64 bit_count = groups_before * 8;
        for (int i = 0; i < 8; i++) {
            if (W30_RES[i] <= rem) bit_count++;
            else break;
        }

        u64 full_words = bit_count / 64;
        int tail_bits = (int)(bit_count % 64);
        u64 count = 0;

        for (u64 w = 0; w < full_words; w++)
            count += __builtin_popcountll(bits[w]);

        if (tail_bits != 0) {
            u64 mask = (1ULL << tail_bits) - 1;
            count += __builtin_popcountll(bits[full_words] & mask);
        }
        return count;
    }
};

static void sieve_segment_w30(W30Sieve& seg) {
    u64 base = seg.base;
    u64 seg_end = base + seg.num_groups * 30;

    for (size_t i = 3; i < g_base_primes.size(); i++) {
        u64 p = g_base_primes[i];
        u64 min_start = (base > p * p) ? base : p * p;
        if (min_start >= seg_end) continue;
        int p_mod = (int)(p % 30);
        u64 step = 30 * p;
        u64 k_min = (min_start + p - 1) / p;
        u64 k_mod_base = k_min % 30;

        for (int r = 0; r < 8; r++) {
            int target_k_mod = g_k_target[p_mod][r];
            u64 k = k_min + ((u64)(target_k_mod - (int)k_mod_base + 30) % 30);
            for (u64 n = k * p; n < seg_end; n += step)
                seg.mark_composite(n);
        }
    }
}

// ============================================================================
// Table data
// ============================================================================

static u64 g_table_end = 0;
static u64 g_checkpoint_interval = 0;
static u64 g_num_checkpoints = 0;
static u64* g_checkpoints = nullptr;

// Last-query cache: reuses the previous sieved partial checkpoint interval.
// This helps mobile when users issue nearby/repeated queries while typing.
static W30Sieve g_cached_seg;
static bool g_cache_valid = false;
static u64 g_cache_ci = 0;
static u64 g_cache_hi = 0;
static u64 g_cache_hits = 0;
static u64 g_cache_misses = 0;

static void clear_cache_internal() {
    g_cache_valid = false;
    g_cache_ci = 0;
    g_cache_hi = 0;
}

// ============================================================================
// Exported functions (called from JavaScript)
// ============================================================================

extern "C" {

// Initialize the table from raw binary data
// data points to the full pi.dat file content (header + checkpoints)
// Returns: 1 on success, 0 on failure
int init_table(const uint8_t* data, int data_len) {
    if (data_len < 64) return 0;

    // Verify magic
    if (memcmp(data, "PIFIN01", 7) != 0) return 0;

    // Read header
    memcpy(&g_table_end, data + 8, 8);
    memcpy(&g_checkpoint_interval, data + 16, 8);
    memcpy(&g_num_checkpoints, data + 24, 8);

    // Verify data size
    u64 expected_size = 64 + g_num_checkpoints * 8;
    if ((u64)data_len < expected_size) return 0;

    // Copy checkpoint data
    if (g_checkpoints) free(g_checkpoints);
    g_checkpoints = (u64*)malloc(g_num_checkpoints * sizeof(u64));
    if (!g_checkpoints) return 0;
    memcpy(g_checkpoints, data + 64, g_num_checkpoints * 8);

    clear_cache_internal();
    g_cache_hits = 0;
    g_cache_misses = 0;

    // Build base primes for sieving partial blocks
    u64 sqrt_end = (u64)sqrt((double)g_table_end) + 1;
    build_base_primes(sqrt_end);
    build_mod_tables();

    return 1;
}

// Query π(x)
// Returns the count of primes <= x
// x_lo and x_hi are the low and high 32 bits of x (since WASM is 32-bit)
double query_pi(double x_dbl) {
    u64 x = (u64)x_dbl;

    if (x < 2) return 0;
    if (x > g_table_end || !g_checkpoints) return -1;

    // Which checkpoint interval?
    u64 ci = x / g_checkpoint_interval;
    u64 cp_end = (ci + 1) * g_checkpoint_interval - 1;
    if (cp_end > g_table_end) cp_end = g_table_end;

    // Exact checkpoint hit
    if (x == cp_end) return (double)g_checkpoints[ci];

    // Base count from previous checkpoint
    u64 base_count = (ci > 0) ? g_checkpoints[ci - 1] : 0;

    // Small primes for interval 0
    u64 small_prime_count = 0;
    if (ci == 0) {
        if (x >= 2) small_prime_count++;
        if (x >= 3) small_prime_count++;
        if (x >= 5) small_prime_count++;
    }

    // Sieve partial block, with conservative last-query cache.
    u64 interval_start = ci * g_checkpoint_interval;
    u64 seg_count = 0;

    if (g_cache_valid && g_cache_ci == ci && x <= g_cache_hi) {
        ++g_cache_hits;
        seg_count = g_cached_seg.count_primes_to(x);
    } else {
        ++g_cache_misses;
        g_cached_seg.init(interval_start, x);
        sieve_segment_w30(g_cached_seg);
        g_cache_valid = true;
        g_cache_ci = ci;
        g_cache_hi = x;
        seg_count = g_cached_seg.count_primes();
    }

    return (double)(base_count + small_prime_count + seg_count);
}

// Get the table's maximum value
double get_table_end() {
    return (double)g_table_end;
}

// Get checkpoint interval
double get_checkpoint_interval() {
    return (double)g_checkpoint_interval;
}

// Get number of checkpoints
double get_num_checkpoints() {
    return (double)g_num_checkpoints;
}

void clear_query_cache() {
    clear_cache_internal();
}

double get_cache_hits() {
    return (double)g_cache_hits;
}

double get_cache_misses() {
    return (double)g_cache_misses;
}

}  // extern "C"
