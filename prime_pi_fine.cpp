// ============================================================================
// prime_pi_fine.cpp
// π(x) calculator with fine-grained prefix sum table
//
// Stores cumulative π at every 30,000-number interval (1000 groups of 30).
// Query: look up nearest checkpoint, sieve at most 30K numbers. Microseconds.
//
// Build uses large sieve blocks (configurable) for cache efficiency, but
// samples the count at fine intervals within each block.
//
// Modes:
//   build  <end> [threads] [table_file]    — sieve and build table
//   extend <new_end> [threads] [table_file] — extend existing table
//   query  <x> [table_file]                — instant π(x) lookup
//   info   [table_file]                    — show table metadata
//
// Build:
//   g++ -O3 -march=native -fopenmp -std=c++17 -o prime_pi_fine prime_pi_fine.cpp
//
// Examples:
//   ./prime_pi_fine build  1000000000000 12 pi.dat
//   ./prime_pi_fine extend 2000000000000 12 pi.dat
//   ./prime_pi_fine query  123456789012 pi.dat
//   ./prime_pi_fine info   pi.dat
// ============================================================================

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <omp.h>

using u64 = uint64_t;
using u32 = uint32_t;
using i64 = int64_t;

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
// Base primes
// ============================================================================

static std::vector<u32> g_base_primes;

static void build_base_primes(u64 limit) {
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
}

// ============================================================================
// Wheel-30 sieve (with fine-grained count sampling)
// ============================================================================

struct W30Sieve {
    std::vector<u64> bits;
    u64 base, lo, hi, num_groups, num_words;

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
        bits.assign(num_words, ~0ULL);
        u64 tail = total_bits % 64;
        if (tail != 0)
            bits[num_words - 1] &= (1ULL << tail) - 1;

        u64 first_group_start = base;
        if (first_group_start < lo) {
            for (int b = 0; b < 8; b++) {
                u64 n = first_group_start + W30_RES[b];
                if (n < lo) {
                    u64 flat_bit = b;
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
        if (base == 0) bits[0] &= ~1ULL;
    }

    // Count primes in entire sieve
    u64 count_primes() const {
        u64 count = 0;
        for (u64 w = 0; w < num_words; w++)
            count += __builtin_popcountll(bits[w]);
        return count;
    }

    // Count primes in groups [0, end_group) — for partial counting
    u64 count_primes_groups(u64 end_group) const {
        if (end_group >= num_groups) return count_primes();
        u64 end_bit = end_group * 8;
        u64 full_words = end_bit / 64;
        u64 count = 0;
        for (u64 w = 0; w < full_words; w++)
            count += __builtin_popcountll(bits[w]);
        // Partial last word
        u64 remaining_bits = end_bit % 64;
        if (remaining_bits > 0)
            count += __builtin_popcountll(bits[full_words] & ((1ULL << remaining_bits) - 1));
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
// Fine-grained prefix table
//
// Checkpoint every CHECKPOINT_GROUPS groups of 30 = every CHECKPOINT_INTERVAL numbers.
// checkpoint[i] = π(i * CHECKPOINT_INTERVAL - 1) = cumulative primes up to that point.
//
// For block 0, this includes primes 2, 3, 5.
// ============================================================================

static constexpr u64 CHECKPOINT_GROUPS = 100000;     // 1000000 groups of 30
static constexpr u64 CHECKPOINT_INTERVAL = CHECKPOINT_GROUPS * 30;  // = 300000 numbers
static constexpr u64 BLOCK_MULTIPLIER = 1;      // tune this
// Sieve block size for cache efficiency (separate from checkpoint interval)
// Keep this an exact multiple of CHECKPOINT_GROUPS so checkpoint intervals
// never cross sieve-block boundaries during build/extend.
static constexpr u64 SIEVE_BLOCK_GROUPS = CHECKPOINT_GROUPS * BLOCK_MULTIPLIER;  // 256 checkpoints/block = 7.68M numbers
static constexpr u64 SIEVE_BLOCK_SIZE = SIEVE_BLOCK_GROUPS * 30;

static_assert(SIEVE_BLOCK_GROUPS % CHECKPOINT_GROUPS == 0,
              "SIEVE_BLOCK_GROUPS must be an exact multiple of CHECKPOINT_GROUPS");

static inline bool is_checkpoint_aligned_end(u64 x) {
    return (x % CHECKPOINT_INTERVAL) == (CHECKPOINT_INTERVAL - 1);
}

static inline u64 next_checkpoint_aligned_end(u64 x) {
    u64 rem = x % CHECKPOINT_INTERVAL;
    if (rem == CHECKPOINT_INTERVAL - 1) return x;
    return x + (CHECKPOINT_INTERVAL - 1 - rem);
}

struct PiTableFine {
    u64 end;                       // table covers [0, end]
    u64 num_checkpoints;           // number of checkpoint entries
    std::vector<u64> checkpoints;  // checkpoints[i] = π(checkpoint_end(i))

    // The number that checkpoint i covers up to (inclusive)
    u64 checkpoint_end(u64 i) const {
        return std::min((i + 1) * CHECKPOINT_INTERVAL - 1, end);
    }

    // Save to binary file
    bool save(const char* filename) const {
        FILE* f = fopen(filename, "wb");
        if (!f) return false;

        const char magic[] = "PIFIN01";
        fwrite(magic, 1, 8, f);
        fwrite(&end, sizeof(u64), 1, f);
        u64 ci = CHECKPOINT_INTERVAL;
        fwrite(&ci, sizeof(u64), 1, f);
        fwrite(&num_checkpoints, sizeof(u64), 1, f);

        // Pad header to 64 bytes
        char pad[64 - 8 - 3 * 8] = {0};
        fwrite(pad, 1, sizeof(pad), f);

        // Write checkpoint array
        fwrite(checkpoints.data(), sizeof(u64), num_checkpoints, f);
        fclose(f);
        return true;
    }

    bool load(const char* filename) {
        FILE* f = fopen(filename, "rb");
        if (!f) return false;

        char magic[8];
        if (fread(magic, 1, 8, f) != 8 || memcmp(magic, "PIFIN01", 8) != 0) {
            fclose(f);
            return false;
        }

        u64 ci;
        if (fread(&end, sizeof(u64), 1, f) != 1) { fclose(f); return false; }
        if (fread(&ci, sizeof(u64), 1, f) != 1) { fclose(f); return false; }
        if (fread(&num_checkpoints, sizeof(u64), 1, f) != 1) { fclose(f); return false; }

        if (ci != CHECKPOINT_INTERVAL) {
            fprintf(stderr, "Error: table has checkpoint interval %llu, expected %llu\n",
                    (unsigned long long)ci, (unsigned long long)CHECKPOINT_INTERVAL);
            fclose(f);
            return false;
        }

        char pad[64 - 8 - 3 * 8];
        if (fread(pad, 1, sizeof(pad), f) != sizeof(pad)) { fclose(f); return false; }

        checkpoints.resize(num_checkpoints);
        if (fread(checkpoints.data(), sizeof(u64), num_checkpoints, f) != num_checkpoints) {
            fclose(f);
            return false;
        }

        fclose(f);
        return true;
    }

    // Query π(x) — O(1) lookup + tiny sieve
    u64 query(u64 x) const {
        if (x < 2) return 0;
        if (x > end) {
            fprintf(stderr, "Error: x=%llu exceeds table range %llu\n",
                    (unsigned long long)x, (unsigned long long)end);
            return 0;
        }

        // Which checkpoint interval does x fall in?
        u64 ci = x / CHECKPOINT_INTERVAL;

        // If x is exactly at a checkpoint boundary
        if (x == checkpoint_end(ci)) return checkpoints[ci];

        // Base count from previous checkpoint
        u64 base_count = (ci > 0) ? checkpoints[ci - 1] : 0;

        // Start of this checkpoint interval
        u64 interval_start = ci * CHECKPOINT_INTERVAL;

        // Need to count primes in [interval_start, x]
        // For interval 0, also count 2, 3, 5
        u64 small_prime_count = 0;
        if (ci == 0) {
            if (x >= 2) small_prime_count++;
            if (x >= 3) small_prime_count++;
            if (x >= 5) small_prime_count++;
        }

        // Tiny sieve — at most 30K numbers
        W30Sieve seg;
        seg.init(interval_start, x);
        sieve_segment_w30(seg);
        u64 seg_count = seg.count_primes();

        return base_count + small_prime_count + seg_count;
    }
};

// ============================================================================
// Build: sieve in large blocks, sample checkpoints within each block
// ============================================================================

static void do_build(u64 end, int threads, const char* table_file) {
    u64 requested_end = end;
    if (!is_checkpoint_aligned_end(end)) {
        end = next_checkpoint_aligned_end(end);
        printf("Requested end %llu is not checkpoint-aligned; using %llu instead.\n",
               (unsigned long long)requested_end,
               (unsigned long long)end);
    }

    omp_set_num_threads(threads);

    u64 num_checkpoints = end / CHECKPOINT_INTERVAL + 1;

    printf("BUILD MODE\n");
    printf("Target: [0, %llu]\n", (unsigned long long)end);
    printf("Checkpoint interval: %llu numbers\n", (unsigned long long)CHECKPOINT_INTERVAL);
    printf("Sieve block: %llu numbers (%.0f KB)\n",
           (unsigned long long)SIEVE_BLOCK_SIZE,
           (double)(SIEVE_BLOCK_GROUPS) / 1024.0);
    printf("Checkpoints: %llu\n", (unsigned long long)num_checkpoints);
    printf("Table size: %.1f MB\n", (double)(num_checkpoints * 8) / (1024.0 * 1024.0));
    printf("Threads: %d\n", threads);
    printf("Output: %s\n\n", table_file);

    u64 sqrt_end = (u64)sqrt((double)end) + 1;
    printf("Building base primes up to %llu... ", (unsigned long long)sqrt_end);
    fflush(stdout);
    build_base_primes(sqrt_end);
    printf("done. %zu base primes.\n\n", g_base_primes.size());
    build_mod_tables();

    // Each sieve block covers SIEVE_BLOCK_SIZE numbers.
    // Within each sieve block, we extract checkpoint counts at every
    // CHECKPOINT_INTERVAL boundary.
    u64 num_sieve_blocks = (end + SIEVE_BLOCK_SIZE) / SIEVE_BLOCK_SIZE;
    u64 checkpoints_per_sieve_block = SIEVE_BLOCK_SIZE / CHECKPOINT_INTERVAL;
    // = 7680000 / 30000 = 256 checkpoints per sieve block exactly

    // Allocate: each sieve block produces its own local counts, then we prefix-sum
    // Store per-block local checkpoint counts
    // local_counts[sieve_block][checkpoint_within_block] = count in that interval
    std::vector<std::vector<u64>> local_counts(num_sieve_blocks);

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_last_report = t_start;

    u64 batch_size = (u64)threads * 2;

    for (u64 batch_lo = 0; batch_lo < num_sieve_blocks; batch_lo += batch_size) {
        u64 batch_hi = std::min(batch_lo + batch_size, num_sieve_blocks);
        u64 actual = batch_hi - batch_lo;

        #pragma omp parallel for schedule(dynamic)
        for (u64 bi = 0; bi < actual; bi++) {
            u64 sb = batch_lo + bi;
            u64 sb_start = sb * SIEVE_BLOCK_SIZE;
            u64 sb_end = std::min(sb_start + SIEVE_BLOCK_SIZE - 1, end);

            // Sieve this block
            W30Sieve seg;
            seg.init(sb_start, sb_end);
            sieve_segment_w30(seg);

            // Extract checkpoint counts within this sieve block
            // How many checkpoint intervals overlap with this sieve block?
            u64 first_cp = sb_start / CHECKPOINT_INTERVAL;
            u64 last_cp = sb_end / CHECKPOINT_INTERVAL;
            u64 n_cp = last_cp - first_cp + 1;

            local_counts[sb].resize(n_cp);

            u64 prev_count = 0;
            for (u64 c = 0; c < n_cp; c++) {
                u64 cp_idx = first_cp + c;
                u64 cp_end_num = std::min((cp_idx + 1) * CHECKPOINT_INTERVAL - 1, end);

                // How many groups from sieve base to cp_end_num?
                u64 end_group;
                if (cp_end_num >= seg.base) {
                    end_group = (cp_end_num - seg.base) / 30 + 1;
                    if (end_group > seg.num_groups) end_group = seg.num_groups;
                } else {
                    end_group = 0;
                }

                u64 cum = seg.count_primes_groups(end_group);
                local_counts[sb][c] = cum - prev_count;
                prev_count = cum;
            }

            // Block 0: add primes 2, 3, 5 to first checkpoint
            if (sb == 0 && n_cp > 0) {
                u64 cp0_end = std::min(CHECKPOINT_INTERVAL - 1, end);
                u64 small = 0;
                if (cp0_end >= 2) small++;
                if (cp0_end >= 3) small++;
                if (cp0_end >= 5) small++;
                local_counts[0][0] += small;
            }
        }

        auto t_now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t_now - t_last_report).count();
        if (elapsed >= 5.0 || batch_hi == num_sieve_blocks) {
            double total_elapsed = std::chrono::duration<double>(t_now - t_start).count();
            double progress = (double)batch_hi / num_sieve_blocks * 100.0;
            u64 numbers_done = std::min(batch_hi * SIEVE_BLOCK_SIZE, end + 1);
            double rate = (double)numbers_done / total_elapsed / 1e6;
            double frac = (double)batch_hi / num_sieve_blocks;
            double eta = (frac > 0.001) ? total_elapsed * (1.0 - frac) / frac : 0;

            printf("\r  [%.1f%%] Sieve block %llu/%llu | %.1f M/s | elapsed %.0fs | ETA %.0fs   ",
                   progress,
                   (unsigned long long)batch_hi,
                   (unsigned long long)num_sieve_blocks,
                   rate, total_elapsed, eta);
            fflush(stdout);
            t_last_report = t_now;
        }
    }

    // Phase 2: flatten local counts and prefix-sum
    printf("\n  Building prefix sum (%llu entries)... ",
           (unsigned long long)num_checkpoints);
    fflush(stdout);

    PiTableFine table;
    table.end = end;
    table.num_checkpoints = num_checkpoints;
    table.checkpoints.resize(num_checkpoints);

    u64 running = 0;
    u64 cp_idx = 0;
    for (u64 sb = 0; sb < num_sieve_blocks && cp_idx < num_checkpoints; sb++) {
        for (u64 c = 0; c < local_counts[sb].size() && cp_idx < num_checkpoints; c++) {
            running += local_counts[sb][c];
            table.checkpoints[cp_idx] = running;
            cp_idx++;
        }
    }

    printf("done.\n");

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(t_end - t_start).count();

    printf("  Saving to %s... ", table_file);
    fflush(stdout);
    if (table.save(table_file))
        printf("done.\n");
    else
        printf("FAILED!\n");

    printf("\n");
    printf("════════════════════════════════════════\n");
    printf("  π(%llu) = %llu\n", (unsigned long long)end, (unsigned long long)running);
    printf("  Build time: %.3f s\n", total_time);
    printf("  Rate: %.1f M/s\n", (double)end / total_time / 1e6);
    printf("  Table: %llu checkpoints, %.1f MB\n",
           (unsigned long long)num_checkpoints,
           (double)(num_checkpoints * 8) / (1024.0 * 1024.0));
    printf("════════════════════════════════════════\n");
}

// ============================================================================
// Extend: load existing table, sieve only new range, append checkpoints
// ============================================================================

static void do_extend(u64 new_end, int threads, const char* table_file) {
    PiTableFine table;
    if (!table.load(table_file)) {
        fprintf(stderr, "Error: cannot load table from %s\n", table_file);
        return;
    }

    u64 old_end = table.end;
    if (!is_checkpoint_aligned_end(old_end)) {
        fprintf(stderr, "Error: existing table end %llu is not checkpoint-aligned; cannot safely extend it.\n",
                (unsigned long long)old_end);
        fprintf(stderr, "Please rebuild with this corrected version.\n");
        return;
    }

    u64 requested_new_end = new_end;
    if (!is_checkpoint_aligned_end(new_end)) {
        new_end = next_checkpoint_aligned_end(new_end);
        printf("Requested new_end %llu is not checkpoint-aligned; using %llu instead.\n",
               (unsigned long long)requested_new_end,
               (unsigned long long)new_end);
    }

    if (new_end <= old_end) {
        printf("Table already covers [0, %llu]. Nothing to extend.\n",
               (unsigned long long)old_end);
        return;
    }

    omp_set_num_threads(threads);

    printf("EXTEND MODE\n");
    printf("Existing: [0, %llu] — %llu checkpoints, π = %llu\n",
           (unsigned long long)old_end,
           (unsigned long long)table.num_checkpoints,
           (unsigned long long)table.checkpoints[table.num_checkpoints - 1]);
    printf("Extending to: %llu\n", (unsigned long long)new_end);
    printf("Threads: %d\n\n", threads);

    u64 sqrt_new = (u64)sqrt((double)new_end) + 1;
    printf("Building base primes up to %llu... ", (unsigned long long)sqrt_new);
    fflush(stdout);
    build_base_primes(sqrt_new);
    printf("done. %zu base primes.\n", g_base_primes.size());
    build_mod_tables();

    // New checkpoints needed
    u64 total_checkpoints = new_end / CHECKPOINT_INTERVAL + 1;
    u64 old_checkpoints = table.num_checkpoints;
    u64 new_checkpoints = total_checkpoints - old_checkpoints;

    // New sieve range starts right after last complete checkpoint
    u64 new_start = old_checkpoints * CHECKPOINT_INTERVAL;
    u64 new_range = new_end - new_start + 1;

    printf("New checkpoints: %llu\n", (unsigned long long)new_checkpoints);
    printf("New range: [%llu, %llu] (%llu numbers)\n",
           (unsigned long long)new_start,
           (unsigned long long)new_end,
           (unsigned long long)new_range);
    printf("\n");

    // Sieve new range in blocks, extract checkpoint counts
    u64 num_sieve_blocks = (new_range + SIEVE_BLOCK_SIZE - 1) / SIEVE_BLOCK_SIZE;
    std::vector<std::vector<u64>> local_counts(num_sieve_blocks);

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_last_report = t_start;

    u64 batch_size = (u64)threads * 2;

    for (u64 batch_lo = 0; batch_lo < num_sieve_blocks; batch_lo += batch_size) {
        u64 batch_hi = std::min(batch_lo + batch_size, num_sieve_blocks);
        u64 actual = batch_hi - batch_lo;

        #pragma omp parallel for schedule(dynamic)
        for (u64 bi = 0; bi < actual; bi++) {
            u64 sb = batch_lo + bi;
            u64 sb_start = new_start + sb * SIEVE_BLOCK_SIZE;
            u64 sb_end = std::min(sb_start + SIEVE_BLOCK_SIZE - 1, new_end);

            W30Sieve seg;
            seg.init(sb_start, sb_end);
            sieve_segment_w30(seg);

            u64 first_cp = sb_start / CHECKPOINT_INTERVAL;
            u64 last_cp = sb_end / CHECKPOINT_INTERVAL;
            u64 n_cp = last_cp - first_cp + 1;

            local_counts[sb].resize(n_cp);

            u64 prev_count = 0;
            for (u64 c = 0; c < n_cp; c++) {
                u64 cp_idx = first_cp + c;
                u64 cp_end_num = std::min((cp_idx + 1) * CHECKPOINT_INTERVAL - 1, new_end);

                u64 end_group;
                if (cp_end_num >= seg.base) {
                    end_group = (cp_end_num - seg.base) / 30 + 1;
                    if (end_group > seg.num_groups) end_group = seg.num_groups;
                } else {
                    end_group = 0;
                }

                u64 cum = seg.count_primes_groups(end_group);
                local_counts[sb][c] = cum - prev_count;
                prev_count = cum;
            }
        }

        auto t_now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t_now - t_last_report).count();
        if (elapsed >= 5.0 || batch_hi == num_sieve_blocks) {
            double total_elapsed = std::chrono::duration<double>(t_now - t_start).count();
            double progress = (double)batch_hi / num_sieve_blocks * 100.0;
            double rate = (double)std::min(batch_hi * SIEVE_BLOCK_SIZE, new_range)
                          / total_elapsed / 1e6;
            double frac = (double)batch_hi / num_sieve_blocks;
            double eta = (frac > 0.001) ? total_elapsed * (1.0 - frac) / frac : 0;

            printf("\r  [%.1f%%] Block %llu/%llu | %.1f M/s | elapsed %.0fs | ETA %.0fs   ",
                   progress,
                   (unsigned long long)batch_hi,
                   (unsigned long long)num_sieve_blocks,
                   rate, total_elapsed, eta);
            fflush(stdout);
            t_last_report = t_now;
        }
    }

    // Prefix-sum new checkpoints, continuing from last stored value
    printf("\n  Extending prefix sum... ");
    fflush(stdout);

    u64 running = table.checkpoints[old_checkpoints - 1];
    table.checkpoints.resize(total_checkpoints);

    u64 cp_idx = old_checkpoints;
    for (u64 sb = 0; sb < num_sieve_blocks && cp_idx < total_checkpoints; sb++) {
        for (u64 c = 0; c < local_counts[sb].size() && cp_idx < total_checkpoints; c++) {
            running += local_counts[sb][c];
            table.checkpoints[cp_idx] = running;
            cp_idx++;
        }
    }

    table.end = new_end;
    table.num_checkpoints = total_checkpoints;
    printf("done.\n");

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(t_end - t_start).count();

    printf("  Saving to %s... ", table_file);
    fflush(stdout);
    if (table.save(table_file))
        printf("done.\n");
    else
        printf("FAILED!\n");

    printf("\n");
    printf("════════════════════════════════════════\n");
    printf("  π(%llu) = %llu\n", (unsigned long long)new_end, (unsigned long long)running);
    printf("  Extended from π(%llu) = %llu\n",
           (unsigned long long)old_end,
           (unsigned long long)table.checkpoints[old_checkpoints - 1]);
    printf("  New primes found: %llu\n",
           (unsigned long long)(running - table.checkpoints[old_checkpoints - 1]));
    printf("  Extend time: %.3f s\n", total_time);
    printf("  Rate: %.1f M/s\n", (double)new_range / total_time / 1e6);
    printf("  Table: %llu checkpoints, %.1f MB\n",
           (unsigned long long)total_checkpoints,
           (double)(total_checkpoints * 8) / (1024.0 * 1024.0));
    printf("════════════════════════════════════════\n");
}

// ============================================================================
// Query
// ============================================================================

static void do_query(u64 x, const char* table_file) {
    PiTableFine table;
    if (!table.load(table_file)) {
        fprintf(stderr, "Error: cannot load table from %s\n", table_file);
        return;
    }
    if (x > table.end) {
        fprintf(stderr, "Error: x=%llu exceeds table range [0, %llu]\n",
                (unsigned long long)x, (unsigned long long)table.end);
        return;
    }

    u64 sqrt_x = (u64)sqrt((double)x) + 1;
    build_base_primes(sqrt_x);
    build_mod_tables();

    auto t0 = std::chrono::high_resolution_clock::now();
    u64 result = table.query(x);
    auto t1 = std::chrono::high_resolution_clock::now();
    double us = std::chrono::duration<double>(t1 - t0).count() * 1e6;

    printf("π(%llu) = %llu\n", (unsigned long long)x, (unsigned long long)result);
    printf("Query time: %.1f µs\n", us);
}

// ============================================================================
// Info
// ============================================================================

static void do_info(const char* table_file) {
    PiTableFine table;
    if (!table.load(table_file)) {
        fprintf(stderr, "Error: cannot load table from %s\n", table_file);
        return;
    }

    printf("Pi Table (fine): %s\n", table_file);
    printf("  Range:              [0, %llu]\n", (unsigned long long)table.end);
    printf("  Checkpoint interval: %llu numbers\n", (unsigned long long)CHECKPOINT_INTERVAL);
    printf("  Checkpoints:        %llu\n", (unsigned long long)table.num_checkpoints);
    printf("  Table size:         %.1f MB\n",
           (double)(table.num_checkpoints * 8) / (1024.0 * 1024.0));
    printf("  π(%llu) = %llu\n",
           (unsigned long long)table.end,
           (unsigned long long)table.checkpoints[table.num_checkpoints - 1]);

    // Sample queries
    printf("\n  Sample values:\n");
    u64 sqrt_end = (u64)sqrt((double)table.end) + 1;
    build_base_primes(sqrt_end);
    build_mod_tables();

    u64 samples[] = {100, 1000, 1000000, 1000000000, 1000000000000ULL};
    for (auto s : samples) {
        if (s <= table.end) {
            u64 result = table.query(s);
            printf("    π(%llu) = %llu\n", (unsigned long long)s, (unsigned long long)result);
        }
    }
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s build  <end> [threads] [table_file]\n", argv[0]);
        fprintf(stderr, "  %s extend <new_end> [threads] [table_file]\n", argv[0]);
        fprintf(stderr, "  %s query  <x> [table_file]\n", argv[0]);
        fprintf(stderr, "  %s info   [table_file]\n", argv[0]);
        return 1;
    }

    const char* mode = argv[1];

    if (strcmp(mode, "build") == 0) {
        if (argc < 3) { fprintf(stderr, "Need <end>\n"); return 1; }
        u64 end = strtoull(argv[2], nullptr, 10);
        int threads = (argc >= 4) ? atoi(argv[3]) : 1;
        const char* tf = (argc >= 5) ? argv[4] : "pi_table.dat";
        do_build(end, threads, tf);

    } else if (strcmp(mode, "extend") == 0) {
        if (argc < 3) { fprintf(stderr, "Need <new_end>\n"); return 1; }
        u64 new_end = strtoull(argv[2], nullptr, 10);
        int threads = (argc >= 4) ? atoi(argv[3]) : 1;
        const char* tf = (argc >= 5) ? argv[4] : "pi_table.dat";
        do_extend(new_end, threads, tf);

    } else if (strcmp(mode, "query") == 0) {
        if (argc < 3) { fprintf(stderr, "Need <x>\n"); return 1; }
        u64 x = strtoull(argv[2], nullptr, 10);
        const char* tf = (argc >= 4) ? argv[3] : "pi_table.dat";
        do_query(x, tf);

    } else if (strcmp(mode, "info") == 0) {
        const char* tf = (argc >= 3) ? argv[2] : "pi_table.dat";
        do_info(tf);

    } else {
        fprintf(stderr, "Unknown mode: %s\n", mode);
        return 1;
    }

    return 0;
}