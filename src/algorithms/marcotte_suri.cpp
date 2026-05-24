// Implementation of the Find_Matching procedure of Marcotte & Suri (1991),
// "Fast Matching Algorithms for Points on a Polygon", SIAM J. Comput. 20(3),
// pp. 405-422. The structure of this file mirrors the procedure in Figure 4
// (page 417). Inline comments tagged with `Paper line N` reference the
// numbered lines of that figure; comments tagged with `Lemma K` reference the
// numbered lemmas of the paper.
//
// Performance notes (vs. the textbook version of this implementation):
//
//   * All transient buffers used by find_matching_rec, the conquer phase and
//     SMAWK are owned by a single `Scratch` struct allocated once at the top
//     of marcotte_suri_matching and threaded through the recursion. Recursion
//     levels and SMAWK levels each reuse pre-grown vectors via .clear() /
//     .assign() / .resize() rather than reallocating. This eliminates the
//     per-level operator new / free traffic that dominated the profile of a
//     straightforward implementation.
//   * The conquer phase builds a single flat `combined_pts[]` array of points
//     in the canonical CCW labelling x_0, y_0, x_1, y_1, ..., y_{r-1}, so the
//     inner D_orig() lookups become two array accesses + two scalar adds (no
//     R vs. L branch and no double indirection through the original index
//     vector).
//   * find_matching_rec passes views into the original `indices` array down
//     the recursion (just a pointer + length pair); P1 / P2 are never copied.
//     Residual / extensible buffers for both children at a given depth live
//     in a per-depth scratch pool, so each level reuses its own pair of
//     buffers across calls.
//   * The conquer main loop tracks `p1` / `p2` directly as running indices
//     descending / ascending through R and L; no `list1` / `list2` vectors
//     are needed.

#include "mwpm/algorithms.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mwpm {

namespace {

// =========================================================================
// Scratch buffers
// =========================================================================
//
// One Scratch instance is created at the top of marcotte_suri_matching and
// passed by reference through the entire computation. Every transient vector
// inside the conquer phase and SMAWK lives here. Within one call, buffers
// are reused freely; nothing is freed until the Scratch goes out of scope.

struct Scratch {
    // --- conquer phase ---
    // CCW-merged inputs: combined_orig[k] is the original-input index of the
    // point whose canonical conquer-phase label is k (0..N-1). combined_pts[k]
    // is the corresponding mwpm::Point. Built once per conquer call.
    std::vector<std::size_t> combined_orig;
    std::vector<Point> combined_pts;

    std::vector<double> u_buf;       // size r = N/2
    std::vector<double> v_buf;       // size r = N/2
    std::vector<int> g1_partner;     // size m
    std::vector<int> g2_partner;     // size r - m
    std::vector<std::uint8_t> alive; // size N
    std::vector<int> survivors;      // scratch for the breakthrough loops

    // --- SMAWK ---
    // The recursion in solve_smawk_impl uses one (rows_odd, cols_reduced,
    // M_odd) triple per recursion depth. Depth never exceeds log2(N), so we
    // grow these vector-of-vectors lazily.
    std::vector<std::vector<int>> smawk_rows_odd;
    std::vector<std::vector<int>> smawk_cols_reduced;
    std::vector<std::vector<int>> smawk_M_odd;
    std::vector<int> smawk_reduce_stack;
    // Caller-supplied initial rows / cols and final result for the top-level
    // SMAWK call inside each conquer phase. Kept separate from the per-level
    // pools so the level-0 invocation does not alias its own rows_odd buffer.
    std::vector<int> smawk_init_rows;
    std::vector<int> smawk_init_cols;
    std::vector<int> smawk_init_result;

    // --- find_matching_rec ---
    // Per-depth residual / extensible buffers for the two children of a
    // recursion node. residual_pool[2 * depth + child] is the residual that
    // the child at the given depth will write into; same for extensible_pool.
    std::vector<std::vector<std::size_t>> residual_pool;
    std::vector<std::vector<Pair>> extensible_pool;

    // Both pools are sized once at the top of marcotte_suri_matching. Once
    // sized, references / pointers taken into the outer vectors stay stable
    // for the duration of the call. Growing mid-recursion would invalidate
    // every such reference, so we never do that.
    void ensure_smawk_depth(std::size_t d) {
        if (smawk_rows_odd.size() < d) {
            smawk_rows_odd.resize(d);
            smawk_cols_reduced.resize(d);
            smawk_M_odd.resize(d);
        }
    }

    void ensure_recursion_depth(std::size_t levels) {
        const std::size_t need = 2 * levels;
        if (residual_pool.size() < need) {
            residual_pool.resize(need);
            extensible_pool.resize(need);
        }
    }
};

// =========================================================================
// SMAWK row-minima for a totally monotone matrix
// =========================================================================
//
// `solve_smawk(rows, cols, A, out)` writes, for each row index in `rows`, the
// column index from `cols` that achieves the row minimum. Both `rows` and
// `cols` must be sorted in increasing order. `A(row_idx, col_idx)` evaluates
// the (row, col) entry on demand. Ties are broken in favour of the leftmost
// column (matching the strict-< comparison used by the simple-scan reference
// path so the two MS variants produce the same matching).
//
// Reference: Aggarwal, Klawe, Moran, Shor, Wilber, "Geometric applications of
// a matrix-searching algorithm", Algorithmica 2 (1987), pp. 209-233.

template <typename F>
void smawk_reduce_into(const std::vector<int>& rows,
                       const std::vector<int>& cols,
                       F& A,
                       std::vector<int>& stack,
                       std::vector<int>& out) {
    // REDUCE: keep at most |rows| candidate columns. The stack holds the
    // current candidates in the order they were accepted; at any point, if
    // `stack` has k entries then entry k - 1 beats every later-popped column
    // at row rows[k - 1].
    //
    // We pop on strict-> ("top is strictly worse than the new column at the
    // diagnostic row") so that ties keep the earlier column on the stack,
    // matching the leftmost-min tie-break used downstream.
    stack.clear();
    const std::size_t cap = rows.size();
    stack.reserve(cap);
    for (const int c : cols) {
        while (!stack.empty()) {
            const std::size_t k = stack.size() - 1;
            if (A(rows[k], stack.back()) > A(rows[k], c)) {
                stack.pop_back();
            } else {
                break;
            }
        }
        if (stack.size() < cap) {
            stack.push_back(c);
        }
        // Otherwise we already have |rows| candidates and no later column
        // could displace any of them, so c is dropped.
    }
    out.assign(stack.begin(), stack.end());
}

template <typename F>
void solve_smawk_impl(Scratch& s,
                      std::size_t level,
                      const std::vector<int>& rows,
                      const std::vector<int>& cols,
                      F& A,
                      std::vector<int>& result_out) {
    if (rows.empty()) {
        result_out.clear();
        return;
    }
    // SMAWK pool is pre-grown at the top of marcotte_suri_matching.
    assert(level < s.smawk_rows_odd.size());

    // Step 1: REDUCE when there are more cols than rows. Cp is either the
    // freshly reduced column list (owned by this level's scratch) or just
    // the caller's `cols` if no reduction was needed.
    std::vector<int>& Cp_buf = s.smawk_cols_reduced[level];
    const std::vector<int>* Cp_ptr;
    if (cols.size() > rows.size()) {
        smawk_reduce_into(rows, cols, A, s.smawk_reduce_stack, Cp_buf);
        Cp_ptr = &Cp_buf;
    } else {
        Cp_ptr = &cols;
    }
    const std::vector<int>& Cp = *Cp_ptr;

    // Step 2: build odd-indexed rows and recurse. Each level uses its own
    // rows_odd buffer so the const reference `rows` (which may alias the
    // parent's rows_odd buffer) stays valid for the INTERPOLATE step.
    std::vector<int>& rows_odd = s.smawk_rows_odd[level];
    rows_odd.clear();
    rows_odd.reserve(rows.size() / 2);
    for (std::size_t k = 1; k < rows.size(); k += 2) {
        rows_odd.push_back(rows[k]);
    }

    std::vector<int>& M_odd = s.smawk_M_odd[level];
    solve_smawk_impl(s, level + 1, rows_odd, Cp, A, M_odd);

    // Step 3: INTERPOLATE for even-indexed rows. Their minimum columns are
    // bracketed in Cp by the neighbouring odd rows' minima, so a single
    // forward sweep with a running `cp_pos` finds the bracket in amortised
    // O(|Cp|) total work across all rows.
    result_out.assign(rows.size(), 0);
    int* const result_data = result_out.data();
    const int* const Cp_data = Cp.data();
    const int Cp_size = static_cast<int>(Cp.size());

    int cp_pos = 0;
    for (std::size_t k = 0; k < rows.size(); ++k) {
        if ((k & 1u) == 1u) {
            // Odd row already solved.
            const int m_col = M_odd[k >> 1];
            result_data[k] = m_col;
            while (cp_pos < Cp_size && Cp_data[cp_pos] != m_col) ++cp_pos;
            assert(cp_pos < Cp_size);
        } else {
            const int lo = (k > 0) ? cp_pos : 0;
            int hi;
            if (k + 1 < rows.size()) {
                int t = lo;
                const int target = M_odd[k >> 1];
                while (t < Cp_size && Cp_data[t] != target) ++t;
                assert(t < Cp_size);
                hi = t;
            } else {
                hi = Cp_size - 1;
            }
            assert(lo <= hi);
            int best_col = Cp_data[lo];
            double best_val = A(rows[k], Cp_data[lo]);
            for (int idx = lo + 1; idx <= hi; ++idx) {
                const double v = A(rows[k], Cp_data[idx]);
                if (v < best_val) {
                    best_val = v;
                    best_col = Cp_data[idx];
                }
            }
            result_data[k] = best_col;
        }
    }
}

// =========================================================================
// Conquer phase (lines 5..23 of Figure 4)
// =========================================================================
//
// We model the combined sequence R then L as a single CCW labelling
//
//   x_0 = R[0], y_0 = R[1], x_1 = R[2], ..., y_{m-1} = R[2m-1],
//   x_m = L[0], y_m = L[1], ..., y_{r-1} = L[|L|-1]
//
// where m = |R|/2 and r = (|R| + |L|)/2. After this point we index entirely
// in this re-labelled space; combined_orig and combined_pts give us O(1)
// access to the original-input index and the Point of any conquer-phase
// label without further branching.
//
// `out_extensible` is appended to (not cleared). `out_residual` is overwritten.

void run_conquer(Scratch& s,
                 const std::vector<Point>& points,
                 const std::size_t* R, std::size_t R_size,
                 const std::size_t* L, std::size_t L_size,
                 bool simple_scan,
                 std::vector<Pair>& out_extensible,
                 std::vector<std::size_t>& out_residual) {
    const std::size_t N = R_size + L_size;
    assert(R_size % 2 == 0 && L_size % 2 == 0);
    out_residual.clear();
    if (N == 0) return;
    if (N == 2) {
        // Degenerate two-point conquer: hand the pair up to the parent
        // untouched (line 24 of the outermost call will close it off).
        out_residual.reserve(2);
        out_residual.insert(out_residual.end(), R, R + R_size);
        out_residual.insert(out_residual.end(), L, L + L_size);
        return;
    }

    const int m = static_cast<int>(R_size / 2);
    const int r = static_cast<int>(N / 2);

    // ---- Build the flat combined arrays. ----
    auto& combined_orig = s.combined_orig;
    auto& combined_pts = s.combined_pts;
    combined_orig.resize(N);
    combined_pts.resize(N);
    std::memcpy(combined_orig.data(), R, R_size * sizeof(std::size_t));
    std::memcpy(combined_orig.data() + R_size, L, L_size * sizeof(std::size_t));
    for (std::size_t k = 0; k < N; ++k) {
        combined_pts[k] = points[combined_orig[k]];
    }

    const Point* const cp = combined_pts.data();
    auto X = [cp](int i) -> const Point& { return cp[2 * i]; };
    auto Y = [cp](int i) -> const Point& { return cp[2 * i + 1]; };

    // ---- u, v vertex weights (paper §3, eq. (4)). ----
    auto& u = s.u_buf;
    auto& v = s.v_buf;
    u.assign(static_cast<std::size_t>(r), 0.0);
    v.assign(static_cast<std::size_t>(r), 0.0);
    // u[0] = 0; v[i] = d(x_i, y_i) - u[i]; u[i+1] = d(y_i, x_{i+1}) - v[i].
    for (int i = 0; i < r; ++i) {
        v[i] = distance(X(i), Y(i)) - u[i];
        if (i + 1 < r) {
            u[i + 1] = distance(Y(i), X(i + 1)) - v[i];
        }
    }

    // Original-weights weighted distance D(x_i, y_j) = d(x_i, y_j) - u_i - v_j.
    const double* const u_data = u.data();
    const double* const v_data = v.data();
    auto D_orig = [cp, u_data, v_data](int i, int j) -> double {
        const Point& a = cp[2 * i];
        const Point& b = cp[2 * j + 1];
        const double dx = a.x - b.x;
        const double dy = a.y - b.y;
        return std::sqrt(dx * dx + dy * dy) - u_data[i] - v_data[j];
    };

    // ---- Paper line 5: compute G_1 and G_2. ----
    auto& g1_partner = s.g1_partner;
    auto& g2_partner = s.g2_partner;
    g1_partner.assign(static_cast<std::size_t>(m), -1);
    g2_partner.assign(static_cast<std::size_t>(r - m), -1);

    if (simple_scan) {
        // O(|R| * |L|) reference path.
        for (int i = 0; i < m; ++i) {
            double best = std::numeric_limits<double>::infinity();
            int best_j = -1;
            for (int j = m; j < r; ++j) {
                const double dij = D_orig(i, j);
                if (dij < best) {
                    best = dij;
                    best_j = j;
                }
            }
            g1_partner[static_cast<std::size_t>(i)] = best_j;
        }
        for (int i = m; i < r; ++i) {
            double best = std::numeric_limits<double>::infinity();
            int best_j = -1;
            for (int j = 0; j < m; ++j) {
                const double dij = D_orig(i, j);
                if (dij < best) {
                    best = dij;
                    best_j = j;
                }
            }
            g2_partner[static_cast<std::size_t>(i - m)] = best_j;
        }
    } else if (m > 0 && r - m > 0) {
        // SMAWK path. For each matrix we orient rows / columns top-to-bottom
        // so the weighted-distance matrix is totally monotone in the standard
        // (row-minima-shift-right) form.

        // G_1: rows index R top-to-bottom (SMAWK row k <-> x_{m-1-k}). Cols
        // index L top-to-bottom (already the natural order).
        auto& rows1 = s.smawk_init_rows;
        auto& cols1 = s.smawk_init_cols;
        rows1.resize(static_cast<std::size_t>(m));
        cols1.resize(static_cast<std::size_t>(r - m));
        std::iota(rows1.begin(), rows1.end(), 0);
        std::iota(cols1.begin(), cols1.end(), 0);
        const int m_local = m;
        auto A_g1 = [&D_orig, m_local](int row, int col) -> double {
            const int x_idx = m_local - 1 - row;
            const int y_idx = m_local + col;
            return D_orig(x_idx, y_idx);
        };
        auto& mins1 = s.smawk_init_result;
        solve_smawk_impl(s, 0, rows1, cols1, A_g1, mins1);
        for (std::size_t k = 0; k < mins1.size(); ++k) {
            const int x_idx = m - 1 - static_cast<int>(k);
            const int y_idx = m + mins1[k];
            g1_partner[static_cast<std::size_t>(x_idx)] = y_idx;
        }

        // G_2: rows index L top-to-bottom (natural). Cols index R top-to-
        // bottom (SMAWK col l <-> y_{m-1-l}).
        auto& rows2 = s.smawk_init_rows;
        auto& cols2 = s.smawk_init_cols;
        rows2.resize(static_cast<std::size_t>(r - m));
        cols2.resize(static_cast<std::size_t>(m));
        std::iota(rows2.begin(), rows2.end(), 0);
        std::iota(cols2.begin(), cols2.end(), 0);
        auto A_g2 = [&D_orig, m_local](int row, int col) -> double {
            const int x_idx = m_local + row;
            const int y_idx = m_local - 1 - col;
            return D_orig(x_idx, y_idx);
        };
        auto& mins2 = s.smawk_init_result;
        solve_smawk_impl(s, 0, rows2, cols2, A_g2, mins2);
        for (std::size_t k = 0; k < mins2.size(); ++k) {
            const int y_idx = m - 1 - mins2[k];
            g2_partner[k] = y_idx;
        }
    }

    // ---- alive[] tracks Q membership; we remove breakthrough endpoints
    // lazily by checking alive[] when advancing p1 / p2. ----
    auto& alive = s.alive;
    alive.assign(N, 1);

    // Combined positions for x_i, y_j (in canonical labelling): x_i is at
    // position 2 * i, y_j at position 2 * j + 1.
    auto x_alive = [&alive](int i) -> bool {
        return alive[static_cast<std::size_t>(2 * i)] != 0;
    };
    auto y_alive = [&alive](int j) -> bool {
        return alive[static_cast<std::size_t>(2 * j + 1)] != 0;
    };

    // ---- Paper line 6: List_1, List_2 enumerated as running pointers. ----
    //
    // List_1 walks G_1 edges from "highest" to "lowest", which in canonical
    // labelling means x_i for i descending from m - 1 down to 0. List_2 walks
    // G_2 edges from highest to lowest, which means x_p for p ascending from
    // m up to r - 1. `p1` and `p2` are the current x-indices into R and L
    // respectively; an x-index is "consumed" by incrementing the pointer.
    int p1 = m - 1;         // current G_1 x_i (in R; descending)
    int p2 = m;             // current G_2 x_p (in L; ascending)

    auto advance_g1 = [&]() {
        while (p1 >= 0) {
            const int j = g1_partner[static_cast<std::size_t>(p1)];
            if (j >= 0 && x_alive(p1) && y_alive(j)) break;
            --p1;
        }
    };
    auto advance_g2 = [&]() {
        while (p2 < r) {
            const int j = g2_partner[static_cast<std::size_t>(p2 - m)];
            if (j >= 0 && x_alive(p2) && y_alive(j)) break;
            ++p2;
        }
    };

    // Two convex chords on the boundary cross iff their endpoints interleave
    // cyclically. With combined positions x_i = 2 i (R), y_j = 2 j + 1 (L),
    // x_p = 2 p (L), y_q = 2 q + 1 (R), this becomes:
    auto chords_cross = [](int i, int j, int p, int q) -> bool {
        const bool d_on_ccw_arc = (q >= i);
        const bool c_on_ccw_arc = (p <= j);
        return d_on_ccw_arc != c_on_ccw_arc;
    };
    // For non-crossing chords, the G_1 edge precedes (is closer to the top)
    // exactly when both endpoints of the G_2 edge lie "below" the G_1 edge:
    // q < i and p > j. The two conditions are equivalent under non-crossing,
    // so a single comparison (g2_p > g1_j) is enough -- see the main loop.

    double delta = 0.0;  // Paper line 7.

    // ---- Paper lines 8..23: main loop. ----
    while (true) {
        advance_g1();
        advance_g2();
        const bool list1_empty = (p1 < 0);
        const bool list2_empty = (p2 >= r);
        if (list1_empty && list2_empty) break;  // Paper line 23.

        const int g1_i = list1_empty ? -1 : p1;
        const int g1_j = list1_empty ? -1 : g1_partner[static_cast<std::size_t>(p1)];
        const int g2_p = list2_empty ? -1 : p2;
        const int g2_q = list2_empty ? -1
                                     : g2_partner[static_cast<std::size_t>(p2 - m)];

        bool found_type1 = false;
        bool found_type2 = false;
        int crit1_i = -1, crit1_j = -1;
        int crit2_p = -1, crit2_q = -1;

        const bool intersect = !list1_empty && !list2_empty &&
                               chords_cross(g1_i, g1_j, g2_p, g2_q);

        if (intersect) {
            // Paper line 10: candidates cross. Both are simultaneously
            // eligible for evaluation against the delta-shifted weighted
            // distance. At most one will end up being a critical edge.
            found_type1 = (D_orig(g1_i, g1_j) + delta < 0.0);
            found_type2 = (D_orig(g2_p, g2_q) - delta < 0.0);
            crit1_i = g1_i; crit1_j = g1_j;
            crit2_p = g2_p; crit2_q = g2_q;
            --p1;
            ++p2;
        } else {
            // Paper line 12: take whichever list yields the higher edge.
            bool process_g1;
            if (list2_empty) {
                process_g1 = true;
            } else if (list1_empty) {
                process_g1 = false;
            } else {
                // For non-crossing chords G_1 precedes G_2 iff g2_p > g1_j.
                process_g1 = (g2_p > g1_j);
            }
            if (process_g1) {
                found_type1 = (D_orig(g1_i, g1_j) + delta < 0.0);
                crit1_i = g1_i;
                crit1_j = g1_j;
                --p1;
            } else {
                found_type2 = (D_orig(g2_p, g2_q) - delta < 0.0);
                crit2_p = g2_p;
                crit2_q = g2_q;
                ++p2;
            }
        }

        // ---- Paper lines 15-18: G_1 breakthrough (Lemma 5 with k < l). ----
        if (found_type1) {
            const int i = crit1_i;
            const int j = crit1_j;
            assert(i < j);
            // Pair up alive points strictly between x_i and y_j (inclusive
            // of y_i and x_j, the literal extensible-set endpoints).
            auto& sb = s.survivors;
            sb.clear();
            sb.reserve(static_cast<std::size_t>(2 * (j - i)));
            for (int pos = 2 * i + 1; pos <= 2 * j; ++pos) {
                if (alive[static_cast<std::size_t>(pos)]) sb.push_back(pos);
            }
            assert(sb.size() % 2 == 0);
            const std::size_t* const cog = combined_orig.data();
            for (std::size_t k = 0; k + 1 < sb.size(); k += 2) {
                const int a = sb[k];
                const int b = sb[k + 1];
                out_extensible.emplace_back(cog[a], cog[b]);
                alive[static_cast<std::size_t>(a)] = 0;
                alive[static_cast<std::size_t>(b)] = 0;
            }
            // Paper line 17: List_2 edges with deleted endpoints are skipped
            // lazily inside advance_g2(); nothing to do here.
            // Paper line 18: delta := u_i + v_j - d(x_i, y_j).
            delta = u[static_cast<std::size_t>(i)] +
                    v[static_cast<std::size_t>(j)] -
                    distance(X(i), Y(j));
        }

        // ---- Paper lines 19-22: G_2 breakthrough (Lemma 5 with k > l + 1). ----
        if (found_type2) {
            const int p = crit2_p;
            const int q = crit2_q;
            assert(p > q + 1);
            auto& sb = s.survivors;
            sb.clear();
            sb.reserve(static_cast<std::size_t>(2 * (p - q)));
            for (int pos = 2 * q + 2; pos <= 2 * p - 1; ++pos) {
                if (alive[static_cast<std::size_t>(pos)]) sb.push_back(pos);
            }
            assert(sb.size() % 2 == 0);
            const std::size_t* const cog = combined_orig.data();
            for (std::size_t k = 0; k + 1 < sb.size(); k += 2) {
                const int a = sb[k];
                const int b = sb[k + 1];
                out_extensible.emplace_back(cog[a], cog[b]);
                alive[static_cast<std::size_t>(a)] = 0;
                alive[static_cast<std::size_t>(b)] = 0;
            }
            // Paper line 22: delta := -u_p - v_q + d(x_p, y_q).
            delta = -u[static_cast<std::size_t>(p)] -
                    v[static_cast<std::size_t>(q)] +
                    distance(X(p), Y(q));
        }
    }

    // ---- Build residual: alive points in CCW order. ----
    out_residual.reserve(N);
    const std::size_t* const cog = combined_orig.data();
    for (std::size_t pos = 0; pos < N; ++pos) {
        if (alive[pos]) out_residual.push_back(cog[pos]);
    }
}

// =========================================================================
// find_matching_rec
// =========================================================================
//
// Operates on a CCW-contiguous slice of the original `indices` array via the
// (indices_begin, n) view. The descent never copies P1 / P2. Each level
// writes its residual and extensible into caller-provided buffers
// (`out_residual`, `out_extensible`); the buffers at depth+1 used to hold
// the two children's outputs come from the per-depth scratch pool.

void find_matching_rec(Scratch& s,
                       std::size_t depth,
                       const std::vector<Point>& points,
                       const std::size_t* indices_begin,
                       std::size_t n,
                       bool simple_scan,
                       std::vector<std::size_t>& out_residual,
                       std::vector<Pair>& out_extensible) {
    assert(n % 2 == 0);
    out_residual.clear();
    out_extensible.clear();
    if (n == 0) return;
    if (n == 2) {
        out_residual.assign(indices_begin, indices_begin + n);
        return;
    }

    // Paper line 1: split into nearly-equal even halves. We use n_pairs = n/2;
    // p1_pairs = ceil(n_pairs / 2), with clamping so both halves are non-empty.
    const std::size_t n_pairs = n / 2;
    std::size_t p1_pairs = (n_pairs + 1) / 2;
    if (p1_pairs == 0) p1_pairs = 1;
    if (p1_pairs == n_pairs) p1_pairs = n_pairs - 1;
    const std::size_t split = 2 * p1_pairs;

    // Pool is pre-grown at the top of marcotte_suri_matching; no growth here.
    const std::size_t child_slot = 2 * (depth + 1);
    assert(child_slot + 1 < s.residual_pool.size());
    auto& r1_residual   = s.residual_pool[child_slot + 0];
    auto& r2_residual   = s.residual_pool[child_slot + 1];
    auto& r1_extensible = s.extensible_pool[child_slot + 0];
    auto& r2_extensible = s.extensible_pool[child_slot + 1];

    // Paper line 2: recurse on the two halves.
    find_matching_rec(s, depth + 1, points,
                      indices_begin, split,
                      simple_scan,
                      r1_residual, r1_extensible);
    find_matching_rec(s, depth + 1, points,
                      indices_begin + split, n - split,
                      simple_scan,
                      r2_residual, r2_extensible);

    // Paper line 3: M starts as the union of every extensible-set edge
    // collected below this level. The conquer step (line 4 onwards) will
    // append further extensible edges into out_extensible.
    out_extensible.reserve(r1_extensible.size() + r2_extensible.size());
    out_extensible.insert(out_extensible.end(),
                          r1_extensible.begin(), r1_extensible.end());
    out_extensible.insert(out_extensible.end(),
                          r2_extensible.begin(), r2_extensible.end());

    // Paper line 4: conquer the union of the two residuals. R is from P_1,
    // L from P_2, both already in CCW order.
    run_conquer(s, points,
                r1_residual.data(), r1_residual.size(),
                r2_residual.data(), r2_residual.size(),
                simple_scan,
                out_extensible,
                out_residual);
}

}  // namespace

MatchingResult marcotte_suri_matching(const std::vector<Point>& points,
                                      bool simple_scan) {
    if (points.size() % 2 != 0) {
        throw std::invalid_argument("marcotte_suri_matching requires an even number of points");
    }
    MatchingResult result;
    if (points.empty()) return result;

    // The very top of the recursion: provide an `indices` view 0, 1, ..., N-1.
    // Scratch is thread_local so its inner buffer *capacities* survive
    // between calls on the same thread; this amortises the one-time pool
    // setup cost away from each individual call. The library remains
    // thread-safe (each thread sees its own Scratch).
    static thread_local Scratch scratch;
    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});

    // Pre-grow every per-depth pool to its worst-case length once, before any
    // references into these pools are taken. After this point, no buffer
    // grows the outer vector-of-vectors, so references / pointers into the
    // pool stay stable for the rest of the call.
    //
    // The matching recursion halves the slice length at each step, so its
    // depth is at most ceil(log2(N / 2)) + 1; padding to 64 levels handles
    // any input size up to 2^64. SMAWK has the same log2(N) bound on
    // recursion depth. Both pools are vectors of empty vectors, so the
    // reserved-but-unused entries cost only a few bytes each. With a
    // thread_local Scratch, this resize is a no-op after the first call on a
    // given thread.
    constexpr std::size_t kMaxDepth = 64;
    scratch.ensure_recursion_depth(kMaxDepth);
    scratch.ensure_smawk_depth(kMaxDepth);

    // Final output buffers live at the bottom of the per-depth pool (depth 0).
    auto& top_residual   = scratch.residual_pool[0];
    auto& top_extensible = scratch.extensible_pool[0];

    find_matching_rec(scratch, 0, points,
                      indices.data(), indices.size(),
                      simple_scan,
                      top_residual, top_extensible);

    // Paper line 24 at the top level: pair the surviving residual points as
    // the boundary matching M_1 = {x_0 y_0, x_1 y_1, ...}.
    assert(top_residual.size() % 2 == 0);
    result.pairs.reserve(top_extensible.size() + top_residual.size() / 2);
    result.pairs = std::move(top_extensible);
    for (std::size_t i = 0; i + 1 < top_residual.size(); i += 2) {
        result.pairs.emplace_back(top_residual[i], top_residual[i + 1]);
    }

    double cost = 0.0;
    for (const auto& [a, b] : result.pairs) {
        cost += distance(points[a], points[b]);
    }
    result.cost = cost;
    return result;
}

}  // namespace mwpm
