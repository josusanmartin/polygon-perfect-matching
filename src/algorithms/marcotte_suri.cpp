// Implementation of the Find_Matching procedure of Marcotte & Suri (1991),
// "Fast Matching Algorithms for Points on a Polygon", SIAM J. Comput. 20(3),
// pp. 405-422. The structure of this file mirrors the procedure in Figure 4
// (page 417). Inline comments tagged with `Paper line N` reference the
// numbered lines of that figure; comments tagged with `Lemma K` reference the
// numbered lemmas of the paper.

#include "mwpm/algorithms.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace mwpm {

namespace {

// Result of a recursive call to find_matching_rec. `extensible` collects all
// edges that have been placed into the matching via Lemma 5 breakthroughs at
// this level or below. `residual` contains the original-input indices of the
// points that survived the conquer phase (i.e. those that would be closed off
// by the line-24 boundary matching at the top level), in CCW order.
struct RecResult {
    std::vector<Pair> extensible;
    std::vector<std::size_t> residual;
};

// Compute the result of a single conquer step: given two CCW-ordered residual
// halves R (from P_1) and L (from P_2), assign weights, build G_1 and G_2,
// process the negative edges in their total order, and return any extensible
// edges found together with the points that survived. When `simple_scan` is
// true, G_1 and G_2 are built via direct O(|R| * |L|) weighted-distance
// comparisons; otherwise they are built via the linear-time SMAWK matrix-
// searching algorithm.
RecResult conquer_phase(const std::vector<Point>& points,
                        const std::vector<std::size_t>& R,
                        const std::vector<std::size_t>& L,
                        bool simple_scan);

// -------------------------------------------------------------------------
// SMAWK row-minima for a totally monotone matrix.
//
// `solve_smawk(rows, cols, A)` returns, for each row index in `rows`, the
// column index from `cols` that achieves the row minimum. Both `rows` and
// `cols` must be sorted in increasing order. `A(row_idx, col_idx)` evaluates
// the (row, col) entry of the underlying matrix on demand. Ties are broken
// in favour of the leftmost column (matching the strict-< comparison used by
// the simple-scan reference path).
//
// Reference: Aggarwal, Klawe, Moran, Shor, Wilber, "Geometric applications of
// a matrix-searching algorithm", Algorithmica 2 (1987), pp. 209-233.
// -------------------------------------------------------------------------

template <typename F>
std::vector<int> smawk_reduce(const std::vector<int>& rows,
                              const std::vector<int>& cols,
                              F&& A) {
    // REDUCE: prune `cols` down to a subset of size <= |rows| containing
    // every row-minimising column. The stack holds candidate columns in the
    // order they were accepted; at any point, if `stack` has k entries then
    // entry k - 1 beats every later-popped column at row rows[k - 1].
    //
    // We pop on strict-> ("top is strictly worse than the new column at the
    // diagnostic row") so that ties keep the earlier column on the stack,
    // matching the leftmost-min tie-break used downstream.
    std::vector<int> stack;
    stack.reserve(rows.size());
    for (int c : cols) {
        while (!stack.empty()) {
            const std::size_t k = stack.size() - 1;
            if (A(rows[k], stack.back()) > A(rows[k], c)) {
                stack.pop_back();
            } else {
                break;
            }
        }
        if (stack.size() < rows.size()) {
            stack.push_back(c);
        }
        // Otherwise we already have |rows| candidates and no later column
        // could displace any of them, so c is dropped.
    }
    return stack;
}

template <typename F>
std::vector<int> solve_smawk(const std::vector<int>& rows,
                             const std::vector<int>& cols,
                             F&& A) {
    if (rows.empty()) return {};

    // Step 1: REDUCE columns when there are more columns than rows.
    const std::vector<int> Cp =
        cols.size() > rows.size() ? smawk_reduce(rows, cols, A) : cols;

    // Step 2: recurse on the rows at odd 0-based positions (rows[1], rows[3],
    // ...). Their row-minima are computed with the reduced column set Cp.
    std::vector<int> rows_odd;
    rows_odd.reserve(rows.size() / 2);
    for (std::size_t k = 1; k < rows.size(); k += 2) {
        rows_odd.push_back(rows[k]);
    }
    const std::vector<int> M_odd = solve_smawk(rows_odd, Cp, A);

    // Step 3: INTERPOLATE for the rows at even positions. Their minimum
    // columns are bracketed by the minimum columns of the neighbouring odd
    // rows (by total monotonicity). We walk through Cp once with a single
    // pointer to find these bracket positions in amortised O(|Cp|).
    std::vector<int> result(rows.size());
    int cp_pos = 0;  // running position in Cp
    for (std::size_t k = 0; k < rows.size(); ++k) {
        if (k % 2 == 1) {
            // Odd row already solved.
            result[k] = M_odd[k / 2];
            // Advance cp_pos to the position of M_odd[k/2] in Cp.
            while (cp_pos < static_cast<int>(Cp.size()) &&
                   Cp[cp_pos] != M_odd[k / 2]) {
                ++cp_pos;
            }
            assert(cp_pos < static_cast<int>(Cp.size()));
        } else {
            // Even row: search Cp[lo..hi] inclusive.
            const int lo = (k > 0) ? cp_pos : 0;
            int hi;
            if (k + 1 < rows.size()) {
                // Find the position of M_odd[k/2] starting from lo. The
                // bracket positions are non-decreasing across rows so this
                // forward scan is amortised O(|Cp|) total.
                int t = lo;
                const int target = M_odd[k / 2];
                while (t < static_cast<int>(Cp.size()) && Cp[t] != target) {
                    ++t;
                }
                assert(t < static_cast<int>(Cp.size()));
                hi = t;
            } else {
                hi = static_cast<int>(Cp.size()) - 1;
            }
            assert(lo <= hi);
            int best_col = Cp[lo];
            double best_val = A(rows[k], Cp[lo]);
            for (int idx = lo + 1; idx <= hi; ++idx) {
                const double v = A(rows[k], Cp[idx]);
                if (v < best_val) {
                    best_val = v;
                    best_col = Cp[idx];
                }
            }
            result[k] = best_col;
        }
    }
    return result;
}

RecResult find_matching_rec(const std::vector<Point>& points,
                            const std::vector<std::size_t>& indices,
                            bool simple_scan) {
    RecResult out;
    const std::size_t N = indices.size();
    // The recursion always receives an even-sized contiguous run.
    assert(N % 2 == 0);
    if (N == 0) {
        return out;
    }
    if (N == 2) {
        // Base case: the two points become R or L for the parent's conquer
        // step; they will be paired by line 24 if no critical edge consumes
        // them earlier.
        out.residual = indices;
        return out;
    }

    // Paper line 1: split P into nearly-equal halves with even counts. We use
    // n_pairs = N/2 (each pair = (x_i, y_i)), keep ceil(n/2) pairs in P1 so
    // that for n >= 2 both halves are non-empty.
    const std::size_t n_pairs = N / 2;
    std::size_t p1_pairs = (n_pairs + 1) / 2;
    if (p1_pairs == 0) p1_pairs = 1;
    if (p1_pairs == n_pairs) p1_pairs = n_pairs - 1;  // ensure P_2 non-empty
    const std::size_t split = 2 * p1_pairs;
    std::vector<std::size_t> P1(indices.begin(), indices.begin() + split);
    std::vector<std::size_t> P2(indices.begin() + split, indices.end());

    // Paper line 2: solve the two halves recursively.
    RecResult r1 = find_matching_rec(points, P1, simple_scan);
    RecResult r2 = find_matching_rec(points, P2, simple_scan);

    // Paper line 3: M starts as the union of every extensible-set edge found
    // below this level.
    out.extensible = std::move(r1.extensible);
    out.extensible.insert(out.extensible.end(),
                          r2.extensible.begin(),
                          r2.extensible.end());

    // Paper line 4: Q is what remains after deleting endpoints of M; here it
    // is exactly the union of the two children's residuals. R is from P_1, L
    // from P_2, both in CCW order.
    RecResult merged = conquer_phase(points, r1.residual, r2.residual,
                                     simple_scan);

    out.extensible.insert(out.extensible.end(),
                          merged.extensible.begin(),
                          merged.extensible.end());
    out.residual = std::move(merged.residual);
    return out;
}

// ---------------------------------------------------------------------------
// Conquer phase: this is the body of Find_Matching from line 5 to line 23 of
// Figure 4. We model the combined sequence R then L as a single CCW labelling
//
//   x_0 = R[0], y_0 = R[1], x_1 = R[2], ..., y_{m-1} = R[2m-1],
//   x_m = L[0], y_m = L[1], ..., y_{r-1} = L[|L|-1]
//
// where m = |R|/2 and r = (|R| + |L|)/2. Indices i, j etc. below refer to
// this re-labelled space (so x_i sits at "position" 2*i and y_i at 2*i+1).
// ---------------------------------------------------------------------------

RecResult conquer_phase(const std::vector<Point>& points,
                        const std::vector<std::size_t>& R,
                        const std::vector<std::size_t>& L,
                        bool simple_scan) {
    RecResult out;
    const std::size_t N = R.size() + L.size();
    assert(R.size() % 2 == 0 && L.size() % 2 == 0);
    if (N == 0) return out;
    if (N == 2) {
        // Two-point conquer is degenerate: just hand the pair up to the
        // parent untouched.
        out.residual.reserve(2);
        out.residual.insert(out.residual.end(), R.begin(), R.end());
        out.residual.insert(out.residual.end(), L.begin(), L.end());
        return out;
    }

    const int m = static_cast<int>(R.size() / 2);
    const int r = static_cast<int>(N / 2);

    // Maps a position in the combined sequence (0..N-1) to the original index
    // in `points`.
    auto orig = [&](int pos) -> std::size_t {
        return pos < static_cast<int>(R.size()) ? R[pos]
                                                : L[pos - static_cast<int>(R.size())];
    };
    auto X = [&](int i) -> const Point& { return points[orig(2 * i)]; };
    auto Y = [&](int i) -> const Point& { return points[orig(2 * i + 1)]; };

    // Paper §3, eq. (4): vertex weights u_i, v_i. By construction the only
    // adjacent boundary pair whose weights need not sum to their distance is
    // (x_0, y_{r-1}); every other adjacency satisfies u + v = d.
    std::vector<double> u(static_cast<std::size_t>(r), 0.0);
    std::vector<double> v(static_cast<std::size_t>(r), 0.0);
    u[0] = 0.0;
    for (int i = 0; i < r; ++i) {
        v[i] = distance(X(i), Y(i)) - u[i];
        if (i + 1 < r) {
            u[i + 1] = distance(Y(i), X(i + 1)) - v[i];
        }
    }
    // Original-weights weighted distance D(x_i, y_j) = d(x_i, y_j) - u_i - v_j.
    auto D_orig = [&](int i, int j) -> double {
        return distance(X(i), Y(j)) - u[i] - v[j];
    };

    // Paper line 5: compute G_1 and G_2.
    //
    // G_1: for each x_i in R (i in [0, m)), find a y_j in L (j in [m, r))
    //      minimising D(x_i, y_j).
    // G_2: symmetrically, for each x_p in L (p in [m, r)), find a y_q in R.
    //
    // Two implementations:
    //   * simple_scan == true:  O(|R| * |L|) direct comparisons.
    //   * simple_scan == false: SMAWK row-minima on the totally-monotone
    //                           weighted-distance matrices, O(|R| + |L|).
    //
    // For SMAWK to apply with the standard (row-minima-shift-right) form of
    // total monotonicity, we orient each matrix top-to-bottom: G_1 walks R
    // from top (x_{m-1}) to bottom (x_0); G_2 walks R columns from top
    // (y_{m-1}) to bottom (y_0). The wrappers below translate (row, col)
    // SMAWK indices back into the canonical x_i / y_j indices.
    std::vector<int> g1_partner(static_cast<std::size_t>(m), -1);
    std::vector<int> g2_partner(static_cast<std::size_t>(r - m), -1);

    if (simple_scan) {
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
    } else {
        // G_1 via SMAWK. Rows index R top-to-bottom: SMAWK row k <-> x_{m-1-k}.
        // Cols index L top-to-bottom (already the natural order).
        if (m > 0 && r - m > 0) {
            std::vector<int> rows1(static_cast<std::size_t>(m));
            std::vector<int> cols1(static_cast<std::size_t>(r - m));
            std::iota(rows1.begin(), rows1.end(), 0);
            std::iota(cols1.begin(), cols1.end(), 0);
            auto A_g1 = [&](int row, int col) -> double {
                const int x_idx = m - 1 - row;
                const int y_idx = m + col;
                return D_orig(x_idx, y_idx);
            };
            const std::vector<int> mins1 = solve_smawk(rows1, cols1, A_g1);
            for (std::size_t k = 0; k < mins1.size(); ++k) {
                const int x_idx = m - 1 - static_cast<int>(k);
                const int y_idx = m + mins1[k];
                g1_partner[static_cast<std::size_t>(x_idx)] = y_idx;
            }
        }
        // G_2 via SMAWK. Rows index L top-to-bottom (natural). Cols index R
        // top-to-bottom: SMAWK col l <-> y_{m-1-l} (R reversed).
        if (r - m > 0 && m > 0) {
            std::vector<int> rows2(static_cast<std::size_t>(r - m));
            std::vector<int> cols2(static_cast<std::size_t>(m));
            std::iota(rows2.begin(), rows2.end(), 0);
            std::iota(cols2.begin(), cols2.end(), 0);
            auto A_g2 = [&](int row, int col) -> double {
                const int x_idx = m + row;
                const int y_idx = m - 1 - col;
                return D_orig(x_idx, y_idx);
            };
            const std::vector<int> mins2 = solve_smawk(rows2, cols2, A_g2);
            for (std::size_t k = 0; k < mins2.size(); ++k) {
                const int y_idx = m - 1 - mins2[k];
                g2_partner[k] = y_idx;
            }
        }
    }

    // Paper line 6: List_1 (resp. List_2) is the list of G_1 (resp. G_2) edges
    // in their total order from the y_{m-1} x_m interface ("top") downwards.
    // For G_1, "higher" means larger x-index in R; for G_2 it means smaller
    // x-index in L. We expose these as front-pointers `p1`, `p2` together
    // with the underlying lists list1[], list2[].
    std::vector<int> list1;
    list1.reserve(static_cast<std::size_t>(m));
    for (int i = m - 1; i >= 0; --i) list1.push_back(i);
    std::vector<int> list2;
    list2.reserve(static_cast<std::size_t>(r - m));
    for (int i = m; i < r; ++i) list2.push_back(i);

    // `alive[pos]` tracks whether the point at combined position `pos` still
    // belongs to Q. Lines 17 and 21 logically remove edges with deleted
    // endpoints; we implement the removal lazily by skipping over dead edges
    // when advancing the front pointers.
    std::vector<char> alive(static_cast<std::size_t>(N), 1);

    auto g1_partner_of = [&](int i) -> int {
        return g1_partner[static_cast<std::size_t>(i)];
    };
    auto g2_partner_of = [&](int i) -> int {
        return g2_partner[static_cast<std::size_t>(i - m)];
    };

    auto edge_alive_g1 = [&](int i) -> bool {
        const int j = g1_partner_of(i);
        if (j < 0) return false;
        return alive[2 * i] && alive[2 * j + 1];
    };
    auto edge_alive_g2 = [&](int i) -> bool {
        const int j = g2_partner_of(i);
        if (j < 0) return false;
        return alive[2 * i] && alive[2 * j + 1];
    };

    int p1 = 0, p2 = 0;
    auto advance_g1 = [&]() {
        while (p1 < static_cast<int>(list1.size()) && !edge_alive_g1(list1[p1])) ++p1;
    };
    auto advance_g2 = [&]() {
        while (p2 < static_cast<int>(list2.size()) && !edge_alive_g2(list2[p2])) ++p2;
    };

    // Two convex chords on the boundary cross iff their endpoints interleave
    // cyclically. With combined positions for x_i (=2i, R), y_j (=2j+1, L),
    // x_p (=2p, L), y_q (=2q+1, R), this works out to:
    auto chords_cross = [&](int i, int j, int p, int q) -> bool {
        const bool d_on_ccw_arc = (q >= i);
        const bool c_on_ccw_arc = (p <= j);
        return d_on_ccw_arc != c_on_ccw_arc;
    };
    // For non-crossing chords, the G_1 edge precedes (is closer to the top)
    // exactly when both endpoints of the G_2 edge lie "below" it: q < i and
    // p > j. (The two conditions are equivalent under non-crossing.)
    auto g1_precedes_g2 = [&](int i, int j, int p, int q) -> bool {
        (void)i; (void)q;
        return p > j;
    };

    double delta = 0.0;  // Paper line 7.

    // Paper lines 8-23: main loop. The loop invariant is "the union of M and
    // an optimal matching of Q is an optimal matching of P". At each iteration
    // we examine the smallest (highest) negative edge among List_1 ∪ List_2.
    while (true) {
        advance_g1();
        advance_g2();
        const bool list1_empty = p1 >= static_cast<int>(list1.size());
        const bool list2_empty = p2 >= static_cast<int>(list2.size());
        if (list1_empty && list2_empty) break;  // Paper line 23.

        const int g1_i = list1_empty ? -1 : list1[p1];
        const int g1_j = list1_empty ? -1 : g1_partner_of(g1_i);
        const int g2_p = list2_empty ? -1 : list2[p2];
        const int g2_q = list2_empty ? -1 : g2_partner_of(g2_p);

        bool found_type1 = false;
        bool found_type2 = false;
        int crit1_i = -1, crit1_j = -1;  // For G_1 breakthrough.
        int crit2_p = -1, crit2_q = -1;  // For G_2 breakthrough.

        // Paper line 9: do the candidate G_1 / G_2 edges intersect?
        const bool intersect = !list1_empty && !list2_empty &&
                               chords_cross(g1_i, g1_j, g2_p, g2_q);

        if (intersect) {
            // Paper line 10: when the candidate edges cross, both candidates
            // are simultaneously eligible for evaluation against the current
            // delta-shifted weighted distance. The paper guarantees at most
            // one of them will end up being the critical edge.
            found_type1 = (D_orig(g1_i, g1_j) + delta < 0.0);
            found_type2 = (D_orig(g2_p, g2_q) - delta < 0.0);
            crit1_i = g1_i; crit1_j = g1_j;
            crit2_p = g2_p; crit2_q = g2_q;
            // Paper line 11: remove both consumed edges from their lists.
            ++p1;
            ++p2;
        } else {
            // Paper line 12: scan whichever of List_1, List_2 yields the
            // smaller (higher) edge in the total order. If only one list has
            // surviving edges, take it; otherwise compare against
            // g1_precedes_g2.
            bool process_g1;
            if (list2_empty) {
                process_g1 = true;
            } else if (list1_empty) {
                process_g1 = false;
            } else {
                process_g1 = g1_precedes_g2(g1_i, g1_j, g2_p, g2_q);
            }

            if (process_g1) {
                // Paper line 13.
                found_type1 = (D_orig(g1_i, g1_j) + delta < 0.0);
                crit1_i = g1_i;
                crit1_j = g1_j;
                ++p1;
            } else {
                // Paper line 14.
                found_type2 = (D_orig(g2_p, g2_q) - delta < 0.0);
                crit2_p = g2_p;
                crit2_q = g2_q;
                ++p2;
            }
        }

        // Paper line 15-18: G_1 breakthrough (Lemma 5 with k < l). The
        // critical edge is x_i y_j with i < j; the literal extensible set is
        // { y_i x_{i+1}, y_{i+1} x_{i+2}, ..., y_{j-1} x_j }. Lemma 5 is
        // stated in terms of the *current* Q, however, so when prior
        // breakthroughs have already removed points strictly between x_i and
        // y_j we must pair up whichever points remain alive in CCW order
        // rather than blindly applying the index pattern. The same Q-relative
        // reading holds for the G_2 branch below.
        if (found_type1) {
            const int i = crit1_i;
            const int j = crit1_j;
            assert(i < j);
            std::vector<int> survivors_between;
            survivors_between.reserve(static_cast<std::size_t>(2 * (j - i)));
            for (int pos = 2 * i + 1; pos <= 2 * j; ++pos) {
                if (alive[static_cast<std::size_t>(pos)]) survivors_between.push_back(pos);
            }
            assert(survivors_between.size() % 2 == 0);
            for (std::size_t s = 0; s + 1 < survivors_between.size(); s += 2) {
                const int a = survivors_between[s];
                const int b = survivors_between[s + 1];
                out.extensible.emplace_back(orig(a), orig(b));
                alive[static_cast<std::size_t>(a)] = 0;
                alive[static_cast<std::size_t>(b)] = 0;
            }
            // Paper line 17: edges x_r y_s in List_2 with r <= j or s >= i are
            // exactly those whose endpoints were just deleted. We deal with
            // this lazily inside advance_g2().
            // Paper line 18: delta := u_i + v_j - d(x_i, y_j). With the
            // accumulated-shift identity Delta = -D_orig(x_i, y_j) (see the
            // "update lemma" derivation), the new delta evaluated with the
            // *original* weights is exactly this quantity.
            delta = u[static_cast<std::size_t>(i)] +
                    v[static_cast<std::size_t>(j)] -
                    distance(X(i), Y(j));
        }

        // Paper line 19-22: G_2 breakthrough (Lemma 5 with k > l + 1). The
        // critical edge is x_p y_q with p > q + 1; the literal extensible set
        // is { x_{q+1} y_{q+1}, ..., x_{p-1} y_{p-1} }. As in the G_1 case
        // above we pair up the alive survivors strictly between y_q and x_p,
        // since prior breakthroughs may have already consumed some pairs.
        if (found_type2) {
            const int p = crit2_p;
            const int q = crit2_q;
            assert(p > q + 1);
            std::vector<int> survivors_between;
            survivors_between.reserve(static_cast<std::size_t>(2 * (p - q)));
            for (int pos = 2 * q + 2; pos <= 2 * p - 1; ++pos) {
                if (alive[static_cast<std::size_t>(pos)]) survivors_between.push_back(pos);
            }
            assert(survivors_between.size() % 2 == 0);
            for (std::size_t s = 0; s + 1 < survivors_between.size(); s += 2) {
                const int a = survivors_between[s];
                const int b = survivors_between[s + 1];
                out.extensible.emplace_back(orig(a), orig(b));
                alive[static_cast<std::size_t>(a)] = 0;
                alive[static_cast<std::size_t>(b)] = 0;
            }
            // Paper line 21: lazy removal handled by advance_g1().
            // Paper line 22: delta := -u_p - v_q + d(x_p, y_q). This makes
            // delta negative (D_orig(x_p, y_q) < 0 since the edge is
            // critical), as expected by the §4 description.
            delta = -u[static_cast<std::size_t>(p)] -
                    v[static_cast<std::size_t>(q)] +
                    distance(X(p), Y(q));
        }
    }

    // Build the residual: the points still alive, kept in CCW order. The
    // parent (or the top-level wrapper, if this is the outermost call) will
    // close them off via line 24.
    out.residual.reserve(static_cast<std::size_t>(N));
    for (int pos = 0; pos < static_cast<int>(N); ++pos) {
        if (alive[static_cast<std::size_t>(pos)]) {
            out.residual.push_back(orig(pos));
        }
    }
    return out;
}

}  // namespace

MatchingResult marcotte_suri_matching(const std::vector<Point>& points,
                                      bool simple_scan) {
    if (points.size() % 2 != 0) {
        throw std::invalid_argument("marcotte_suri_matching requires an even number of points");
    }
    MatchingResult result;
    if (points.empty()) return result;

    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});

    RecResult rec = find_matching_rec(points, indices, simple_scan);

    // Paper line 24 at the top level: pair the surviving residual points as
    // the boundary matching M_1 = {x_0 y_0, x_1 y_1, ...}.
    auto& res = rec.residual;
    assert(res.size() % 2 == 0);
    for (std::size_t i = 0; i + 1 < res.size(); i += 2) {
        rec.extensible.emplace_back(res[i], res[i + 1]);
    }

    result.pairs = std::move(rec.extensible);
    result.cost = 0.0;
    for (const auto& [a, b] : result.pairs) {
        result.cost += distance(points[a], points[b]);
    }
    return result;
}

}  // namespace mwpm
