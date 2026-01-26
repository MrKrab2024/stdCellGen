#include "stdcell/folding.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <vector>

namespace stdcell {

static void compositions_rec(int rem, int max_part, std::vector<int>& cur, std::vector<std::vector<int>>& out) {
    if (rem == 0) { out.push_back(cur); return; }
    const int up = std::min(rem, max_part);
    for (int v = 1; v <= up; ++v) {
        cur.push_back(v);
        compositions_rec(rem - v, max_part, cur, out);
        cur.pop_back();
    }
}

static std::vector<std::vector<int>> compositions(int total, int max_part) {
    std::vector<std::vector<int>> out; if (total <= 0) return out; std::vector<int> cur; compositions_rec(total, max_part, cur, out); return out;
}

static void row_assign_rec(int k, int num_rows, int i, std::vector<int>& cur, std::vector<std::vector<int>>& out) {
    if (i == k) { out.push_back(cur); return; }
    for (int r = 0; r < num_rows; ++r) { cur.push_back(r); row_assign_rec(k, num_rows, i+1, cur, out); cur.pop_back(); }
}

static std::vector<std::vector<int>> row_assignments(int k, int num_rows) {
    std::vector<std::vector<int>> out; if (k <= 0 || num_rows <= 0) { if (k==0) out.push_back({}); return out; }
    std::vector<int> cur; row_assign_rec(k, num_rows, 0, cur, out); return out;
}

static std::vector<Orient> pattern_for_len(int count, bool start_sd) {
    std::vector<Orient> v; v.reserve(count);
    for (int i = 0; i < count; ++i) {
        bool even = (i % 2 == 0);
        bool sd = start_sd ? even : !even;
        v.push_back(sd ? Orient::SD : Orient::DS);
    }
    return v;
}

static std::vector<std::vector<Orient>> row_orient_patterns(int count) {
    if (count <= 0) return { { } };
    if (count == 1) return { {Orient::SD}, {Orient::DS} };
    return { pattern_for_len(count, true), pattern_for_len(count, false) };
}

static inline int phys_y_for(bool is_nmos, int logical_row) {
    if ((logical_row % 2) == 0) {
        return is_nmos ? (2*logical_row) : (2*logical_row + 1);
    } else {
        return is_nmos ? (2*logical_row + 1) : (2*logical_row);
    }
}

int effective_nfin(const Transistor& t, const TechRules& tr) {
    if (t.nfin > 0) return t.nfin;
    if (tr.default_w_per_fin_um > 0.0) {
        int nf = (int)std::max(1.0, std::round(t.w_um / tr.default_w_per_fin_um));
        return nf;
    }
    // Fallback: 1 fin if unknown
    return 1;
}

std::vector<PairVariant> enumerate_pair_variants(
    const Transistor& n,
    const Transistor& p,
    int max_fins_per_row,
    int num_logical_rows,
    const TechRules& tr
) {
    std::vector<PairVariant> out;
    if (num_logical_rows <= 0) return out;
    if (n.type != TransType::NMOS || p.type != TransType::PMOS) return out;

    auto minimal_split = [](int total, int maxp){ std::vector<int> v; while (total>0){ int t = std::min(maxp, total); v.push_back(t); total -= t; } return v; };

    const int nf_n = effective_nfin(n, tr);
    const int nf_p = effective_nfin(p, tr);

    std::vector<int> fn = minimal_split(nf_n, std::max(1, max_fins_per_row));
    std::vector<int> fp = minimal_split(nf_p, std::max(1, max_fins_per_row));

    const int N = (int)fn.size();
    const int P = (int)fp.size();
    const int K = std::min(N, P);
    const int M = std::max(N, P);

    struct Slot { bool hasN=false; int idxN=-1; bool hasP=false; int idxP=-1; };

    auto build_slots = [&](bool singles_prefix)->std::vector<Slot> {
        std::vector<Slot> slots; slots.resize(M);
        int iN=0, iP=0;
        if (singles_prefix) {
            int S = M - K; // singles first
            if (P > N) { // P-only first
                for (int s=0; s<S; ++s) { slots[s].hasP=true; slots[s].idxP=iP++; }
            } else if (N > P) {
                for (int s=0; s<S; ++s) { slots[s].hasN=true; slots[s].idxN=iN++; }
            }
            for (int s=S; s<M; ++s) { slots[s].hasN=true; slots[s].hasP=true; slots[s].idxN=iN++; slots[s].idxP=iP++; }
        } else {
            // pairs first, singles suffix
            for (int s=0; s<K; ++s) { slots[s].hasN=true; slots[s].hasP=true; slots[s].idxN=iN++; slots[s].idxP=iP++; }
            if (P > N) { for (int s=K; s<M; ++s) { slots[s].hasP=true; slots[s].idxP=iP++; } }
            else if (N > P) { for (int s=K; s<M; ++s) { slots[s].hasN=true; slots[s].idxN=iN++; } }
        }
        return slots;
    };

    struct Shape { int rows; int l0; int l1; int dx; }; // rows=1 or 2
    auto enumerate_shapes = [&](int M)->std::vector<Shape> {
        std::vector<Shape> sh;
        // 1-row shape
        sh.push_back({1, M, 0, 0});
        // 2-row shapes
        for (int l0=1; l0<=M-1; ++l0) {
            int l1 = M - l0;
            for (int dx = -(l1-1); dx <= (l0-1); ++dx) {
                sh.push_back({2, l0, l1, dx});
            }
        }
        return sh;
    };

    auto shapes = enumerate_shapes(M);

    // orientation patterns per row length
    auto row_patterns = [&](int count){ return row_orient_patterns(count); };

    for (int r=0; r<num_logical_rows; ++r) {
        int yN = phys_y_for(true, r);
        int yP = phys_y_for(false, r);
        for (bool singles_prefix : {false, true}) {
            auto slots = build_slots(singles_prefix);
            int singlesCount = M - K;
            for (const auto& sh : shapes) {
                // map slot s -> (row,col)
                std::vector<int> row_of(M, 0), col_of(M, 0);
                int cmin = 0, cmax = 0;
                if (sh.rows == 1) {
                    for (int s=0; s<M; ++s) { row_of[s]=0; col_of[s]=s; }
                    cmin = 0; cmax = M-1;
                } else {
                    int l0 = sh.l0, l1 = sh.l1, dx = sh.dx;
                    for (int s=0; s<M; ++s) {
                        if (s < l0) { row_of[s]=0; col_of[s]=s; }
                        else { row_of[s]=1; col_of[s]=dx + (s - l0); }
                    }
                    cmin = std::min(0, dx); cmax = std::max(l0-1, dx + l1 - 1);
                }
                // singles only at endpoints AND contiguous prefix/suffix in slot order
                if (singlesCount > 0) {
                    // slot-order contiguous already guaranteed by build_slots; additionally require column-endpoints
                    // collect columns of singles in slot order
                    std::vector<int> singlesCols; singlesCols.reserve(singlesCount);
                    if (singles_prefix) {
                        for (int s=0; s<singlesCount; ++s) singlesCols.push_back(col_of[s]);
                        // Must occupy left endpoint columns [cmin .. cmin+S-1]
                        bool ok=true; for (int j=0;j<singlesCount;++j) if (singlesCols[j] != (cmin + j)) { ok=false; break; }
                        if (!ok) continue;
                    } else {
                        for (int s=M-singlesCount; s<M; ++s) singlesCols.push_back(col_of[s]);
                        // Must occupy right endpoint columns [cmax-S+1 .. cmax]
                        bool ok=true; for (int j=0;j<singlesCount;++j) if (singlesCols[j] != (cmax - singlesCount + 1 + j)) { ok=false; break; }
                        if (!ok) continue;
                    }
                }
                // Build per-physical-row lists for orientation
                std::map<int, std::vector<int>> colsN, colsP; // y -> list of cols
                for (int s=0; s<M; ++s) {
                    if (slots[s].hasN) colsN[yN].push_back(col_of[s]);
                    if (slots[s].hasP) colsP[yP].push_back(col_of[s]);
                }
                // sort cols per row
                for (auto& kv : colsN) { auto& v = kv.second; std::sort(v.begin(), v.end()); }
                for (auto& kv : colsP) { auto& v = kv.second; std::sort(v.begin(), v.end()); }
                // enumerate starting orientation per row (two choices if row non-empty)
                std::vector<std::pair<int,std::vector<Orient>>> choices; // (y, pattern)
                std::vector<std::pair<int,std::vector<Orient>>> rows_all;
                for (auto& kv : colsN) { rows_all.push_back({kv.first, std::vector<Orient>()}); }
                for (auto& kv : colsP) { if (!colsN.count(kv.first)) rows_all.push_back({kv.first, std::vector<Orient>()}); }
                // rows_all has at most 2 entries (yN and yP)
                std::function<void(int)> rec_rows;
                rec_rows = [&](int idx){
                    if (idx == (int)rows_all.size()) {
                        // emit variant
                        PairVariant var; var.fingers.reserve(N+P);
                        auto get_orient = [&](int y, int col){
                            // find this row vector and its columns vector
                            const std::vector<Orient>* pat = nullptr; const std::vector<int>* cols=nullptr;
                            if (y == yN) { cols = &colsN[y]; }
                            else { cols = &colsP[y]; }
                            // locate col index
                            int pos = (int)(std::lower_bound(cols->begin(), cols->end(), col) - cols->begin());
                            for (const auto& yc : choices) if (yc.first==y) { pat = &yc.second; break; }
                            if (!pat || pat->empty()) return Orient::SD; return (*pat)[pos];
                        };
                        // x mapping: normalize by cmin to start from 0
                        auto to_x = [&](int col){ return (col - cmin) * 2; };
                        // Emit in slot order to keep stability
                        for (int s=0; s<M; ++s) {
                            int col = col_of[s];
                            if (slots[s].hasN) {
                                FingerPlacement fnp; fnp.base_id=n.id; fnp.finger_index=slots[s].idxN; fnp.nfin=fn[slots[s].idxN]; fnp.y=yN; fnp.x=to_x(col); fnp.orient=get_orient(yN, col); fnp.type=TransType::NMOS; var.fingers.push_back(fnp);
                            }
                            if (slots[s].hasP) {
                                FingerPlacement fpp; fpp.base_id=p.id; fpp.finger_index=slots[s].idxP; fpp.nfin=fp[slots[s].idxP]; fpp.y=yP; fpp.x=to_x(col); fpp.orient=get_orient(yP, col); fpp.type=TransType::PMOS; var.fingers.push_back(fpp);
                            }
                        }
                        std::sort(var.fingers.begin(), var.fingers.end(), [](const FingerPlacement& a, const FingerPlacement& b){ if (a.y!=b.y) return a.y<b.y; if (a.x!=b.x) return a.x<b.x; if (a.type!=b.type) return a.type<b.type; if (a.base_id!=b.base_id) return a.base_id<b.base_id; return a.finger_index<b.finger_index; });
                        out.push_back(std::move(var));
                        return;
                    }
                    int y = rows_all[idx].first;
                    const auto& cols = (y==yN ? colsN[y] : colsP[y]);
                    auto pats = row_patterns((int)cols.size());
                    for (const auto& ptn : pats) { choices.push_back({y, ptn}); rec_rows(idx+1); choices.pop_back(); }
                };
                rec_rows(0);
            }
        }
    }
    return out;
}

} // namespace stdcell

