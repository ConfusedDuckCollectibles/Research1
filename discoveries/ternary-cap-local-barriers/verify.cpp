// Self-contained replay for local exchange barriers for the classical 236-cap.
// C++17 standard library only. Optional --certificate prints line packings and colors.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr int dim = 7, n = 2187, m = 236;
using Vec = std::array<int, dim>;
using Edge = std::pair<int, int>;
std::array<Vec, n> v;
std::array<std::array<std::uint16_t, n>, n> third;
std::array<int, n> baseline_index;
std::array<std::vector<Edge>, n> secants;

void require(bool ok, const char* reason) {
    if (!ok) throw std::runtime_error(reason);
}

int pack(const Vec& x) {
    int out = 0;
    for (int digit : x) out = 3 * out + digit;
    return out;
}

int block(int i, int j, int k) {
    return (1 << (i - 1)) | (1 << (j - 1)) | (1 << (k - 1));
}

std::vector<int> baseline() {
    // Generate the classical design via two orbits of (2 3 5 6 4).
    std::array<bool, 64> D{};
    const std::array<int, 7> perm = {0, 1, 3, 5, 2, 6, 4};
    for (int seed : {block(1, 2, 3), block(2, 3, 6)}) {
        int support = seed;
        for (int r = 0; r < 5; ++r) {
            require(!D[support], "support orbits overlap");
            D[support] = true;
            int next = 0;
            for (int i = 1; i <= 6; ++i)
                if ((support >> (i - 1)) & 1) next |= 1 << (perm[i] - 1);
            support = next;
        }
        require(support == seed, "orbit did not close");
    }
    std::vector<int> out;
    for (int x = 0; x < n; ++x) {
        int support = 0, weight = 0, negatives = 0;
        for (int i = 1; i <= 6; ++i) {
            if (v[x][i]) { support |= 1 << (i - 1); ++weight; }
            negatives += v[x][i] == 2;
        }
        const int t = v[x][0];
        const bool in = t == 0 ? weight == 1
            : (weight == 6 && negatives % 2 == 0)
                || (weight == 3 && D[t == 1 ? support : 63 ^ support]);
        if (in) out.push_back(x);
    }
    return out;
}

int exact_pair_cost(int a, int b) {
    const int forced = baseline_index[third[a][b]];
    // Disjoint-set unions gather components of the two residual matchings.
    std::array<int, m> parent, extent, touched;
    parent.fill(-1);
    int count = 0;
    auto make = [&](int x) {
        if (parent[x] < 0) { parent[x] = x; extent[x] = 1; touched[count++] = x; }
    };
    auto root = [&](int x) {
        while (x != parent[x]) x = parent[x];
        return x;
    };
    for (int z : {a, b}) for (auto [u, w] : secants[z]) {
        if (u == forced || w == forced) continue;
        make(u);
        make(w);
        u = root(u);
        w = root(w);
        if (u == w) continue;
        if (extent[u] < extent[w]) std::swap(u, w);
        parent[w] = u;
        extent[u] += extent[w];
    }
    int cost = forced >= 0;
    for (int i = 0; i < count; ++i) {
        const int x = touched[i];
        if (parent[x] == x) cost += extent[x] / 2;
    }
    return cost;
}

void audit_exchange(bool certificate) {
    std::vector<int> outside;
    int min_single = n, attaining = 0;
    for (int i = 0; i < n; ++i) if (baseline_index[i] < 0) {
        outside.push_back(i);
        const int d = static_cast<int>(secants[i].size());
        if (d < min_single) { min_single = d; attaining = 0; }
        attaining += d == min_single;
    }
    require(outside.size() == 1951 && min_single == 10 && attaining == 24,
            "single-point counts disagree");
    std::size_t pairs = 0, edges = 0;
    int minimum = n;
    std::array<std::vector<int>, n> adj;
    std::vector<int> eligible;
    for (int x : outside) if (secants[x].size() <= 16) eligible.push_back(x);
    for (int x : outside) for (int y : outside) if (x < y) {
        ++pairs;
        const int cost = exact_pair_cost(x, y);
        minimum = std::min(minimum, cost);
        if (cost <= 16 && secants[x].size() <= 16 && secants[y].size() <= 16) {
            ++edges;
            adj[x].push_back(y);
            adj[y].push_back(x);
        }
    }
    require(pairs == 1902225 && minimum == 16, "pair counts disagree");
    require(eligible.size() == 1616 && edges == 888, "compatibility graph disagrees");

    auto color_in_order = [&](const std::vector<int>& order) {
        std::array<int, n> color;
        color.fill(-1);
        int used = 0;
        for (int x : order) {
            std::vector<bool> blocked(used, false);
            for (int y : adj[x]) if (color[y] >= 0) blocked[color[y]] = true;
            const int c = static_cast<int>(std::find(blocked.begin(), blocked.end(), false)
                                            - blocked.begin());
            color[x] = c;
            used += c == used;
        }
        return std::make_pair(used, color);
    };
    std::vector<int> order = eligible;
    auto coloring = color_in_order(order);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return adj[a].size() != adj[b].size() ? adj[a].size() > adj[b].size() : a < b;
    });
    auto alt = color_in_order(order);
    if (alt.first < coloring.first) coloring = alt;
    require(coloring.first == 7, "coloring bound disagrees");
    for (int x : eligible) {
        require(coloring.second[x] >= 0 && coloring.second[x] < 7,
                "missing or invalid color");
        for (int y : adj[x]) require(coloring.second[x] != coloring.second[y],
                                     "coloring not proper");
    }
    std::cout << "Single outsider: minimum 10 secants, attained by 24 of 1,951.\n"
              << "Outsider pairs: 1,902,225 checked; minimum exact cost 16.\n"
              << "Budget 16: 1,616 vertices, 888 edges, verified 7-coloring.\n";
    if (certificate)
        for (int x : eligible) std::cout << "color " << x << ' ' << coloring.second[x] << '\n';
}

void audit_hyperplanes(const std::vector<int>& S, bool certificate) {
    std::map<int, int> hist;
    for (int normal = 1; normal < n; ++normal) {
        const auto nonzero = std::find_if(v[normal].begin(), v[normal].end(),
                                         [](int x) { return x != 0; });
        if (*nonzero != 1) continue;
        std::array<int, n> val{};
        for (int x = 0; x < n; ++x)
            val[x] = std::inner_product(v[normal].begin(), v[normal].end(),
                                        v[x].begin(), 0) % 3;
        std::array<std::vector<int>, 3> old;
        for (int x : S) old[val[x]].push_back(x);
        for (int h = 0; h < 3; ++h) {
            std::array<bool, n> blocked{}, used{};
            for (int x : old[(h + 1) % 3]) for (int y : old[(h + 2) % 3])
                blocked[third[x][y]] = true;
            for (int x : old[h]) require(!blocked[x], "blocked baseline point");

            int pool_size = 0;
            std::vector<std::array<int, 3>> packing;
            // Match each allowed outsider to two unused old points.
            // Candidate iteration is in reverse coordinate significance order.
            for (int serial = 0; serial < n; ++serial) {
                int x = 0, copy = serial;
                for (int i = 0; i < dim; ++i) { x = 3*x + copy%3; copy /= 3; }
                if (val[x] != h || blocked[x]) continue;
                ++pool_size;
                if (baseline_index[x] >= 0) continue;
                bool found = false;
                for (int order = 0; order < n; ++order) {
                    int y = 0, rep = order;
                    for (int i = 0; i < dim; ++i) { y = 3*y + rep%3; rep /= 3; }
                    if (baseline_index[y] < 0 || val[y] != h || used[y]) continue;
                    const int z = third[x][y];
                    if (val[z] != h || baseline_index[z] < 0 || used[z]) continue;
                    require(x != y && y != z && x != z && !blocked[y] && !blocked[z],
                            "invalid proposed packing triple");
                    used[y] = used[z] = true;
                    packing.push_back({x, y, z});
                    found = true;
                    break;
                }
                require(found, "could not generate disjoint packing");
            }
            const int size = static_cast<int>(old[h].size());
            const int slack = pool_size - size;
            require(static_cast<int>(packing.size()) == slack, "wrong packing size");
            std::array<bool, n> verify_used{};
            for (const auto& triple : packing) {
                for (int x : triple) {
                    require(!verify_used[x] && !blocked[x] && val[x] == h,
                            "packing not disjoint or not in pool");
                    verify_used[x] = true;
                }
                require(third[triple[0]][triple[1]] == triple[2], "noncollinear triple");
            }
            ++hist[slack];
            if (certificate) {
                std::cout << "plane " << normal << ' ' << h << ' ' << size << ' ' << pool_size;
                for (const auto& triple : packing) for (int x : triple) std::cout << ' ' << x;
                std::cout << '\n';
            }
        }
    }
    require(hist == std::map<int, int>{{0,1669},{1,1220},{2,388},{13,2}},
            "hyperplane counts disagree");
    std::cout << "Hyperplanes: all 3,279 verified; slack 0/1/2/13 on "
                 "1,669/1,220/388/2 hyperplanes.\n";
}
} // namespace

int main(int argc, char** argv) try {
    if (argc > 2 || (argc == 2 && std::string(argv[1]) != "--certificate")) return 2;
    const bool certificate = argc == 2;
    for (int i = 0; i < n; ++i) {
        int copy = i;
        for (int j = dim - 1; j >= 0; --j) { v[i][j] = copy % 3; copy /= 3; }
    }
    for (int x = 0; x < n; ++x) for (int y = x; y < n; ++y) {
        Vec z{};
        for (int j = 0; j < dim; ++j) z[j] = (6 - v[x][j] - v[y][j]) % 3;
        third[x][y] = third[y][x] = static_cast<std::uint16_t>(pack(z));
    }
    const auto S = baseline();
    require(S.size() == m, "baseline has wrong size");
    baseline_index.fill(-1);
    for (int i = 0; i < m; ++i) baseline_index[S[i]] = i;
    int source_pairs = 0;
    for (int i = 0; i < m; ++i) for (int j = i + 1; j < m; ++j) {
        const int z = third[S[i]][S[j]];
        require(baseline_index[z] < 0, "baseline has a forbidden triple");
        secants[z].push_back({i, j});
        ++source_pairs;
    }
    require(source_pairs == 27730, "baseline pair count disagrees");
    for (const auto& rows : secants) {
        std::array<bool, m> used{};
        for (auto [i, j] : rows) {
            require(i != j && !used[i] && !used[j], "secants not a matching");
            used[i] = used[j] = true;
        }
    }
    std::cout << "Baseline: 236 points, all 27,730 pairs checked.\n";
    audit_exchange(certificate);
    audit_hyperplanes(S, certificate);
    std::cout << "PASS: both scoped local barriers verified.\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
}
