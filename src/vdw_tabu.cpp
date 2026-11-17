// Whole-colouring tabu search for 3-term van der Waerden certificates.
// Build: g++ -O3 -std=c++17 vdw_tabu.cpp -o vdw_tabu
// Run:   vdw_tabu N r seconds seed output.json [warm.json] [prepend]
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct Edge { array<int, 3> v; };

static string slurp(const string& path) {
    ifstream f(path);
    if (!f) throw runtime_error("cannot open " + path);
    return string((istreambuf_iterator<char>(f)), istreambuf_iterator<char>());
}

static string json_string(const string& text, const string& key) {
    size_t p = text.find("\"" + key + "\"");
    if (p == string::npos) return {};
    p = text.find(':', p);
    p = text.find('"', p);
    if (p == string::npos) throw runtime_error("missing JSON string " + key);
    string out;
    for (++p; p < text.size() && text[p] != '"'; ++p) out.push_back(text[p]);
    return out;
}

static int digit_value(char c) {
    if (c >= 'a' && c <= 'z') return c - 'a';
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'Z') return c - 'A' + 10;
    throw runtime_error("invalid colour character");
}

static vector<int> read_colouring(const string& path) {
    string text = slurp(path), s = json_string(text, "coloring");
    if (s.empty()) s = json_string(text, "palette");
    if (s.empty()) throw runtime_error("warm JSON needs coloring or palette");
    vector<int> out;
    for (char c : s) out.push_back(digit_value(c));
    return out;
}

struct Search {
    int N, r;
    vector<Edge> edges;
    vector<vector<int>> incident;
    mt19937 rng;

    Search(int N_, int r_, int seed) : N(N_), r(r_), incident(N_), rng(seed) {
        if (r < 2 || r > 26 || N < 3) throw runtime_error("require 2<=r<=26 and N>=3");
        for (int d = 1; 2 * d < N; ++d) {
            for (int a = 0; a + 2 * d < N; ++a) {
                Edge e{{a, a + d, a + 2 * d}};
                int ei = edges.size();
                edges.push_back(e);
                for (int v : e.v) incident[v].push_back(ei);
            }
        }
    }

    int direct_bad(const vector<int>& col) const {
        int bad = 0;
        for (const auto& e : edges)
            bad += col[e.v[0]] == col[e.v[1]] && col[e.v[1]] == col[e.v[2]];
        return bad;
    }

    void save(const string& path, const vector<int>& col, int bad,
              long long steps, double seconds) const {
        string text;
        for (int c : col) text.push_back(char('a' + c));
        ofstream f(path);
        f << "{\"kind\":\"van-der-Waerden-coloring\",\"r\":" << r
          << ",\"k\":3,\"N\":" << N << ",\"positions\":\"1,...,N\","
          << "\"alphabet\":\"a=0, b=1, ...\",\"coloring\":\"" << text << "\","
          << "\"monochromatic_progressions\":" << bad << ",\"steps\":" << steps
          << ",\"seconds\":" << seconds << "}\n";
    }

    int run(double seconds, const string& output, vector<int> col) {
        auto start = chrono::steady_clock::now();
        auto elapsed = [&]() {
            return chrono::duration<double>(chrono::steady_clock::now() - start).count();
        };
        if ((int)col.size() != N) throw runtime_error("wrong initial colouring length");
        vector<int> best = col;
        int bestbad = direct_bad(best);
        long long steps = 0, stale = 0;
        int restarts = 0;
        save(output, best, bestbad, 0, 0);
        cerr << "INITIAL N=" << N << " r=" << r << " APs=" << edges.size()
             << " bad=" << bestbad << "\n";
        if (bestbad == 0) { cout << "FOUND initial\n"; return 0; }

        while (elapsed() < seconds) {
            if (restarts) {
                col = best;
                int perturb = max(2, N / 100);
                for (int q = 0; q < perturb; ++q) col[rng() % N] = rng() % r;
            }
            vector<unsigned char> cnt((size_t)edges.size() * r);
            vector<int> bad, where(edges.size(), -1), tabu_until(N, 0);
            auto is_mono = [&](int e) {
                for (int c = 0; c < r; ++c) if (cnt[(size_t)e * r + c] == 3) return true;
                return false;
            };
            auto set_bad = [&](int e, bool yes) {
                if (yes && where[e] < 0) {
                    where[e] = bad.size(); bad.push_back(e);
                } else if (!yes && where[e] >= 0) {
                    int j = where[e], z = bad.back();
                    bad[j] = z; where[z] = j; bad.pop_back(); where[e] = -1;
                }
            };
            for (int e = 0; e < (int)edges.size(); ++e) {
                for (int v : edges[e].v) ++cnt[(size_t)e * r + col[v]];
                if (is_mono(e)) set_bad(e, true);
            }

            while (elapsed() < seconds && stale < 200000 && !bad.empty()) {
                ++steps; ++stale;
                int e = bad[rng() % bad.size()];
                int chosen_v = -1, chosen_c = -1, best_gain = 1000000000;
                // With probability 2%, make a random non-tabu move in a bad AP.
                if (rng() % 100 < 2) {
                    vector<int> candidates;
                    for (int v : edges[e].v) if (tabu_until[v] <= steps) candidates.push_back(v);
                    if (!candidates.empty()) {
                        chosen_v = candidates[rng() % candidates.size()];
                        chosen_c = rng() % (r - 1);
                        if (chosen_c >= col[chosen_v]) ++chosen_c;
                    }
                }
                if (chosen_v < 0) {
                    for (int v : edges[e].v) {
                        int old = col[v];
                        for (int nw = 0; nw < r; ++nw) if (nw != old) {
                            int gain = 0;
                            for (int z : incident[v]) {
                                gain -= is_mono(z) ? 1 : 0;
                                gain += cnt[(size_t)z * r + nw] == 2 ? 1 : 0;
                            }
                            bool tabu = tabu_until[v] > steps;
                            if (tabu && (int)bad.size() + gain >= bestbad) continue; // aspiration
                            if (gain < best_gain ||
                                (gain == best_gain && chosen_v >= 0 && (rng() & 1))) {
                                best_gain = gain; chosen_v = v; chosen_c = nw;
                            }
                        }
                    }
                }
                if (chosen_v < 0) continue;
                int old = col[chosen_v];
                for (int z : incident[chosen_v]) {
                    set_bad(z, false);
                    --cnt[(size_t)z * r + old];
                    ++cnt[(size_t)z * r + chosen_c];
                    if (is_mono(z)) set_bad(z, true);
                }
                col[chosen_v] = chosen_c;
                tabu_until[chosen_v] = steps + 2 + rng() % 7; // tenure 2,...,8
                if ((int)bad.size() < bestbad) {
                    int tracked = bad.size(), rebuilt = direct_bad(col);
                    if (tracked != rebuilt) {
                        cerr << "INCREMENTAL_MISMATCH " << tracked << " " << rebuilt << "\n";
                        return 3;
                    }
                    bestbad = tracked; best = col; stale = 0;
                    save(output, best, bestbad, steps, elapsed());
                    cerr << "BEST bad=" << bestbad << " steps=" << steps
                         << " restarts=" << restarts << " sec=" << elapsed() << "\n";
                    if (bestbad == 0) { cout << "FOUND steps=" << steps << "\n"; return 0; }
                }
            }
            ++restarts;
        }
        save(output, best, bestbad, steps, elapsed());
        cout << "DONE best=" << bestbad << " steps=" << steps << " restarts=" << restarts << "\n";
        return 1;
    }
};

int main(int argc, char** argv) {
    try {
        if (argc < 6) {
            cerr << "args: N r seconds seed output.json [warm.json] [prepend]\n";
            return 2;
        }
        int N = stoi(argv[1]), r = stoi(argv[2]), seed = stoi(argv[4]);
        double seconds = stod(argv[3]);
        string output = argv[5];
        mt19937 initial_rng(seed);
        vector<int> col(N);
        for (int& c : col) c = initial_rng() % r;
        if (argc > 6) {
            vector<int> warm = read_colouring(argv[6]);
            bool prepend = argc > 7 && stoi(argv[7]) != 0;
            if ((int)warm.size() == N) col = warm;
            else if ((int)warm.size() == N - 1) {
                col = vector<int>(N);
                for (int i = 0; i < N - 1; ++i) col[i + prepend] = warm[i];
                col[prepend ? 0 : N - 1] = initial_rng() % r;
            } else throw runtime_error("warm length must be N or N-1");
        }
        for (int c : col) if (c < 0 || c >= r) throw runtime_error("warm colour out of range");
        Search search(N, r, seed);
        return search.run(seconds, output, col);
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}
