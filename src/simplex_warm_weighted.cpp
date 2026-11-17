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
struct Edge { array<int,3> v; int w = 1; };
static int cell_index(int n, int a, int b) { return a*(n+1) - a*(a-1)/2 + b; }
static vector<string> read_rows(const string& path) {
    ifstream f(path); vector<string> rows;
    string text((istreambuf_iterator<char>(f)), istreambuf_iterator<char>()), s;
    size_t key = text.find("\"rows\"");
    if (key != string::npos) {
        size_t p = text.find('[', key), end = text.find(']', p);
        if (p == string::npos || end == string::npos) return {};
        while ((p = text.find('"', p + 1)) != string::npos && p < end) {
            size_t q = text.find('"', p + 1);
            if (q == string::npos || q > end) return {};
            s = text.substr(p + 1, q - p - 1);
            if (!s.empty() && s.find_first_not_of("0123") == string::npos) rows.push_back(s);
            p = q;
        }
    } else {
        istringstream input(text);
        while (input >> s) if (!s.empty() && s.find_first_not_of("0123") == string::npos) rows.push_back(s);
    }
    return rows;
}
struct Search {
    int n, r = 4, N; vector<Edge> edges; vector<vector<int>> inc;
    Search(int n_) : n(n_), N((n_+1)*(n_+2)/2), inc(N) {
        for (int k = 1; k <= n; ++k) for (int x = 0; x <= n-k; ++x) for (int y = 0; y <= n-k-x; ++y) {
            Edge e{{cell_index(n,x+k,y),cell_index(n,x,y+k),cell_index(n,x,y)}};
            int ei = edges.size(); edges.push_back(e);
            for (int v : e.v) inc[v].push_back(ei);
        }
    }
    vector<string> rows_of(const vector<int>& col) const {
        vector<string> rows;
        for (int a = 0; a <= n; ++a) {
            string s;
            for (int b = 0; b <= n-a; ++b) s.push_back(char('0'+col[cell_index(n,a,b)]));
            rows.push_back(s);
        }
        return rows;
    }
    int direct_bad(const vector<int>& col) const {
        int bad = 0;
        for (auto& e : edges) bad += col[e.v[0]] == col[e.v[1]] && col[e.v[1]] == col[e.v[2]];
        return bad;
    }
    void save(const string& path, const vector<int>& col, int bad, long long steps, double sec) const {
        ofstream f(path); auto rows = rows_of(col);
        f << "{\"t\":3,\"r\":4,\"n\":" << n << ",\"monochromatic_corners\":" << bad
          << ",\"corner_triples\":" << edges.size() << ",\"steps\":" << steps << ",\"seconds\":" << sec << ",\"rows\":[";
        for (int i = 0; i <= n; ++i) { if (i) f << ','; f << '"' << rows[i] << '"'; }
        f << "]}\n";
    }
    int run(const vector<int>& warm, double seconds, int seed, const string& out, int initial_perturb) {
        mt19937 rng(seed); auto start = chrono::steady_clock::now();
        auto elapsed = [&](){ return chrono::duration<double>(chrono::steady_clock::now()-start).count(); };
        vector<int> best = warm, col;
        if (initial_perturb > 0) for (int q = 0; q < initial_perturb; ++q) best[rng()%N] = rng()%r;
        int bestbad = direct_bad(best); long long steps = 0; int restart = 0;
        save(out,best,bestbad,0,0); cerr << "INITIAL n=" << n << " cells=" << N << " edges=" << edges.size() << " bad=" << bestbad << "\n";
        if (bestbad == 0) { cout << "FOUND initial\n"; return 0; }
        while (elapsed() < seconds) {
            col = best;
            int perturb = restart == 0 ? 0 : 1 + rng()%max(2,N/100);
            for (int q = 0; q < perturb; ++q) { int v = rng()%N; col[v] = rng()%r; }
            vector<array<unsigned char,4>> cnt(edges.size());
            vector<int> bad, where(edges.size(),-1), tabu(N*r,0);
            for (auto& e : edges) e.w = 1;
            auto setbad = [&](int ei, bool yes) {
                if (yes && where[ei] < 0) { where[ei] = bad.size(); bad.push_back(ei); }
                else if (!yes && where[ei] >= 0) {
                    int j = where[ei], z = bad.back(); bad[j] = z; where[z] = j; bad.pop_back(); where[ei] = -1;
                }
            };
            for (int ei = 0; ei < (int)edges.size(); ++ei) {
                cnt[ei].fill(0); for (int v : edges[ei].v) cnt[ei][col[v]]++;
                if (cnt[ei][col[edges[ei].v[0]]] == 3) setbad(ei,true);
            }
            int stale = 0;
            while (elapsed() < seconds && stale < 2000000) {
                ++steps; ++stale;
                if ((int)bad.size() < bestbad) {
                    int tracked = bad.size(), rebuilt = direct_bad(col);
                    if (tracked != rebuilt) { cerr << "INCREMENTAL_MISMATCH " << tracked << ' ' << rebuilt << "\n"; return 3; }
                    bestbad = tracked; best = col; stale = 0; save(out,best,bestbad,steps,elapsed());
                    cerr << "BEST bad=" << bestbad << " steps=" << steps << " restart=" << restart << " sec=" << elapsed() << "\n";
                    if (bestbad == 0) { cout << "FOUND steps=" << steps << "\n"; return 0; }
                }
                if (bad.empty()) break;
                int ei = bad[rng()%bad.size()];
                for (int q = 0; q < 8; ++q) { int z = bad[rng()%bad.size()]; if (edges[z].w > edges[ei].w) ei = z; }
                int chosenV = -1, chosenC = -1, bestGain = 1000000000;
                for (int v : edges[ei].v) {
                    int old = col[v];
                    for (int nw = 0; nw < r; ++nw) if (nw != old) {
                        int gain = 0;
                        for (int z : inc[v]) {
                            gain -= cnt[z][old] == 3 ? edges[z].w : 0;
                            gain += cnt[z][nw]+1 == 3 ? edges[z].w : 0;
                        }
                        bool tabu_move = tabu[v*r+nw] > steps;
                        if (tabu_move && (int)bad.size()+gain >= bestbad) continue;
                        if (gain < bestGain || (gain == bestGain && (rng()&1))) {
                            bestGain = gain; chosenV = v; chosenC = nw;
                        }
                    }
                }
                if (bestGain >= 0 && (rng()%1000) < 250) {
                    for (int z : bad) edges[z].w = min(edges[z].w+1,1000000);
                    continue;
                }
                if (chosenV < 0 || (bestGain > 0 && (rng()%1000) < 10)) {
                    chosenV = edges[ei].v[rng()%3]; chosenC = rng()%r;
                    if (chosenC == col[chosenV]) chosenC = (chosenC+1)%r;
                }
                int old = col[chosenV];
                for (int z : inc[chosenV]) {
                    setbad(z,false); cnt[z][old]--; cnt[z][chosenC]++;
                    if (cnt[z][chosenC] == 3) setbad(z,true);
                }
                col[chosenV] = chosenC; tabu[chosenV*r+old] = steps + 12 + rng()%29;
            }
            ++restart; cerr << "RESTART " << restart << " best=" << bestbad << " steps=" << steps << " sec=" << elapsed() << "\n";
        }
        save(out,best,bestbad,steps,elapsed()); cout << "DONE best=" << bestbad << " steps=" << steps << "\n"; return 1;
    }
};
int main(int argc, char** argv) {
    if (argc < 7) { cerr << "args n seconds seed output initial_rows direction [initial_perturb]\n"; return 2; }
    int n = stoi(argv[1]), seed = stoi(argv[3]), dir = stoi(argv[6]);
    double seconds = stod(argv[2]); string out = argv[4];
    auto oldrows = read_rows(argv[5]); int oldn = oldrows.size()-1;
    if (oldn != n-1 && oldn != n) { cerr << "need T_(n-1) rows, or T_n for a perturbation control\n"; return 2; }
    vector<int> oldcol;
    for (int a = 0; a <= oldn; ++a) for (int b = 0; b <= oldn-a; ++b) oldcol.push_back(oldrows[a][b]-'0');
    vector<int> warm((n+1)*(n+2)/2);
    if (oldn == n) warm = oldcol;
    else for (int a = 0; a <= n; ++a) for (int b = 0; b <= n-a; ++b) {
        array<int,3> v{{a,b,n-a-b}}, u{}; bool copied = false;
        for (int q = 0; q < 3; ++q) {
            int d = (dir+q)%3; if (v[d] == 0) continue;
            u = v; u[d]--;
            warm[cell_index(n,a,b)] = oldcol[cell_index(oldn,u[0],u[1])]; copied = true; break;
        }
        if (!copied) return 2;
    }
    Search s(n); return s.run(warm,seconds,seed,out,argc > 7 ? stoi(argv[7]) : 0);
}
