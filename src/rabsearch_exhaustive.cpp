#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;

static vector<int> prime_factors(int n) {
    vector<int> f;
    for (int q = 2; q * (long long)q <= n; ++q) if (n % q == 0) {
        f.push_back(q);
        while (n % q == 0) n /= q;
    }
    if (n > 1) f.push_back(n);
    return f;
}
static int modpow(long long a, long long e, int p) {
    long long z = 1;
    while (e) {
        if (e & 1) z = z * a % p;
        a = a * a % p;
        e >>= 1;
    }
    return (int)z;
}
static bool is_prime(int n) {
    if (n < 2) return false;
    for (int q = 2; q * (long long)q <= n; ++q) if (n % q == 0) return false;
    return true;
}
static int primitive_root(int p) {
    auto f = prime_factors(p - 1);
    for (int g = 2; g < p; ++g) {
        bool ok = true;
        for (int q : f) if (modpow(g, (p - 1) / q, p) == 1) { ok = false; break; }
        if (ok) return g;
    }
    return -1;
}

struct Scanner {
    int p, r, g, W;
    vector<int> chi;
    vector<uint64_t> E; // E[c,s] = {y : chi(y+s)=c}, p bits per (c,s).
    long long prefixes = 0, unions = 0;
    Scanner(int p_, int r_) : p(p_), r(r_), g(primitive_root(p_)), W((p_ + 63) / 64), chi(p_) {
        chi[0] = r;
        int x = 1;
        for (int i = 0; i < p - 1; ++i) { chi[x] = i % r; x = (int)((long long)x * g % p); }
        E.assign((size_t)r * p * W, 0);
        for (int s = 0; s < p; ++s) for (int x = 1; x < p; ++x) {
            int y = x - s; if (y < 0) y += p;
            E[((size_t)chi[x] * p + s) * W + (y >> 6)] |= 1ULL << (y & 63);
        }
    }
    inline const uint64_t* bits(int c, int s) const { return &E[((size_t)c * p + s) * W]; }
    inline bool empty_and_complement(uint64_t* a, const uint64_t* b) const {
        uint64_t any = 0;
        for (int w = 0; w < W; ++w) { a[w] &= ~b[w]; any |= a[w]; }
        ++const_cast<long long&>(unions);
        return any == 0;
    }
    // Return the least legal final element greater than b, or -1. Prefix is {0,1,a,b}.
    int final_element(int a, int b, const vector<uint64_t>& I3) {
        vector<uint64_t> allowed(W), I(W);
        for (int w = 0; w < W; ++w) {
            int lo = w * 64;
            uint64_t m = ~0ULL;
            if (lo + 63 >= p) m &= (p % 64 == 0 ? ~0ULL : ((1ULL << (p % 64)) - 1));
            if (lo <= b) {
                int k = b - lo;
                m = k >= 63 ? 0 : ~((1ULL << (k + 1)) - 1);
                if (lo + 63 >= p && p % 64) m &= (1ULL << (p % 64)) - 1;
            }
            allowed[w] = m;
        }
        const int S[4] = {0, 1, a, b};
        for (int c = 0; c < r; ++c) {
            // Failure (C): a prefix value is zero and all remaining values have colour c.
            for (int j = 0; j < 4; ++j) {
                int s0 = S[j]; bool mono = true;
                for (int h = 0; h < 4; ++h) if (h != j) {
                    int d = S[h] - s0; if (d < 0) d += p;
                    if (d == 0 || chi[d] != c) { mono = false; break; }
                }
                if (mono) {
                    int sh = (-s0) % p; if (sh < 0) sh += p;
                    if (empty_and_complement(allowed.data(), bits(c, sh))) return -1;
                }
            }
            // I = {y : all four prefix values have colour c}.
            const uint64_t* eb = bits(c, b);
            for (int w = 0; w < W; ++w) I[w] = I3[c * W + w] & eb[w];
            for (int w = 0; w < W; ++w) {
                uint64_t z = I[w];
                while (z) {
                    int y = w * 64 + __builtin_ctzll(z); z &= z - 1;
                    if (y >= p) continue;
                    // Failure (A): all five values have colour c.
                    if (empty_and_complement(allowed.data(), bits(c, y))) return -1;
                    // Failure (B): final value is zero and prefix has colour c.
                    int x = (-y) % p; if (x < 0) x += p;
                    allowed[x >> 6] &= ~(1ULL << (x & 63));
                }
            }
        }
        for (int w = 0; w < W; ++w) if (allowed[w]) return w * 64 + __builtin_ctzll(allowed[w]);
        return -1;
    }
    vector<int> scan() {
        auto start = chrono::steady_clock::now();
        vector<uint64_t> I3((size_t)r * W);
        for (int a = 2; a < p - 2; ++a) {
            for (int c = 0; c < r; ++c) {
                const uint64_t *e0 = bits(c, 0), *e1 = bits(c, 1), *ea = bits(c, a);
                for (int w = 0; w < W; ++w) I3[c * W + w] = e0[w] & e1[w] & ea[w];
            }
            for (int b = a + 1; b < p - 1; ++b) {
                ++prefixes;
                int z = final_element(a, b, I3);
                if (z > b) {
                    vector<int> S = {0, 1, a, b, z};
                    double sec = chrono::duration<double>(chrono::steady_clock::now() - start).count();
                    cerr << "HIT p=" << p << " prefixes=" << prefixes << " unions=" << unions << " seconds=" << sec << "\n";
                    return S;
                }
            }
            if (a % 250 == 0) {
                double sec = chrono::duration<double>(chrono::steady_clock::now() - start).count();
                cerr << "PROGRESS p=" << p << " a=" << a << " prefixes=" << prefixes << " unions=" << unions << " seconds=" << sec << "\n";
            }
        }
        double sec = chrono::duration<double>(chrono::steady_clock::now() - start).count();
        cerr << "EXHAUSTED p=" << p << " prefixes=" << prefixes << " unions=" << unions << " seconds=" << sec << "\n";
        return {};
    }
    bool direct(const vector<int>& S) const {
        for (int y = 0; y < p; ++y) {
            int only = -1, n = 0; bool ok = false;
            for (int s : S) {
                int x = y + s; if (x >= p) x -= p;
                if (x == 0) continue;
                ++n;
                if (only < 0) only = chi[x];
                else if (only != chi[x]) { ok = true; break; }
            }
            if (!ok && n > 0) return false;
        }
        return true;
    }
};

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "args: p r output.json  (currently exhaustive for t=5; requires p prime)\n";
        return 2;
    }
    int p = stoi(argv[1]), r = stoi(argv[2]); string out = argv[3];
    if (!is_prime(p) || (p - 1) % r) { cerr << "invalid prime or r does not divide p-1\n"; return 2; }
    Scanner s(p, r);
    auto S = s.scan();
    if (S.empty()) return 1;
    if (!s.direct(S)) { cerr << "DIRECT_REJECT\n"; return 3; }
    ofstream f(out);
    f << "{\"p\":" << p << ",\"r\":" << r << ",\"t\":5,\"g\":" << s.g << ",\"bad\":0,\"weights\":[";
    for (int i = 0; i < 5; ++i) { if (i) f << ','; f << S[i]; }
    f << "],\"prefixes\":" << s.prefixes << ",\"unions\":" << s.unions << ",\"engine\":\"exhaustive-bitset-prefix\",\"direct_recheck\":true}\n";
    cout << "FOUND p=" << p << " r=" << r << " S=";
    for (int x : S) cout << x << ',';
    cout << "\n";
    return 0;
}
