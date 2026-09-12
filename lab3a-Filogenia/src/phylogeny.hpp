#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace phylo {
using Matrix = std::vector<std::vector<double>>;
constexpr double EPS = 1e-12;
inline std::string number(double x) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(6) << x;
    return o.str();
}
inline void validate(const Matrix &d) {
    if (d.size() < 2)
        throw std::invalid_argument("Se requieren al menos dos taxones");
    for (const auto &row : d)
        if (row.size() != d.size())
            throw std::invalid_argument("Matriz no cuadrada");
    for (size_t i = 0; i < d.size(); ++i)
        for (size_t j = 0; j < d.size(); ++j)
            if (!std::isfinite(d[i][j]) || d[i][j] < 0 || std::abs(d[i][j] - d[j][i]) > EPS ||
                (i == j && std::abs(d[i][j]) > EPS))
                throw std::invalid_argument("Distancias invalidas: deben ser finitas, simetricas, "
                                            "no negativas y diagonal cero");
}
struct Alignment {
    std::string a, b;
    int score = 0, compared = 0, differences = 0, gaps = 0;
};
// Needleman-Wunsch global. Empates: diagonal, arriba, izquierda.
inline Alignment align(const std::string &a, const std::string &b, int gap = -2) {
    if (a.empty() || b.empty())
        throw std::invalid_argument("Cadena vacia");
    size_t n = a.size(), m = b.size();
    std::vector<std::vector<int>> f(n + 1, std::vector<int>(m + 1));
    for (size_t i = 1; i <= n; ++i)
        f[i][0] = int(i) * gap;
    for (size_t j = 1; j <= m; ++j)
        f[0][j] = int(j) * gap;
    for (size_t i = 1; i <= n; ++i)
        for (size_t j = 1; j <= m; ++j)
            f[i][j] = std::max({f[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 2 : -1),
                                f[i - 1][j] + gap, f[i][j - 1] + gap});
    Alignment r;
    r.score = f[n][m];
    size_t i = n, j = m;
    while (i || j) {
        if (i && j && f[i][j] == f[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 2 : -1)) {
            r.a += a[--i];
            r.b += b[--j];
        } else if (i && f[i][j] == f[i - 1][j] + gap) {
            r.a += a[--i];
            r.b += '-';
        } else {
            r.a += '-';
            r.b += b[--j];
        }
    }
    std::reverse(r.a.begin(), r.a.end());
    std::reverse(r.b.begin(), r.b.end());
    for (size_t k = 0; k < r.a.size(); ++k) {
        if (r.a[k] == '-' || r.b[k] == '-')
            ++r.gaps;
        else {
            ++r.compared;
            if (r.a[k] != r.b[k])
                ++r.differences;
        }
    }
    return r;
}
inline double jukesCantor(double p) {
    if (!std::isfinite(p) || p < 0 || p >= 0.75)
        throw std::domain_error(
            "Jukes-Cantor requiere 0 <= p < 0.75; no se reemplaza por una constante");
    return -0.75 * std::log1p(-4.0 * p / 3.0);
}
inline std::string normalize(std::string s) {
    // Normalizacion UTF-8 explicita para apellidos en espanol.
    const std::vector<std::pair<std::string, std::string>> replacements = {
        {"á", "a"}, {"é", "e"}, {"í", "i"}, {"ó", "o"}, {"ú", "u"}, {"ü", "u"}, {"ñ", "n"},
        {"Á", "a"}, {"É", "e"}, {"Í", "i"}, {"Ó", "o"}, {"Ú", "u"}, {"Ü", "u"}, {"Ñ", "n"}};
    for (const auto &p : replacements) {
        size_t k = 0;
        while ((k = s.find(p.first, k)) != std::string::npos) {
            s.replace(k, p.first.size(), p.second);
            k += p.second.size();
        }
    }
    std::string out;
    for (unsigned char c : s) {
        if (c >= 'A' && c <= 'Z')
            c = static_cast<unsigned char>(c - 'A' + 'a');
        if (c >= 'a' && c <= 'z')
            out += char(c);
        else if (c != ' ' && c != '-' && c != '\'' && c != '\r' && c != '\t')
            throw std::invalid_argument("Caracter no soportado en apellido");
    }
    if (out.empty())
        throw std::invalid_argument("Apellido vacio");
    return out;
}
inline double editDistance(const std::string &a, const std::string &b) {
    if (a.empty() || b.empty())
        throw std::invalid_argument("Apellido vacio");
    std::vector<size_t> p(b.size() + 1), q(b.size() + 1);
    std::iota(p.begin(), p.end(), 0);
    for (size_t i = 1; i <= a.size(); ++i) {
        q[0] = i;
        for (size_t j = 1; j <= b.size(); ++j)
            q[j] = std::min({p[j] + 1, q[j - 1] + 1, p[j - 1] + size_t(a[i - 1] != b[j - 1])});
        p.swap(q);
    }
    return double(p.back()) / double(std::max(a.size(), b.size()));
}
struct Person {
    std::string id, paternal, maternal;
};
inline double surnameDistance(const Person &a, const Person &b) {
    const auto ap = normalize(a.paternal), am = normalize(a.maternal), bp = normalize(b.paternal),
               bm = normalize(b.maternal);
    return 0.5 * std::min(editDistance(ap, bp) + editDistance(am, bm),
                          editDistance(ap, bm) + editDistance(am, bp));
}
struct Edge {
    int child;
    double length;
};
struct Node {
    std::string label;
    std::vector<Edge> children;
    double height = 0;
    size_t size = 1;
};
struct Tree {
    std::vector<Node> nodes;
    int root = -1;
    bool rooted = true;
    std::vector<std::string> trace;
};
inline Tree build(Matrix d, const std::vector<std::string> &labels, bool nj) {
    validate(d);
    if (labels.size() != d.size())
        throw std::invalid_argument("Numero de etiquetas incorrecto");
    Tree t;
    t.rooted = !nj;
    std::vector<int> active;
    for (const auto &s : labels) {
        active.push_back(int(t.nodes.size()));
        t.nodes.push_back({s, {}, 0, 1});
    }
    while (active.size() > 1) {
        size_t n = active.size(), a = 0, b = 1;
        std::vector<double> sums(n, 0);
        if (nj)
            for (size_t i = 0; i < n; ++i)
                sums[i] = std::accumulate(d[i].begin(), d[i].end(), 0.0);
        auto criterion = [&](size_t i, size_t j) {
            return nj && n > 2 ? double(n - 2) * d[i][j] - sums[i] - sums[j] : d[i][j];
        };
        double best = criterion(a, b);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = i + 1; j < n; ++j)
                if (criterion(i, j) < best - EPS) {
                    a = i;
                    b = j;
                    best = criterion(i, j);
                }
        const Node na = t.nodes[active[a]], nb = t.nodes[active[b]];
        double h = d[a][b] / 2.0, la = h - na.height, lb = h - nb.height;
        if (nj) {
            la = h + (n > 2 ? (sums[a] - sums[b]) / (2.0 * double(n - 2)) : 0);
            lb = d[a][b] - la;
        }
        int id = int(t.nodes.size());
        std::string name = "(" + na.label + "," + nb.label + ")";
        t.nodes.push_back({name, {{active[a], la}, {active[b], lb}}, h, na.size + nb.size});
        std::string trace = "n=" + std::to_string(n) + "; unir " + na.label + " + " + nb.label +
                            "; criterio=" + number(best) + "; d=" + number(d[a][b]) +
                            "; ramas=" + number(la) + "," + number(lb);
        if (nj && n > 2)
            trace += "; sumas=" + number(sums[a]) + "," + number(sums[b]);
        if (!nj)
            trace += "; tamanos=" + std::to_string(na.size) + "," + std::to_string(nb.size) +
                     "; altura=" + number(h);
        t.trace.push_back(trace);
        std::vector<size_t> keep;
        for (size_t k = 0; k < n; ++k)
            if (k != a && k != b)
                keep.push_back(k);
        Matrix next(n - 1, std::vector<double>(n - 1));
        std::vector<int> nextActive;
        for (size_t i = 0; i < keep.size(); ++i) {
            size_t k = keep[i];
            nextActive.push_back(active[k]);
            for (size_t j = 0; j < keep.size(); ++j)
                next[i][j] = d[k][keep[j]];
            double v = nj ? (d[a][k] + d[b][k] - d[a][b]) / 2.0
                          : (double(na.size) * d[a][k] + double(nb.size) * d[b][k]) /
                                double(na.size + nb.size);
            next[i].back() = next.back()[i] = v;
        }
        nextActive.push_back(id);
        active = nextActive;
        d = next;
    }
    t.root = active[0];
    return t;
}
inline std::string newick(const Tree &t, int id) {
    const auto &node = t.nodes[id];
    if (node.children.empty()) {
        std::string s = "'";
        for (char c : node.label) {
            s += c;
            if (c == '\'')
                s += '\'';
        }
        return s + "'";
    }
    std::string s = "(";
    for (size_t i = 0; i < node.children.size(); ++i) {
        if (i)
            s += ",";
        auto e = node.children[i];
        s += newick(t, e.child) + ":" + number(e.length);
    }
    return s + ")";
}
inline Matrix treeDistances(const Tree &t, size_t n) {
    std::vector<std::vector<Edge>> adj(t.nodes.size());
    for (size_t i = 0; i < t.nodes.size(); ++i)
        for (auto e : t.nodes[i].children) {
            adj[i].push_back(e);
            adj[e.child].push_back({int(i), e.length});
        }
    Matrix d(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
        std::function<void(int, int, double)> dfs = [&](int u, int parent, double x) {
            if (size_t(u) < n)
                d[i][u] = x;
            for (auto e : adj[u])
                if (e.child != parent)
                    dfs(e.child, u, x + e.length);
        };
        dfs(int(i), -1, 0);
    }
    return d;
}
inline std::string xml(const std::string &s) {
    std::string out;
    for (char c : s) {
        if (c == '&')
            out += "&amp;";
        else if (c == '<')
            out += "&lt;";
        else if (c == '>')
            out += "&gt;";
        else if (c == '\"')
            out += "&quot;";
        else
            out += c;
    }
    return out;
}
inline void svg(const Tree &t, const std::string &path, const std::string &title) {
    std::vector<double> x(t.nodes.size()), y(t.nodes.size());
    size_t leaf = 0;
    double maxX = 0;
    std::function<void(int, int)> layout = [&](int u, int depth) {
        x[u] = 60 + depth * 140;
        maxX = std::max(maxX, x[u]);
        if (t.nodes[u].children.empty())
            y[u] = 110 + 65 * leaf++;
        else {
            double sum = 0;
            for (auto e : t.nodes[u].children) {
                layout(e.child, depth + 1);
                sum += y[e.child];
            }
            y[u] = sum / t.nodes[u].children.size();
        }
    };
    layout(t.root, 0);
    double width = std::max(850.0, maxX + 300), height = 95 + leaf * 65;
    std::ofstream o(path);
    if (!o)
        throw std::runtime_error("No se puede escribir SVG");
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height
      << "\" viewBox=\"0 0 " << width << " " << height
      << "\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><g "
         "font-family=\"Arial,sans-serif\" fill=\"#172b3a\"><text x=\"28\" y=\"32\" "
         "font-size=\"20\">"
      << xml(title)
      << "</text><text x=\"28\" y=\"58\" font-size=\"13\">Cladograma: ramas no a escala; etiquetas "
         "= longitudes calculadas.</text><text x=\"28\" y=\"79\" font-size=\"13\">"
      << (t.rooted
              ? "UPGMA: raiz inferida bajo reloj constante."
              : "NJ: sin raiz biologica; nodo de dibujo en el punto medio de la ultima arista.")
      << "</text>";
    for (size_t u = 0; u < t.nodes.size(); ++u)
        for (auto e : t.nodes[u].children) {
            std::string color = e.length < -EPS ? "#b3261e" : "#167d8d";
            o << "<path d=\"M " << x[u] << " " << y[u] << " V " << y[e.child] << " H " << x[e.child]
              << "\" fill=\"none\" stroke=\"" << color << "\" stroke-width=\"2\"/><text x=\""
              << (x[u] + x[e.child]) / 2 - 28 << "\" y=\"" << y[e.child] - 7
              << "\" font-size=\"14\" fill=\"" << color << "\">" << number(e.length) << "</text>";
        }
    for (size_t u = 0; u < t.nodes.size(); ++u)
        if (t.nodes[u].children.empty())
            o << "<text x=\"" << x[u] + 9 << "\" y=\"" << y[u] + 5 << "\" font-size=\"15\">"
              << xml(t.nodes[u].label) << "</text>";
    o << "</g></svg>";
}
inline void writeMatrix(const Matrix &d, const std::vector<std::string> &labels,
                        const std::string &path) {
    std::ofstream o(path);
    if (!o)
        throw std::runtime_error("No se puede escribir matriz");
    o << "taxon";
    for (auto &s : labels)
        o << "," << s;
    o << "\n";
    for (size_t i = 0; i < d.size(); ++i) {
        o << labels[i];
        for (double x : d[i])
            o << "," << number(x);
        o << "\n";
    }
}
} // namespace phylo
