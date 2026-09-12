#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

// Laboratorio 03a. Basado en el codigo proporcionado por el estudiante.
// Compilar: g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic filogenia.cpp -o filogenia
// Ejecutar: ./filogenia [directorio_salida]   o   ./filogenia --test
// Se conservan las fusiones y actualizaciones de UPGMA y NJ del original.
// Ajustes: huecos excluidos al estimar sustituciones; JC sin truncamiento;
// validacion, visualizacion SVG, apellidos y metricas de reconstruccion.

struct Step {
    string left, right;
    double criterion, distance, leftLength, rightLength, height;
    int leftSize, rightSize;
};
vector<Step> upgmaSteps, njSteps;

string lengthText(double value) {
    ostringstream out;
    out << fixed << setprecision(12) << value;
    return out.str();
}

string groupName(const string &newick) {
    string result;
    for (size_t i = 0; i < newick.size(); ++i) {
        if (newick[i] == ':') {
            while (i + 1 < newick.size() && newick[i + 1] != ',' && newick[i + 1] != ')')
                ++i;
        } else
            result += newick[i];
    }
    return result;
}

void validateMatrix(const vector<vector<double>> &d, const vector<string> &names) {
    if (d.size() < 2 || d.size() != names.size())
        throw invalid_argument("Se requieren al menos dos taxones");
    set<string> ids;
    for (const auto &name : names)
        if (name.empty() ||
            name.find_first_not_of(
                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-") !=
                string::npos ||
            !ids.insert(name).second)
            throw invalid_argument("Identificador invalido o repetido");
    for (const auto &row : d)
        if (row.size() != d.size())
            throw invalid_argument("Matriz no cuadrada");
    for (size_t i = 0; i < d.size(); ++i)
        for (size_t j = 0; j < d.size(); ++j)
            if (!isfinite(d[i][j]) || d[i][j] < 0 || abs(d[i][j] - d[j][i]) > 1e-12 ||
                (i == j && abs(d[i][j]) > 1e-12))
                throw invalid_argument("Matriz de distancias invalida");
}

struct Alignment {
    int differences;
    int length; // Columnas base-base comparables, excluidos los huecos.
    int gaps;
    string first, second;
    int score;
};

Alignment needlemanWunschDiff(const string &a, const string &b) {
    if (a.empty() || b.empty())
        throw invalid_argument("Secuencia vacia");
    if (a.find_first_not_of("ACGT") != string::npos || b.find_first_not_of("ACGT") != string::npos)
        throw invalid_argument("Solo se admiten bases A,C,G,T");
    int n = a.size(), m = b.size();
    const int match = 1, mismatch = -1, gap = -2;
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i <= n; i++)
        dp[i][0] = i * gap;
    for (int j = 0; j <= m; j++)
        dp[0][j] = j * gap;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int diag = dp[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? match : mismatch);
            int up = dp[i - 1][j] + gap;
            int left = dp[i][j - 1] + gap;
            dp[i][j] = max({diag, up, left});
        }
    }

    string alignedA, alignedB;
    int i = n, j = m;
    while (i > 0 && j > 0) {
        int score = dp[i][j];
        int diagScore = dp[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? match : mismatch);
        int upScore = dp[i - 1][j] + gap;
        if (score == diagScore) {
            alignedA += a[i - 1];
            alignedB += b[j - 1];
            i--;
            j--;
        } else if (score == upScore) {
            alignedA += a[i - 1];
            alignedB += '-';
            i--;
        } else {
            alignedA += '-';
            alignedB += b[j - 1];
            j--;
        }
    }
    while (i > 0) {
        alignedA += a[--i];
        alignedB += '-';
    }
    while (j > 0) {
        alignedB += b[--j];
        alignedA += '-';
    }

    reverse(alignedA.begin(), alignedA.end());
    reverse(alignedB.begin(), alignedB.end());

    int diffs = 0, comparables = 0, gaps = 0;
    for (size_t k = 0; k < alignedA.size(); k++) {
        if (alignedA[k] == '-' || alignedB[k] == '-') {
            ++gaps;
            continue;
        }
        ++comparables;
        if (alignedA[k] != alignedB[k])
            ++diffs;
    }
    return {diffs, comparables, gaps, alignedA, alignedB, dp[n][m]};
}

double jukesCantor(int differences, int totalLength) {
    if (totalLength <= 0 || differences < 0 || differences > totalLength)
        throw invalid_argument("Conteos invalidos para Jukes-Cantor");
    double p = (double)differences / (double)totalLength;
    if (p >= 0.75)
        throw domain_error("Jukes-Cantor requiere p < 0.75; distancia no estimable");
    return -0.75 * log(1.0 - (4.0 / 3.0) * p);
}

vector<vector<double>> buildDistanceMatrix(const vector<string> &seqs) {
    int n = seqs.size();
    vector<vector<double>> D(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            Alignment aln = needlemanWunschDiff(seqs[i], seqs[j]);
            double d = jukesCantor(aln.differences, aln.length);
            D[i][j] = D[j][i] = d;
        }
    }
    return D;
}

void printMatrix(const vector<vector<double>> &D, const vector<string> &names) {
    cout << "\t";
    for (auto &n : names)
        cout << n << "\t";
    cout << "\n";
    for (size_t i = 0; i < D.size(); i++) {
        cout << names[i] << "\t";
        for (size_t j = 0; j < D.size(); j++)
            cout << D[i][j] << "\t";
        cout << "\n";
    }
}

struct Cluster {
    string newick;
    int size;
    double height;
    bool active;
};

string upgma(vector<vector<double>> D, vector<string> names) {
    validateMatrix(D, names);
    upgmaSteps.clear();
    int n = names.size();
    vector<Cluster> clusters(n);
    for (int i = 0; i < n; i++)
        clusters[i] = {names[i], 1, 0.0, true};

    int remaining = n;
    while (remaining > 1) {
        double best = numeric_limits<double>::max();
        int bi = -1, bj = -1;
        for (int i = 0; i < (int)clusters.size(); i++) {
            if (!clusters[i].active)
                continue;
            for (int j = i + 1; j < (int)clusters.size(); j++) {
                if (!clusters[j].active)
                    continue;
                if (D[i][j] < best) {
                    best = D[i][j];
                    bi = i;
                    bj = j;
                }
            }
        }

        double newHeight = best / 2.0;
        double branchI = newHeight - clusters[bi].height;
        double branchJ = newHeight - clusters[bj].height;

        string newNewick = "(" + clusters[bi].newick + ":" + lengthText(branchI) + "," +
                           clusters[bj].newick + ":" + lengthText(branchJ) + ")";

        int sizeI = clusters[bi].size, sizeJ = clusters[bj].size;
        upgmaSteps.push_back({groupName(clusters[bi].newick), groupName(clusters[bj].newick), best,
                              best, branchI, branchJ, newHeight, sizeI, sizeJ});
        vector<double> newRow(clusters.size(), 0.0);
        for (int k = 0; k < (int)clusters.size(); k++) {
            if (!clusters[k].active || k == bi || k == bj)
                continue;
            double dk = (sizeI * D[bi][k] + sizeJ * D[bj][k]) / (sizeI + sizeJ);
            newRow[k] = dk;
        }

        clusters[bi].active = false;
        clusters[bj].active = false;
        clusters.push_back({newNewick, sizeI + sizeJ, newHeight, true});

        D.push_back(vector<double>(clusters.size(), 0.0));
        for (int k = 0; k < (int)clusters.size() - 1; k++) {
            D[k].push_back(newRow[k]);
            D.back()[k] = newRow[k];
        }

        remaining--;
    }

    for (auto &c : clusters)
        if (c.active)
            return c.newick + ";";
    return "";
}

string neighborJoining(vector<vector<double>> D, vector<string> names) {
    validateMatrix(D, names);
    njSteps.clear();
    vector<string> labels = names;
    vector<int> activeIdx;
    for (size_t i = 0; i < names.size(); i++)
        activeIdx.push_back(i);

    while (activeIdx.size() > 2) {
        int n = activeIdx.size();
        vector<double> r(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                if (i != j)
                    r[i] += D[activeIdx[i]][activeIdx[j]];
        }

        double best = numeric_limits<double>::max();
        int bi = -1, bj = -1;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double q = (n - 2) * D[activeIdx[i]][activeIdx[j]] - r[i] - r[j];
                if (q < best) {
                    best = q;
                    bi = i;
                    bj = j;
                }
            }
        }

        int idxI = activeIdx[bi], idxJ = activeIdx[bj];
        double dij = D[idxI][idxJ];
        double branchI = dij / 2.0 + (r[bi] - r[bj]) / (2.0 * (n - 2));
        double branchJ = dij - branchI;
        njSteps.push_back({groupName(labels[idxI]), groupName(labels[idxJ]), best, dij, branchI,
                           branchJ, 0, 0, 0});

        string newLabel = "(" + labels[idxI] + ":" + lengthText(branchI) + "," + labels[idxJ] +
                          ":" + lengthText(branchJ) + ")";

        int newIdx = D.size();
        for (auto &row : D)
            row.push_back(0.0);
        D.push_back(vector<double>(D.size() + 1, 0.0));

        for (int k = 0; k < n; k++) {
            if (k == bi || k == bj)
                continue;
            int idxK = activeIdx[k];
            double dNew = (D[idxI][idxK] + D[idxJ][idxK] - dij) / 2.0;
            D[newIdx][idxK] = dNew;
            D[idxK][newIdx] = dNew;
        }

        labels.push_back(newLabel);

        vector<int> nextActive;
        for (int k = 0; k < n; k++)
            if (k != bi && k != bj)
                nextActive.push_back(activeIdx[k]);
        nextActive.push_back(newIdx);
        activeIdx = nextActive;
    }

    int a = activeIdx[0], b = activeIdx[1];
    double d = D[a][b];
    njSteps.push_back({groupName(labels[a]), groupName(labels[b]), d, d, d / 2, d / 2, 0, 0, 0});
    return "(" + labels[a] + ":" + lengthText(d / 2.0) + "," + labels[b] + ":" +
           lengthText(d / 2.0) + ");";
}

namespace support {
using Matrix = vector<vector<double>>;
constexpr double EPS = 1e-12;
string number(double x) {
    ostringstream out;
    out << fixed << setprecision(6) << x;
    return out.str();
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

} // namespace support

// Estas funciones leen la salida Newick para dibujar y medir el arbol.
// La topologia siempre proviene de upgma() o neighborJoining(), arriba.
support::Tree parseTree(const string &text, const vector<string> &names, bool rooted) {
    support::Tree t;
    t.rooted = rooted;
    for (auto &name : names)
        t.nodes.push_back({name, {}, 0, 1});
    size_t pos = 0;
    set<int> seen;
    function<int()> parse = [&]() -> int {
        if (pos >= text.size())
            throw invalid_argument("Newick incompleto");
        if (text[pos] == '(') {
            ++pos;
            vector<support::Edge> edges;
            while (true) {
                int child = parse();
                if (pos >= text.size() || text[pos++] != ':')
                    throw invalid_argument("Rama sin longitud");
                size_t used = 0;
                double len = stod(text.substr(pos), &used);
                pos += used;
                edges.push_back({child, len});
                if (pos >= text.size())
                    throw invalid_argument("Newick incompleto");
                if (text[pos] == ')') {
                    ++pos;
                    break;
                }
                if (text[pos++] != ',')
                    throw invalid_argument("Newick invalido");
            }
            int id = int(t.nodes.size());
            t.nodes.push_back({"", edges, 0, 1});
            return id;
        }
        size_t start = pos;
        while (pos < text.size() && text[pos] != ':' && text[pos] != ',' && text[pos] != ')')
            ++pos;
        auto it = find(names.begin(), names.end(), text.substr(start, pos - start));
        if (it == names.end())
            throw invalid_argument("Taxon desconocido");
        int id = int(it - names.begin());
        if (!seen.insert(id).second)
            throw invalid_argument("Hoja repetida");
        return id;
    };
    t.root = parse();
    if (pos + 1 != text.size() || text[pos] != ';' || seen.size() != names.size())
        throw invalid_argument("Newick invalido");
    return t;
}

void saveSteps(const vector<Step> &steps, const filesystem::path &file) {
    ofstream out(file);
    out << "paso,grupo_i,grupo_j,criterio,distancia,rama_i,rama_j,altura,tamano_i,tamano_j\n";
    int i = 0;
    for (const auto &s : steps)
        out << ++i << "," << quoted(s.left) << "," << quoted(s.right) << ","
            << support::number(s.criterion) << "," << support::number(s.distance) << ","
            << support::number(s.leftLength) << "," << support::number(s.rightLength) << ","
            << support::number(s.height) << "," << s.leftSize << "," << s.rightSize << "\n";
}

void analyze(const vector<vector<double>> &d, const vector<string> &names,
             const filesystem::path &out, const string &title) {
    filesystem::create_directories(out);
    support::writeMatrix(d, names, (out / "distancias.csv").string());
    ofstream metrics(out / "metricas.csv");
    metrics << "metodo,rmse,mae,error_maximo,ramas_negativas\n";
    for (bool nj : {false, true}) {
        string method = nj ? "nj" : "upgma";
        string nw = nj ? neighborJoining(d, names) : upgma(d, names);
        auto t = parseTree(nw, names, !nj);
        auto fitted = support::treeDistances(t, names.size());
        double sq = 0, absolute = 0, maxError = 0;
        size_t pairs = 0, negative = 0;
        for (size_t i = 0; i < d.size(); ++i)
            for (size_t j = i + 1; j < d.size(); ++j) {
                double e = abs(d[i][j] - fitted[i][j]);
                sq += e * e;
                absolute += e;
                maxError = max(maxError, e);
                ++pairs;
            }
        for (auto &node : t.nodes)
            for (auto e : node.children)
                if (e.length < -1e-12)
                    ++negative;
        metrics << method << "," << support::number(sqrt(sq / pairs)) << ","
                << support::number(absolute / pairs) << "," << support::number(maxError) << ","
                << negative << "\n";
        ofstream(out / (method + ".nwk")) << (nj ? "[&U] " : "[&R] ") << nw << "\n";
        support::svg(t, (out / (method + ".svg")).string(),
                     title + " / " + (nj ? "Neighbor Joining" : "UPGMA"));
        support::writeMatrix(fitted, names, (out / (method + "_distancias_arbol.csv")).string());
        saveSteps(nj ? njSteps : upgmaSteps, out / (method + "_pasos.csv"));
        cout << method << ": " << nw << "\nRMSE=" << support::number(sqrt(sq / pairs))
             << ", ramas negativas=" << negative << "\n";
    }
    int triples = 0, ultra = 0, triangle = 0, quartets = 0, four = 0;
    double maxUltra = 0, maxFour = 0;
    for (size_t i = 0; i < d.size(); ++i)
        for (size_t j = i + 1; j < d.size(); ++j)
            for (size_t k = j + 1; k < d.size(); ++k) {
                vector<double> v = {d[i][j], d[i][k], d[j][k]};
                sort(v.begin(), v.end());
                ++triples;
                if (v[2] - v[1] > 1e-9)
                    ++ultra;
                if (v[2] - v[1] - v[0] > 1e-9)
                    ++triangle;
                maxUltra = max(maxUltra, v[2] - v[1]);
            }
    for (size_t i = 0; i < d.size(); ++i)
        for (size_t j = i + 1; j < d.size(); ++j)
            for (size_t k = j + 1; k < d.size(); ++k)
                for (size_t l = k + 1; l < d.size(); ++l) {
                    vector<double> v = {d[i][j] + d[k][l], d[i][k] + d[j][l], d[i][l] + d[j][k]};
                    sort(v.begin(), v.end());
                    ++quartets;
                    if (v[2] - v[1] > 1e-9)
                        ++four;
                    maxFour = max(maxFour, v[2] - v[1]);
                }
    ofstream(out / "diagnostico.txt")
        << "triples=" << triples << "\nviolaciones_ultrametricas=" << ultra
        << "\nviolaciones_triangulares=" << triangle << "\ncuartetos=" << quartets
        << "\nviolaciones_cuatro_puntos=" << four
        << "\nmax_desvio_ultrametrico=" << support::number(maxUltra)
        << "\nmax_desvio_cuatro_puntos=" << support::number(maxFour) << "\n";
}

int tests() {
    int count = 0;
    auto check = [&](bool ok, const string &name) {
        if (!ok)
            throw runtime_error("FAIL " + name);
        ++count;
        cout << "PASS " << name << "\n";
    };
    auto reject = [&](function<void()> f, const string &name) {
        bool ok = false;
        try {
            f();
        } catch (const exception &) {
            ok = true;
        }
        check(ok, name);
    };
    auto exact = [&](const vector<vector<double>> &d, const vector<string> &names, bool nj) {
        auto t = parseTree(nj ? neighborJoining(d, names) : upgma(d, names), names, !nj);
        auto fitted = support::treeDistances(t, names.size());
        for (size_t i = 0; i < d.size(); ++i)
            for (size_t j = 0; j < d.size(); ++j)
                if (abs(d[i][j] - fitted[i][j]) > 1e-9)
                    return false;
        return true;
    };
    check(abs(jukesCantor(1, 4) - 0.30409883108112323) < 1e-12, "JC ejemplo de clase");
    check(abs(jukesCantor(0, 4)) < 1e-12, "JC identidad");
    reject([] { jukesCantor(3, 4); }, "JC rechaza saturacion");
    reject([] { jukesCantor(0, 0); }, "JC rechaza denominador cero");
    auto a = needlemanWunschDiff("ATTGCCATT", "ATGGCCATT");
    check(a.length == 9 && a.differences == 1 && a.gaps == 0, "Alineamiento S1-S2");
    auto b = needlemanWunschDiff("ACGT", "AGT");
    check(b.length == 3 && b.differences == 0 && b.gaps == 1, "Huecos excluidos de sustituciones");
    reject([] { needlemanWunschDiff("ACNX", "ACGT"); }, "Rechazo de bases invalidas");
    vector<string> names = {"A", "B", "C", "D"};
    check(exact({{0, 2, 6, 6}, {2, 0, 6, 6}, {6, 6, 0, 4}, {6, 6, 4, 0}}, names, false),
          "UPGMA recupera matriz ultrametrica");
    upgma({{0, 8, 7, 12}, {8, 0, 9, 14}, {7, 9, 0, 11}, {12, 14, 11, 0}}, names);
    check(abs(upgmaSteps.back().height - 37.0 / 6) < 1e-9,
          "UPGMA pondera grupos de tamanos distintos");
    check(exact({{0, 5, 9, 9, 8},
                 {5, 0, 10, 10, 9},
                 {9, 10, 0, 8, 7},
                 {9, 10, 8, 0, 3},
                 {8, 9, 7, 3, 0}},
                {"A", "B", "C", "D", "E"}, true),
          "NJ recupera matriz aditiva");
    check(exact({{0, 3}, {3, 0}}, {"A", "B"}, true), "NJ cierre con dos hojas");
    auto t = parseTree(neighborJoining({{0, 1, 1}, {1, 0, 3}, {1, 3, 0}}, {"A", "B", "C"}),
                       {"A", "B", "C"}, false);
    bool negative = false;
    for (auto &node : t.nodes)
        for (auto e : node.children)
            negative |= e.length < 0;
    check(negative, "NJ conserva ramas negativas");
    reject([] { upgma({{0, 1}, {2, 0}}, {"A", "B"}); }, "Rechazo de matriz asimetrica");
    reject([] { neighborJoining({{0}}, {"A"}); }, "Rechazo de muestra insuficiente");
    check(support::normalize("DÍAZ-Núñez") == "diaznunez", "Normalizacion de apellidos");
    check(abs(support::editDistance("diaz", "dias") - 0.25) < 1e-12, "Levenshtein normalizado");
    check(support::surnameDistance({"A", "Diaz", "Castro"}, {"B", "Castro", "Diaz"}) == 0,
          "Apellidos cruzados");
    cout << count << " pruebas correctas\n";
    return 0;
}

int main(int argc, char **argv) {
    try {
        if (argc == 2 && string(argv[1]) == "--test")
            return tests();
        if (argc > 2)
            throw invalid_argument("Uso: filogenia [directorio_salida] o filogenia --test");
        filesystem::path base = argc == 2 ? argv[1] : "resultados";
        filesystem::create_directories(base / "adn");
        // Secuencias exactas del laboratorio, conservadas del codigo original.
        vector<string> names = {"S1", "S2", "S3", "S4", "S5"};
        vector<string> seqs = {"ATTGCCATT", "ATGGCCATT", "ATCCAATTTT", "ATCTTCTT", "ACTGACC"};
        auto D = buildDistanceMatrix(seqs);
        auto p = D;
        ofstream pairs(base / "adn/alineamientos.csv");
        pairs << "par,alineamiento_1,alineamiento_2,score,columnas,comparables,diferencias,huecos,"
                 "p,jc69\n";
        for (size_t i = 0; i < seqs.size(); ++i)
            for (size_t j = i + 1; j < seqs.size(); ++j) {
                auto a = needlemanWunschDiff(seqs[i], seqs[j]);
                p[i][j] = p[j][i] = double(a.differences) / a.length;
                pairs << names[i] << "-" << names[j] << "," << a.first << "," << a.second << ","
                      << a.score << "," << a.first.size() << "," << a.length << "," << a.differences
                      << "," << a.gaps << "," << support::number(p[i][j]) << ","
                      << support::number(D[i][j]) << "\n";
            }
        support::writeMatrix(p, names, (base / "adn/distancias_p.csv").string());
        cout << "Matriz Jukes-Cantor:\n";
        printMatrix(D, names);
        analyze(D, names, base / "adn", "Secuencias ADN");
        // Muestra real disponible: estudiante y profesor. No se conocen otros alumnos.
        vector<support::Person> real = {{"A1", "Díaz", "Neyra"}, {"P1", "Túpac", "Valdivia"}};
        // Ejemplo ficticio separado: no representa alumnos ni parentescos reales.
        vector<support::Person> demo = {{"D1", "Diaz", "Neira"},  {"D2", "Diaz", "Castro"},
                                        {"D3", "Castro", "Diaz"}, {"D4", "Rojas", "Silva"},
                                        {"D5", "Rojas", "Salas"}, {"D6", "Silva", "Rojas"}};
        for (bool fictional : {false, true}) {
            auto people = fictional ? demo : real;
            string folder = fictional ? "apellidos_demo" : "apellidos_clase";
            vector<string> ids;
            for (auto &person : people)
                ids.push_back(person.id);
            vector<vector<double>> matrix(people.size(), vector<double>(people.size()));
            for (size_t i = 0; i < people.size(); ++i)
                for (size_t j = i + 1; j < people.size(); ++j)
                    matrix[i][j] = matrix[j][i] = support::surnameDistance(people[i], people[j]);
            analyze(matrix, ids, base / folder,
                    fictional ? "Apellidos ficticios" : "Apellidos disponibles");
            ofstream labels(base / folder / "etiquetas.csv");
            labels << "id,paterno,materno\n";
            for (auto &person : people)
                labels << person.id << "," << person.paternal << "," << person.maternal << "\n";
        }
        cout << "Archivos locales generados en " << base.string() << "\n";
        return 0;
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
