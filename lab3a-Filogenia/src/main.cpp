#include "phylogeny.hpp"
#include <filesystem>
#include <iostream>
#include <set>
using namespace phylo;
namespace fs = std::filesystem;

static std::vector<std::string> split(const std::string &s) {
    std::vector<std::string> v;
    std::stringstream in(s);
    std::string x;
    while (std::getline(in, x, ','))
        v.push_back(x);
    return v;
}
static void checkId(const std::string &s, std::set<std::string> &ids) {
    if (s.empty() ||
        s.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-") !=
            std::string::npos ||
        !ids.insert(s).second)
        throw std::invalid_argument(
            "ID invalido o duplicado; use letras ASCII, numeros, guion o guion bajo");
}
static std::vector<std::pair<std::string, std::string>> fasta(const std::string &path) {
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("No se puede abrir FASTA: " + path);
    std::vector<std::pair<std::string, std::string>> records;
    std::set<std::string> ids;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        if (line[0] == '>') {
            std::string id = line.substr(1);
            checkId(id, ids);
            records.push_back({id, ""});
        } else {
            if (records.empty())
                throw std::invalid_argument("FASTA sin encabezado");
            for (char c : line) {
                if (c == ' ' || c == '\t')
                    continue;
                if (c >= 'a' && c <= 'z')
                    c = char(c - 'a' + 'A');
                if (std::string("ACGT").find(c) == std::string::npos)
                    throw std::invalid_argument("FASTA solo admite A,C,G,T");
                records.back().second += c;
            }
        }
    }
    if (records.size() < 2)
        throw std::invalid_argument("FASTA requiere dos o mas secuencias");
    for (auto &r : records)
        if (r.second.empty())
            throw std::invalid_argument("Secuencia vacia");
    return records;
}
static void analyze(const Matrix &d, const std::vector<std::string> &labels, const fs::path &out,
                    const std::string &title) {
    writeMatrix(d, labels, (out / "distancias.csv").string());
    std::ofstream metrics(out / "metricas.csv");
    metrics << "metodo,rmse,mae,error_maximo,ramas_negativas\n";
    for (bool nj : {false, true}) {
        std::string method = nj ? "nj" : "upgma";
        Tree t = build(d, labels, nj);
        auto fitted = treeDistances(t, labels.size());
        double sq = 0, absolute = 0, maxError = 0;
        size_t pairs = 0, negative = 0;
        for (size_t i = 0; i < d.size(); ++i)
            for (size_t j = i + 1; j < d.size(); ++j) {
                double e = std::abs(d[i][j] - fitted[i][j]);
                sq += e * e;
                absolute += e;
                maxError = std::max(e, maxError);
                ++pairs;
            }
        for (auto &node : t.nodes)
            for (auto edge : node.children)
                if (edge.length < -EPS)
                    ++negative;
        metrics << method << "," << number(std::sqrt(sq / pairs)) << "," << number(absolute / pairs)
                << "," << number(maxError) << "," << negative << "\n";
        std::ofstream nw(out / (method + ".nwk"));
        nw << (nj ? "[&U] " : "[&R] ") << newick(t, t.root) << ";\n";
        svg(t, (out / (method + ".svg")).string(),
            title + " / " + (nj ? "Neighbor Joining" : "UPGMA"));
        writeMatrix(fitted, labels, (out / (method + "_distancias_arbol.csv")).string());
        std::ofstream trace(out / (method + "_pasos.txt"));
        for (auto &row : t.trace)
            trace << row << "\n";
        std::cout << method << ": RMSE=" << number(std::sqrt(sq / pairs))
                  << "; ramas negativas=" << negative << "\n";
    }
    size_t violations = 0, triples = 0, triangle = 0, quartets = 0, nonAdditive = 0;
    double maxUltra = 0, maxFour = 0;
    for (size_t i = 0; i < d.size(); ++i)
        for (size_t j = i + 1; j < d.size(); ++j)
            for (size_t k = j + 1; k < d.size(); ++k) {
                std::vector<double> v = {d[i][j], d[i][k], d[j][k]};
                std::sort(v.begin(), v.end());
                ++triples;
                double delta = v[2] - v[1];
                if (delta > 1e-9)
                    ++violations;
                maxUltra = std::max(maxUltra, delta);
                if (v[2] > v[0] + v[1] + 1e-9)
                    ++triangle;
            }
    for (size_t i = 0; i < d.size(); ++i)
        for (size_t j = i + 1; j < d.size(); ++j)
            for (size_t k = j + 1; k < d.size(); ++k)
                for (size_t l = k + 1; l < d.size(); ++l) {
                    std::vector<double> v = {d[i][j] + d[k][l], d[i][k] + d[j][l],
                                             d[i][l] + d[j][k]};
                    std::sort(v.begin(), v.end());
                    ++quartets;
                    double delta = v[2] - v[1];
                    if (delta > 1e-9)
                        ++nonAdditive;
                    maxFour = std::max(maxFour, delta);
                }
    std::ofstream quality(out / "diagnostico.txt");
    quality << "triples=" << triples << "\nviolaciones_ultrametricas=" << violations
            << "\nmax_desvio_ultrametrico=" << number(maxUltra)
            << "\nviolaciones_triangulares=" << triangle << "\ncuartetos=" << quartets
            << "\nviolaciones_cuatro_puntos=" << nonAdditive
            << "\nmax_desvio_cuatro_puntos=" << number(maxFour) << "\n";
}
int main(int argc, char **argv) {
    try {
        if (argc < 4) {
            std::cerr << "Uso: filogenia dna entrada.fasta salida [penalidad_gap]\n     filogenia "
                         "surnames entrada.csv salida\n";
            return 2;
        }
        std::string mode = argv[1];
        fs::path out = argv[3];
        if (mode != "dna" && mode != "surnames")
            throw std::invalid_argument("Modo desconocido");
        if (argc > 5 || (mode == "surnames" && argc != 4))
            throw std::invalid_argument("Numero de argumentos incorrecto");
        fs::create_directories(out);
        std::vector<std::string> labels;
        if (mode == "dna") {
            int gap = -2;
            if (argc == 5) {
                std::string s = argv[4];
                size_t pos = 0;
                gap = std::stoi(s, &pos);
                if (pos != s.size() || gap >= 0 || gap < -100)
                    throw std::invalid_argument("Gap debe ser entero entre -100 y -1");
            }
            auto records = fasta(argv[2]);
            for (auto &r : records)
                labels.push_back(r.first);
            Matrix d(records.size(), std::vector<double>(records.size())), p = d;
            std::ofstream pairs(out / "alineamientos.csv");
            pairs << "par,alineamiento_1,alineamiento_2,score,columnas,comparables,diferencias,"
                     "huecos,p,jc69\n";
            for (size_t i = 0; i < records.size(); ++i)
                for (size_t j = i + 1; j < records.size(); ++j) {
                    auto a = align(records[i].second, records[j].second, gap);
                    if (a.compared == 0)
                        throw std::domain_error("Sin sitios comparables");
                    double proportion = double(a.differences) / a.compared;
                    double corrected = jukesCantor(proportion);
                    p[i][j] = p[j][i] = proportion;
                    d[i][j] = d[j][i] = corrected;
                    pairs << labels[i] << "-" << labels[j] << "," << a.a << "," << a.b << ","
                          << a.score << "," << a.a.size() << "," << a.compared << ","
                          << a.differences << "," << a.gaps << "," << number(proportion) << ","
                          << number(corrected) << "\n";
                }
            std::ofstream settings(out / "parametros.txt");
            settings << "alineamiento=global Needleman-Wunsch por pares\nmatch=2\nmismatch=-1\ngap="
                     << gap << "\nempates=diagonal,arriba,izquierda\nhuecos=exclusion por par\n";
            writeMatrix(p, labels, (out / "distancias_p.csv").string());
            analyze(d, labels, out, "Secuencias ADN");
        } else {
            std::ifstream f(argv[2]);
            if (!f)
                throw std::runtime_error("No se puede abrir CSV de apellidos");
            std::string line;
            std::getline(f, line);
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            if (line != "id,paterno,materno")
                throw std::invalid_argument("CSV requiere encabezado id,paterno,materno");
            std::vector<Person> people;
            std::set<std::string> ids;
            while (std::getline(f, line)) {
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                if (line.empty())
                    continue;
                auto v = split(line);
                if (v.size() != 3)
                    throw std::invalid_argument("CSV requiere tres columnas sin comas internas");
                checkId(v[0], ids);
                normalize(v[1]);
                normalize(v[2]);
                people.push_back({v[0], v[1], v[2]});
                labels.push_back(v[0]);
            }
            if (people.size() < 2)
                throw std::invalid_argument("Se requieren dos o mas personas");
            Matrix d(people.size(), std::vector<double>(people.size()));
            for (size_t i = 0; i < d.size(); ++i)
                for (size_t j = i + 1; j < d.size(); ++j)
                    d[i][j] = d[j][i] = surnameDistance(people[i], people[j]);
            std::ofstream keys(out / "etiquetas.csv");
            keys << "id,paterno,materno\n";
            for (auto &p : people)
                keys << p.id << "," << p.paternal << "," << p.maternal << "\n";
            analyze(d, labels, out, "Similitud de apellidos");
        }
        std::cout << "Resultados guardados en " << out.string() << "\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
