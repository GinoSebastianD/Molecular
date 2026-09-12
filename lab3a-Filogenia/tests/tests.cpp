#include "../src/phylogeny.hpp"
#include <iostream>
using namespace phylo;
static int count = 0;
static void check(bool value, const std::string &name) {
    if (!value)
        throw std::runtime_error(name);
    ++count;
    std::cout << "PASS " << name << "\n";
}
static bool near(double a, double b) {
    return std::abs(a - b) < 1e-9;
}
template <class F> void throws(F f, const std::string &name) {
    bool ok = false;
    try {
        f();
    } catch (const std::exception &) {
        ok = true;
    }
    check(ok, name);
}
int main() {
    try {
        check(near(jukesCantor(0), 0), "JC identidad");
        check(near(jukesCantor(0.25), 0.30409883108112323), "JC ejemplo de teoria p=0.25");
        throws([] { jukesCantor(0.75); }, "JC rechaza saturacion");
        throws([] { jukesCantor(-0.1); }, "JC rechaza p negativo");
        auto a = align("ATTGCCATT", "ATGGCCATT");
        check(a.compared == 9 && a.differences == 1 && a.gaps == 0, "Alineamiento S1 S2");
        auto b = align("ACGT", "AGT");
        check(b.a == "ACGT" && b.b == "A-GT" && b.gaps == 1 && b.differences == 0,
              "Alineamiento con delecion");
        std::vector<std::string> names = {"A", "B", "C", "D"};
        Matrix ultra = {{0, 2, 6, 6}, {2, 0, 6, 6}, {6, 6, 0, 4}, {6, 6, 4, 0}};
        auto ut = build(ultra, names, false);
        auto ud = treeDistances(ut, 4);
        bool same = true;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                same &= near(ultra[i][j], ud[i][j]);
        check(same, "UPGMA recupera matriz ultrametrica");
        Matrix classMatrix = {{0, 8, 7, 12}, {8, 0, 9, 14}, {7, 9, 0, 11}, {12, 14, 11, 0}};
        auto ct = build(classMatrix, names, false);
        check(near(ct.nodes[ct.root].height, 37.0 / 6.0),
              "UPGMA pondera grupos desiguales: ejemplo de clase");
        Matrix additive = {{0, 5, 9, 9, 8},
                           {5, 0, 10, 10, 9},
                           {9, 10, 0, 8, 7},
                           {9, 10, 8, 0, 3},
                           {8, 9, 7, 3, 0}};
        auto nt = build(additive, {"A", "B", "C", "D", "E"}, true);
        auto nd = treeDistances(nt, 5);
        same = true;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                same &= near(additive[i][j], nd[i][j]);
        check(same, "NJ recupera distancias de arbol aditivo conocido");
        Matrix two = {{0, 3}, {3, 0}};
        check(near(treeDistances(build(two, {"A", "B"}, true), 2)[0][1], 3),
              "NJ cierre de dos taxones");
        Matrix neg = {{0, 1, 1}, {1, 0, 3}, {1, 3, 0}};
        auto ng = build(neg, {"A", "B", "C"}, true);
        bool negative = false;
        for (auto &n : ng.nodes)
            for (auto e : n.children)
                negative |= e.length < 0;
        check(negative, "NJ conserva ramas negativas para diagnostico");
        Matrix zero(4, std::vector<double>(4));
        auto z1 = build(zero, names, true), z2 = build(zero, names, true);
        check(newick(z1, z1.root) == newick(z2, z2.root), "Desempate determinista");
        throws([] { validate({{0, 1}, {2, 0}}); }, "Rechaza matriz asimetrica");
        throws([] { validate({{0}, {1, 0}}); }, "Rechaza matriz no cuadrada");
        throws([] { validate({{0, -1}, {-1, 0}}); }, "Rechaza distancias negativas de entrada");
        check(normalize("DÍAZ-Núñez") == "diaznunez", "Normalizacion de tildes y enie");
        check(near(editDistance("diaz", "dias"), 0.25), "Levenshtein normalizado");
        Person p = {"P1", "Diaz", "Castro"}, q = {"P2", "Castro", "Diaz"};
        check(near(surnameDistance(p, q), 0), "Apellidos cruzados equivalentes");
        check(near(surnameDistance(p, {"P3", "Rojas", "Silva"}),
                   surnameDistance({"P3", "Rojas", "Silva"}, p)),
              "Distancia apellidos simetrica");
        std::cout << count << " pruebas correctas\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL " << e.what() << "\n";
        return 1;
    }
}
