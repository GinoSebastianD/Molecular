#ifndef DISTANCIAS_H
#define DISTANCIAS_H
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

struct Alineamiento {
    string a, b;
    int diferencias = 0, comparables = 0, huecos = 0;
};

// Needleman-Wunsch: coincidencia +1, cambio -1 y hueco -2.
Alineamiento alinear(string a, string b) {
    if (a.empty() || b.empty() || a.find_first_not_of("ACGT") != string::npos ||
        b.find_first_not_of("ACGT") != string::npos)
        throw runtime_error("Se necesitan secuencias no vacias con bases A C G T");
    int n = a.size(), m = b.size();
    vector<vector<int>> tabla(n + 1, vector<int>(m + 1));
    for (int i = 0; i <= n; i++) tabla[i][0] = -2 * i;
    for (int j = 0; j <= m; j++) tabla[0][j] = -2 * j;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int diagonal = tabla[i-1][j-1] + (a[i-1] == b[j-1] ? 1 : -1);
            tabla[i][j] = max(diagonal, max(tabla[i-1][j]-2, tabla[i][j-1]-2));
        }
    }
    Alineamiento r;
    int i = n, j = m;
    // En empates: diagonal, arriba, izquierda.
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && tabla[i][j] == tabla[i-1][j-1] + (a[i-1] == b[j-1] ? 1 : -1)) {
            r.a += a[--i]; r.b += b[--j];
        } else if (i > 0 && tabla[i][j] == tabla[i-1][j] - 2) {
            r.a += a[--i]; r.b += '-';
        } else {
            r.a += '-'; r.b += b[--j];
        }
    }
    reverse(r.a.begin(), r.a.end());
    reverse(r.b.begin(), r.b.end());
    for (int k = 0; k < (int)r.a.size(); k++) {
        if (r.a[k] == '-' || r.b[k] == '-') r.huecos++;
        else {
            r.comparables++;
            if (r.a[k] != r.b[k]) r.diferencias++;
        }
    }
    return r;
}

double jukesCantor(double p) {
    if (!isfinite(p) || p < 0 || p >= 0.75)
        throw runtime_error("Jukes-Cantor requiere 0 <= p < 0.75");
    return p == 0 ? 0 : -0.75 * log(1.0 - 4.0 * p / 3.0);
}

// Los datos de apellidos se escriben en mayusculas ASCII, sin tildes.
double levenshtein(string a, string b) {
    int n = a.size(), m = b.size();
    if (n == 0 || m == 0) throw runtime_error("Apellido vacio");
    vector<vector<int>> tabla(n + 1, vector<int>(m + 1));
    for (int i = 0; i <= n; i++) tabla[i][0] = i;
    for (int j = 0; j <= m; j++) tabla[0][j] = j;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            tabla[i][j] = min(tabla[i-1][j-1] + (a[i-1] != b[j-1]),
                              min(tabla[i-1][j] + 1, tabla[i][j-1] + 1));
    return double(tabla[n][m]) / max(n, m);
}

double distanciaApellidos(string a, string b, string c, string d) {
    double directo = (levenshtein(a,c) + levenshtein(b,d)) / 2;
    double cruzado = (levenshtein(a,d) + levenshtein(b,c)) / 2;
    return min(directo, cruzado);
}
#endif
