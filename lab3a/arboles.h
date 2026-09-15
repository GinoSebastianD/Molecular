#ifndef ARBOLES_H
#define ARBOLES_H
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

struct Nodo {
    string nombre;
    int izq = -1, der = -1, padre = -1, cantidad = 1;
    double altura = 0, rama = 0;
    bool activo = true;
};

void validar(const vector<vector<double>>& d, const vector<string>& nombres) {
    int n = nombres.size();
    if (n < 2 || (int)d.size() != n) throw runtime_error("Matriz de tamano incorrecto");
    for (int i = 0; i < n; i++) {
        if ((int)d[i].size() != n) throw runtime_error("Matriz no cuadrada");
        if (nombres[i].empty() || nombres[i].find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_") != string::npos)
            throw runtime_error("Etiqueta invalida: use letras, numeros o guion bajo");
        for (int j = 0; j < i; j++)
            if (nombres[i] == nombres[j]) throw runtime_error("Etiqueta repetida");
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (!isfinite(d[i][j]) || d[i][j] < 0 || abs(d[i][j]-d[j][i]) > 1e-9 ||
                (i == j && abs(d[i][j]) > 1e-9)) throw runtime_error("Distancia invalida");
}

// Reservamos 2*n-1 posiciones: hojas primero y grupos al final.
vector<Nodo> preparar(vector<vector<double>>& d, vector<string> nombres) {
    validar(d, nombres);
    int n = nombres.size();
    d.resize(2*n-1);
    for (int i = 0; i < 2*n-1; i++) d[i].resize(2*n-1, 0);
    vector<Nodo> t(n);
    for (int i = 0; i < n; i++) t[i].nombre = nombres[i];
    return t;
}

void unir(vector<Nodo>& t, int a, int b, double la, double lb, double altura = 0) {
    Nodo nuevo;
    nuevo.nombre = "(" + t[a].nombre + "," + t[b].nombre + ")";
    nuevo.izq = a; nuevo.der = b; nuevo.altura = altura;
    nuevo.cantidad = t[a].cantidad + t[b].cantidad;
    t[a].padre = t[b].padre = t.size();
    t[a].rama = la; t[b].rama = lb;
    t[a].activo = t[b].activo = false;
    t.push_back(nuevo);
}

vector<Nodo> upgma(vector<vector<double>> d, vector<string> nombres, ostream& pasos) {
    vector<Nodo> t = preparar(d, nombres);
    int n = nombres.size();
    for (int restantes = n; restantes > 1; restantes--) {
        int a = -1, b = -1, u = t.size();
        double menor = INFINITY;
        for (int i = 0; i < u; i++) if (t[i].activo)
            for (int j = i+1; j < u; j++) if (t[j].activo)
                if (d[i][j] < menor - 1e-12) { menor = d[i][j]; a = i; b = j; }
        double h = menor / 2;
        pasos << t[a].nombre << " + " << t[b].nombre << "  d=" << menor << "  h=" << h << '\n';
        for (int k = 0; k < u; k++) if (t[k].activo && k != a && k != b)
            d[u][k] = d[k][u] = (t[a].cantidad*d[a][k] + t[b].cantidad*d[b][k]) /
                                (t[a].cantidad + t[b].cantidad);
        unir(t, a, b, h-t[a].altura, h-t[b].altura, h);
    }
    return t;
}

vector<Nodo> neighborJoining(vector<vector<double>> d, vector<string> nombres, ostream& pasos) {
    vector<Nodo> t = preparar(d, nombres);
    int n = nombres.size();
    for (int restantes = n; restantes > 2; restantes--) {
        int u = t.size(), a = -1, b = -1;
        vector<double> suma(u, 0);
        for (int i = 0; i < u; i++) if (t[i].activo)
            for (int j = 0; j < u; j++) if (t[j].activo) suma[i] += d[i][j];
        double menor = INFINITY;
        for (int i = 0; i < u; i++) if (t[i].activo)
            for (int j = i+1; j < u; j++) if (t[j].activo) {
                double q = (restantes-2)*d[i][j] - suma[i] - suma[j];
                if (q < menor - 1e-12) { menor = q; a = i; b = j; }
            }
        double la = d[a][b]/2 + (suma[a]-suma[b])/(2*(restantes-2));
        double lb = d[a][b] - la;
        pasos << t[a].nombre << " + " << t[b].nombre << "  Q=" << menor
              << "  ramas=" << la << ", " << lb << '\n';
        for (int k = 0; k < u; k++) if (t[k].activo && k != a && k != b)
            d[u][k] = d[k][u] = (d[a][k] + d[b][k] - d[a][b]) / 2;
        unir(t, a, b, la, lb); // Se conservan las longitudes negativas.
    }
    int a = -1, b = -1;
    for (int i = 0; i < (int)t.size(); i++) if (t[i].activo) {
        if (a == -1) a = i; else b = i;
    }
    pasos << "Arista final: " << t[a].nombre << " -- " << t[b].nombre << "  d=" << d[a][b] << '\n';
    // Raiz de dibujo: divide la ultima arista, sin darle sentido ancestral.
    unir(t, a, b, d[a][b]/2, d[a][b]/2);
    return t;
}

double distanciaArbol(const vector<Nodo>& t, int a, int b) {
    vector<bool> visitado(t.size(), false);
    vector<double> desdeA(t.size(), 0);
    double distancia = 0;
    while (a != -1) {
        visitado[a] = true; desdeA[a] = distancia;
        distancia += t[a].rama; a = t[a].padre;
    }
    distancia = 0;
    while (!visitado[b]) { distancia += t[b].rama; b = t[b].padre; }
    return distancia + desdeA[b];
}

double errorRMSE(const vector<Nodo>& t, const vector<vector<double>>& d) {
    double suma = 0;
    int pares = 0;
    for (int i = 0; i < (int)d.size(); i++)
        for (int j = i+1; j < (int)d.size(); j++) {
            double error = distanciaArbol(t,i,j)-d[i][j];
            suma += error*error; pares++;
        }
    return sqrt(suma/pares);
}
#endif
