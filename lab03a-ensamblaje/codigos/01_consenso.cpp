#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <cstdlib>
using namespace std;

vector<vector<string>> fragmentos;
vector<bool> usado;
string mejor, menor;
int objetivo, mejorEnlace = -1;
long long estados = 0, completos = 0;

string complementoReverso(string cadena) {
    reverse(cadena.begin(), cadena.end());
    for (int i = 0; i < (int)cadena.size(); i++) {
        if (cadena[i] == 'A') cadena[i] = 'T';
        else if (cadena[i] == 'T') cadena[i] = 'A';
        else if (cadena[i] == 'C') cadena[i] = 'G';
        else cadena[i] = 'C';
    }
    return cadena;
}

void buscar(string cadena, int cantidad, int enlaceMinimo) {
    estados++;
    if (!mejor.empty()) {
        int limite = objetivo + abs((int)mejor.size() - objetivo);
        if ((int)cadena.size() > limite) return;
    }
    if (cantidad == (int)fragmentos.size()) {
        completos++;
        if (menor.empty() || cadena.size() < menor.size()) menor = cadena;
        int distancia = abs((int)cadena.size() - objetivo);
        int anterior = abs((int)mejor.size() - objetivo);
        if (mejor.empty() || distancia < anterior ||
            (distancia == anterior && cadena.size() < mejor.size()) ||
            (distancia == anterior && cadena.size() == mejor.size() &&
             enlaceMinimo > mejorEnlace)) {
            mejor = cadena;
            mejorEnlace = enlaceMinimo;
        }
        return;
    }
    for (int i = 0; i < (int)fragmentos.size(); i++) {
        if (usado[i]) continue;
        usado[i] = true;
        for (int orientacion = 0; orientacion < 2; orientacion++) {
            string siguiente = fragmentos[i][orientacion];
            if (cadena.empty()) {
                buscar(siguiente, cantidad + 1, 1000000);
            } else if (cadena.find(siguiente) != string::npos) {
                buscar(cadena, cantidad + 1, enlaceMinimo);
            } else {
                int limite = min(cadena.size(), siguiente.size());
                for (int k = limite; k >= 1; k--) {
                    if (cadena.substr(cadena.size() - k) == siguiente.substr(0, k)) {
                        buscar(cadena + siguiente.substr(k), cantidad + 1,
                               min(enlaceMinimo, k));
                    }
                }
            }
        }
        usado[i] = false;
    }
}

void mostrar(string cadena) {
    cout << cadena << "\nLongitud: " << cadena.size() << " bases\n";
    int encontrados = 0;
    for (int i = 0; i < (int)fragmentos.size(); i++) {
        int orientacion = 0;
        size_t posicion = cadena.find(fragmentos[i][0]);
        if (posicion == string::npos) {
            orientacion = 1;
            posicion = cadena.find(fragmentos[i][1]);
        }
        if (posicion == string::npos) continue;
        encontrados++;
        string parte = fragmentos[i][orientacion];
        cout << "f" << i + 1 << (orientacion == 0 ? "+ " : "- ");
        cout << string(posicion, '-') << parte;
        cout << string(cadena.size() - posicion - parte.size(), '-');
        cout << "  [" << posicion + 1 << ", " << posicion + parte.size() << "]\n";
    }
    cout << "Verificacion exacta: " << encontrados << "/";
    cout << fragmentos.size() << " fragmentos\n";
}

int main() {
    int n;
    if (!(cin >> n >> objetivo) || n < 1 || n > 10 || objetivo < 1) {
        cerr << "Entrada invalida: n entre 1 y 10, y longitud positiva.\n";
        return 1;
    }
    fragmentos.resize(n, vector<string>(2));
    usado.assign(n, false);
    for (int i = 0; i < n; i++) {
        if (!(cin >> fragmentos[i][0]) ||
            fragmentos[i][0].find_first_not_of("ACGT") != string::npos) {
            cerr << "Cada fragmento debe contener solamente A, C, G y T.\n";
            return 1;
        }
        fragmentos[i][1] = complementoReverso(fragmentos[i][0]);
    }
    auto inicio = chrono::steady_clock::now();
    buscar("", 0, 1000000);
    auto fin = chrono::steady_clock::now();
    cout << "CONSENSO SIN ERRORES\nFragmentos: " << n;
    cout << "\nLongitud objetivo: " << objetivo << " bases\n";
    cout << "Orientaciones: + directa, - complemento reverso\n";
    if (mejor.empty()) {
        cout << "No existe un ensamblaje conectado con solapamientos positivos.\n";
    } else {
        cout << "\nCONSENSO MAS CERCANO AL OBJETIVO\n";
        mostrar(mejor);
        cout << "Diferencia con el objetivo: " << abs((int)mejor.size() - objetivo);
        cout << " bases\n";
        if (mejorEnlace != 1000000)
            cout << "Menor enlace de la construccion elegida: " << mejorEnlace << " bases\n";
        cout << "\nCONSENSO MAS CORTO\n";
        mostrar(menor);
    }
    cout << "\nEstados visitados: " << estados;
    cout << "\nConstrucciones completas evaluadas: " << completos;
    cout << fixed << setprecision(3) << "\nTiempo de busqueda: ";
    cout << chrono::duration<double, milli>(fin - inicio).count() << " ms\n";
    return 0;
}
