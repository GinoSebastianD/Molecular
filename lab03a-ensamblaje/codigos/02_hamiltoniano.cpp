#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <cstdlib>
using namespace std;

vector<string> cadenas;
vector<vector<int>> enlace;
vector<bool> usado;
vector<int> camino, mejorCamino;
int n, objetivo, t, mejorLongitud = 0;
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

int solapamiento(string a, string b) {
    int limite = (int)min(a.size(), b.size());
    for (int k = limite; k >= 1; k--) {
        if (a.substr(a.size() - k) == b.substr(0, k)) return k;
    }
    return 0;
}

string nombre(int vertice) {
    return "f" + to_string(vertice / 2 + 1) + (vertice % 2 == 0 ? "+" : "-");
}

void buscar(int longitud) {
    estados++;
    if ((int)camino.size() == n) {
        completos++;
        int diferencia = abs(longitud - objetivo);
        int anterior = abs(mejorLongitud - objetivo);
        if (mejorCamino.empty() || diferencia < anterior ||
            (diferencia == anterior && longitud < mejorLongitud)) {
            mejorCamino = camino;
            mejorLongitud = longitud;
        }
        return;
    }
    for (int siguiente = 0; siguiente < 2 * n; siguiente++) {
        int fragmento = siguiente / 2;
        if (usado[fragmento]) continue;
        int k = 0;
        if (!camino.empty()) {
            k = enlace[camino.back()][siguiente];
            if (k < t) continue;
        }
        usado[fragmento] = true;
        camino.push_back(siguiente);
        buscar(longitud + (int)cadenas[siguiente].size() - k);
        camino.pop_back();
        usado[fragmento] = false;
    }
}

bool guardarGrafo() {
    ofstream archivo("grafo.dot");
    if (!archivo) return false;
    archivo << "digraph Ensamblaje {\n  rankdir=LR;\n";
    for (int i = 0; i < 2 * n; i++) {
        archivo << "  v" << i << " [label=\"" << nombre(i) << "\"];\n";
    }
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2 * n; j++) {
            if (i / 2 == j / 2 || enlace[i][j] < t) continue;
            bool elegido = false;
            for (int k = 1; k < (int)mejorCamino.size(); k++) {
                if (mejorCamino[k - 1] == i && mejorCamino[k] == j) elegido = true;
            }
            archivo << "  v" << i << " -> v" << j;
            archivo << " [label=\"" << enlace[i][j] << "\"";
            if (elegido) archivo << ", color=blue, penwidth=3";
            archivo << "];\n";
        }
    }
    archivo << "}\n";
    return true;
}

int main() {
    if (!(cin >> n >> objetivo >> t) || n < 1 || n > 10 || objetivo < 1 || t < 1) {
        cerr << "Entrada invalida: n entre 1 y 10; longitud y t positivos.\n";
        return 1;
    }
    for (int i = 0; i < n; i++) {
        string cadena;
        if (!(cin >> cadena) || cadena.find_first_not_of("ACGT") != string::npos) {
            cerr << "Cada fragmento debe contener solamente A, C, G y T.\n";
            return 1;
        }
        cadenas.push_back(cadena);
        cadenas.push_back(complementoReverso(cadena));
    }
    usado.assign(n, false);
    enlace.assign(2 * n, vector<int>(2 * n, 0));
    auto inicio = chrono::steady_clock::now();
    int aristas = 0;
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < 2 * n; j++) {
            if (i / 2 == j / 2) continue;
            enlace[i][j] = solapamiento(cadenas[i], cadenas[j]);
            if (enlace[i][j] >= t) aristas++;
        }
    }
    buscar(0);
    auto fin = chrono::steady_clock::now();
    cout << "CAMINO HAMILTONIANO\nFragmentos: " << n;
    cout << "\nLongitud objetivo: " << objetivo << " bases\nLinkage t: " << t;
    cout << "\nVertices orientados: " << 2 * n << "\nAristas: " << aristas;
    cout << "\nCaminos completos: " << completos << "\n";
    if (mejorCamino.empty()) {
        cout << "No existe camino hamiltoniano para este valor de t.\n";
    } else {
        string resultado = cadenas[mejorCamino[0]];
        vector<int> posicion(n, 0);
        cout << "\nCamino: " << nombre(mejorCamino[0]);
        for (int i = 1; i < n; i++) {
            int a = mejorCamino[i - 1], b = mejorCamino[i];
            int k = enlace[a][b];
            posicion[i] = (int)resultado.size() - k;
            resultado += cadenas[b].substr(k);
            cout << " -[" << k << "]-> " << nombre(b);
        }
        cout << "\n\nSECUENCIAS ENLAZADAS\n";
        for (int i = 0; i < n; i++) {
            int v = mejorCamino[i];
            cout << nombre(v) << " " << string(posicion[i], '-') << cadenas[v];
            cout << string(resultado.size() - posicion[i] - cadenas[v].size(), '-');
            cout << "  [" << posicion[i] + 1 << ", ";
            cout << posicion[i] + cadenas[v].size() << "]\n";
        }
        cout << "ADN " << resultado << "\nLongitud: " << resultado.size() << " bases";
        cout << "\nDiferencia con el objetivo: " << abs((int)resultado.size() - objetivo);
        int verificados = 0;
        for (int i = 0; i < n; i++) {
            if (resultado.find(cadenas[2 * i]) != string::npos ||
                resultado.find(cadenas[2 * i + 1]) != string::npos) verificados++;
        }
        cout << " bases\nVerificacion exacta: " << verificados << "/" << n << " fragmentos\n";
    }
    cout << "\nEstados visitados: " << estados;
    cout << fixed << setprecision(3) << "\nTiempo de grafo y busqueda: ";
    cout << chrono::duration<double, milli>(fin - inicio).count() << " ms\n";
    if (!guardarGrafo()) {
        cerr << "No se pudo guardar grafo.dot.\n";
        return 1;
    }
    cout << "Grafo guardado en grafo.dot\n";
    return 0;
}
