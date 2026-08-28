
#include <algorithm>
#include "dotplot.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct pares {
    string secuenciaAlineadaA;
    string secuenciaAlineadaB;
};

class algoritmo {
public:
    algoritmo(int coincidencia = 1, int noCoincidencia = -1, int espacio = -2)
        : coincidencia_(coincidencia), noCoincidencia_(noCoincidencia), espacio_(espacio) {
    }

    void computeMatrix(const string& secuenciaA, const string& secuenciaB) {
        secuenciaA_ = secuenciaA;
        secuenciaB_ = secuenciaB;
        int filas = secuenciaA.size(), columnas = secuenciaB.size();
        matrizPuntuacion_.assign(filas + 1, vector<int>(columnas + 1, 0));

        for (int fila = 0; fila <= filas; fila++) {
            matrizPuntuacion_[fila][0] = fila * espacio_;
        }
        for (int columna = 0; columna <= columnas; columna++) {
            matrizPuntuacion_[0][columna] = columna * espacio_;
        }
        for (int fila = 1; fila <= filas; fila++) {
            for (int columna = 1; columna <= columnas; columna++) {
                int diagonal = matrizPuntuacion_[fila - 1][columna - 1] + (secuenciaA[fila - 1] == secuenciaB[columna - 1] ? coincidencia_ : noCoincidencia_);
                int arriba = matrizPuntuacion_[fila - 1][columna] + espacio_;
                int izquierda = matrizPuntuacion_[fila][columna - 1] + espacio_;
                matrizPuntuacion_[fila][columna] = max({ diagonal, arriba, izquierda });
            }
        }
    }

    int score() const { return matrizPuntuacion_[secuenciaA_.size()][secuenciaB_.size()]; }

    vector<pares> allOptimalAlignments(size_t maxResultados = 100000) {
        vector<pares> resultados;
        string alineacionA, alineacionB;
        traceback(secuenciaA_.size(), secuenciaB_.size(), alineacionA, alineacionB, resultados, maxResultados);
        for (auto& resultado : resultados) {
            reverse(resultado.secuenciaAlineadaA.begin(), resultado.secuenciaAlineadaA.end());
            reverse(resultado.secuenciaAlineadaB.begin(), resultado.secuenciaAlineadaB.end());
        }
        return resultados;
    }

    static int countGapBreaks(const pares& alineacion) {
        int rupturas = 0;
        bool enEspacioA = false, enEspacioB = false;
        for (size_t posicion = 0; posicion < alineacion.secuenciaAlineadaA.size(); posicion++) {
            bool espacioA = (alineacion.secuenciaAlineadaA[posicion] == '-');
            bool espacioB = (alineacion.secuenciaAlineadaB[posicion] == '-');
            if (espacioA && !enEspacioA) rupturas++;
            if (espacioB && !enEspacioB) rupturas++;
            enEspacioA = espacioA;
            enEspacioB = espacioB;
        }
        return rupturas;
    }

    static size_t bestByFewestBreaks(const vector<pares>& alineaciones) {
        size_t mejor = 0;
        int menorNumeroRupturas = -1;
        for (size_t posicion = 0; posicion < alineaciones.size(); posicion++) {
            int rupturas = countGapBreaks(alineaciones[posicion]);
            if (menorNumeroRupturas == -1 || rupturas < menorNumeroRupturas) {
                menorNumeroRupturas = rupturas;
                mejor = posicion;
            }
        }
        return mejor;
    }

    const vector<vector<int>>& matrix() const { return matrizPuntuacion_; }

private:
    void traceback(int fila, int columna, string& alineacionA,
        string& alineacionB, vector<pares>& resultados,
        size_t maxResultados) {
        if (resultados.size() >= maxResultados) return;

        if (fila == 0 && columna == 0) {
            resultados.push_back({ alineacionA, alineacionB });
            return;
        }

        if (fila > 0 && columna > 0) {
            int puntuacionDiagonal = matrizPuntuacion_[fila - 1][columna - 1] +
                (secuenciaA_[fila - 1] == secuenciaB_[columna - 1] ? coincidencia_ : noCoincidencia_);
            if (matrizPuntuacion_[fila][columna] == puntuacionDiagonal) {
                alineacionA.push_back(secuenciaA_[fila - 1]);
                alineacionB.push_back(secuenciaB_[columna - 1]);
                traceback(fila - 1, columna - 1, alineacionA, alineacionB, resultados, maxResultados);
                alineacionA.pop_back();
                alineacionB.pop_back();
            }
        }
        if (fila > 0 && matrizPuntuacion_[fila][columna] == matrizPuntuacion_[fila - 1][columna] + espacio_) {
            alineacionA.push_back(secuenciaA_[fila - 1]);
            alineacionB.push_back('-');
            traceback(fila - 1, columna, alineacionA, alineacionB, resultados, maxResultados);
            alineacionA.pop_back();
            alineacionB.pop_back();
        }
        if (columna > 0 && matrizPuntuacion_[fila][columna] == matrizPuntuacion_[fila][columna - 1] + espacio_) {
            alineacionA.push_back('-');
            alineacionB.push_back(secuenciaB_[columna - 1]);
            traceback(fila, columna - 1, alineacionA, alineacionB, resultados, maxResultados);
            alineacionA.pop_back();
            alineacionB.pop_back();
        }
    }

    int coincidencia_, noCoincidencia_, espacio_;
    string secuenciaA_, secuenciaB_;
    vector<vector<int>> matrizPuntuacion_;
};

vector<pair<string, string>> readSequences(
    const string& ruta) {
    ifstream archivo(ruta);
    

    vector<pair<string, string>> secuencias;
    string linea, encabezado, secuenciaActual;
    bool tieneEncabezado = false;
    bool tieneSecuenciaActual = false;

    auto guardarSecuencia = [&]() {
        if (!secuenciaActual.empty()) {
            secuencias.push_back({ tieneEncabezado ? encabezado : ("seq" + to_string(secuencias.size() + 1)),
                             secuenciaActual });
        }
        secuenciaActual.clear();
        tieneEncabezado = false;
        tieneSecuenciaActual = false;
        };

    auto esLineaNumerada = [](const string& texto) {
        for (char caracter : texto) {
            if (isdigit((unsigned char)caracter)) return true;
        }
        return false;
        };

    while (getline(archivo, linea)) {
        while (!linea.empty() && (linea.back() == '\r' || linea.back() == '\n' ||
            isspace((unsigned char)linea.back()))) {
            linea.pop_back();
        }
        if (linea.empty()) continue;

        if (linea[0] == '>') {
            guardarSecuencia();
            encabezado = linea.substr(1);
            tieneEncabezado = true;
        }
        else if (esLineaNumerada(linea)) {
            for (char caracter : linea) {
                if (isalpha((unsigned char)caracter)) secuenciaActual.push_back(toupper(caracter));
            }
            tieneSecuenciaActual = true;
        }
        else {
            if (tieneSecuenciaActual || !secuenciaActual.empty()) {
                guardarSecuencia();
            }
            encabezado = linea;
            tieneEncabezado = true;
        }
    }
    guardarSecuencia();
    return secuencias;
}

void printAlignment(const pares& alineacion, ostream& salida = cout) {
    salida << alineacion.secuenciaAlineadaA << "\n" << alineacion.secuenciaAlineadaB << "\n";
}

void runCase(const string& secuenciaA, const string& secuenciaB,
    const string& nombreA, const string& nombreB,
    bool mostrarTodas = true, size_t maximoParaImprimir = 20) {
    algoritmo alineador;
    alineador.computeMatrix(secuenciaA, secuenciaB);
   
    cout << "Score optimo: " << alineador.score() << "\n";

    auto alineaciones = alineador.allOptimalAlignments();
    cout << "Numero de alineamientos optimos: " << alineaciones.size() << "\n";

    if (mostrarTodas) {
        size_t limite = min(alineaciones.size(), maximoParaImprimir);
        for (size_t posicion = 0; posicion < limite; posicion++) {
            cout << "Solucion " << (posicion + 1) << " ---\n";
            printAlignment(alineaciones[posicion]);
        }
        if (alineaciones.size() > maximoParaImprimir) {
            cout << "... (" << alineaciones.size() - maximoParaImprimir
                << " soluciones mas omitidas)\n";
        }
    }

    if (!alineaciones.empty()) {
        size_t indice = algoritmo::bestByFewestBreaks(alineaciones);
        cout << "Mejor solucion\n";
        printAlignment(alineaciones[indice]);
    }

    if (secuenciaA.size() > 200 || secuenciaB.size() > 200) {
        const string rutaImagen = "dotplot_" + nombreA + "_" + nombreB + ".bmp";
        generarDotPlot(secuenciaA, secuenciaB, rutaImagen);
        cout << "se guardoen: " << rutaImagen << "\n";
    }
    cout << "\n";
}

void handleFileOption(const string& ruta) {
    auto secuencias = readSequences(ruta);
    if (secuencias.empty()) {
        cerr << "No se encontraron secuencias en el archivo\n";
        return;
    }
    cout << "Se leyeron " << secuencias.size() << " secuencias:\n";
    for (auto& [nombre, secuencia] : secuencias) {
        cout << "  " << nombre << ": " << secuencia.size() << " nucleotidos\n";
    }
    cout << "\n";

    if (secuencias.size() >= 3) {
        cout << "Se detectaron >= 3 secuencias. Deseas correr las C(3,2)\n"
            << "combinaciones automaticamente? (s/n): ";
        string resp;
        getline(cin, resp);
        if (!resp.empty() && (resp[0] == 's' || resp[0] == 'S')) {
            for (size_t indiceA = 0; indiceA < secuencias.size(); indiceA++) {
                for (size_t indiceB = indiceA + 1; indiceB < secuencias.size(); indiceB++) {
                    runCase(secuencias[indiceA].second, secuencias[indiceB].second, secuencias[indiceA].first,
                        secuencias[indiceB].first, false);
                }
            }
            return;
        }
    }
    if (secuencias.size() < 2) {
        cerr << "Se necesitan al menos 2 secuencias para alinear\n";
        return;
    }
    runCase(secuencias[0].second, secuencias[1].second, secuencias[0].first, secuencias[1].first);
}

void runMenu() {
    while (true) {
        cout << "1) Ingresar dos cadenas manualmente\n"
            << "2) Leer secuencias desde un archivo\n"
            << "3) Salir\n"
            << "Opcion: ";
        string opcion;
        getline(cin, opcion);

        if (opcion == "1") {
            string secuenciaA, secuenciaB;
            cout << "Cadena A: ";
            getline(cin, secuenciaA);
            cout << "Cadena B: ";
            getline(cin, secuenciaB);
            if (secuenciaA.empty() || secuenciaB.empty()) {
                cout << "Las cadenas no pueden estar vacias.\n\n";
                continue;
            }
            runCase(secuenciaA, secuenciaB, "SeqA", "SeqB");
        }
        else if (opcion == "2") {
            cout << "Ruta del archivo (Sequencias.txt): ";
            string ruta;
            getline(cin, ruta);
            try {
                handleFileOption(ruta);
            }
            catch (const exception& e) {
                cerr << "Error: " << e.what() << "\n\n";
            }
        }
        else if (opcion == "3") {
            break;
        }
        else {
            cout << "Opcion invalida.\n\n";
        }
    }
}

int main(int cantidadArgumentos, char** argumentos) {
   
    vector<string> opciones(argumentos + 1, argumentos + cantidadArgumentos);

    if (opciones.empty()) {
        runMenu();
        return 0;
    }

    if (opciones[0] == "--file") {
        if (opciones.size() < 2) {
            cerr << "Falta el nombre del archivo\n";
            return 1;
        }
        auto secuencias = readSequences(opciones[1]);
        if (secuencias.empty()) {
            cerr << "No se encontraron secuencias en el archivo\n";
            return 1;
        }
        cout << "Se leyeron " << secuencias.size() << " secuencias:\n";
        for (auto& [nombre, secuencia] : secuencias) {
            cout << "  " << nombre << ": " << secuencia.size() << " nucleotidos\n";
        }
        cout << "\n";

        if (secuencias.size() < 2) {
            cerr << "Se necesitan al menos 2 secuencias para alinear\n";
            return 1;
        }

        bool procesarTodas = (opciones.size() >= 3 && opciones[2] == "--all");
        if (procesarTodas && secuencias.size() >= 3) {
            for (size_t indiceA = 0; indiceA < secuencias.size(); indiceA++) {
                for (size_t indiceB = indiceA + 1; indiceB < secuencias.size(); indiceB++) {
                    runCase(secuencias[indiceA].second, secuencias[indiceB].second, secuencias[indiceA].first,
                        secuencias[indiceB].first, false);
                }
            }
        }
        else {
            runCase(secuencias[0].second, secuencias[1].second, secuencias[0].first,
                secuencias[1].first);
        }
    }
    else {
        if (opciones.size() < 2) {
            cerr << "Se requieren dos secuencias\n";
            return 1;
        }
        runCase(opciones[0], opciones[1], "SeqA", "SeqB");
    }

    return 0;
}