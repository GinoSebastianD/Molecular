#include "distancias.h"
#include "arboles.h"
#include <sstream>
using namespace std;

int total = 0;
void comprobar(bool correcto, string nombre) {
    if (!correcto) throw runtime_error("FALLO: " + nombre);
    cout << "OK  " << nombre << '\n';
    total++;
}

int main() {
    try {
        ostringstream pasos;
        comprobar(abs(jukesCantor(0.25)-0.304098831081123) < 1e-12, "JC con p=0.25");
        comprobar(jukesCantor(0) == 0, "JC con secuencias iguales");
        bool rechaza = false;
        try { jukesCantor(0.75); } catch (const exception&) { rechaza = true; }
        comprobar(rechaza, "JC rechaza p=0.75");
        Alineamiento a = alinear("ATTGCCATT","ATGGCCATT");
        comprobar(a.diferencias == 1 && a.comparables == 9 && a.huecos == 0, "Conteos de S1 y S2");
        a = alinear("ACGT","AGT");
        comprobar(a.diferencias == 0 && a.comparables == 3 && a.huecos == 1, "Huecos separados de sustituciones");
        rechaza = false;
        try { alinear("ACNX","ACGT"); } catch (const exception&) { rechaza = true; }
        comprobar(rechaza, "Rechazo de bases invalidas");
        vector<string> nombres = {"A","B","C","D"};
        vector<vector<double>> ultra = {{0,2,6,6},{2,0,6,6},{6,6,0,4},{6,6,4,0}};
        vector<Nodo> t = upgma(ultra,nombres,pasos);
        comprobar(errorRMSE(t,ultra) < 1e-12 && abs(t.back().altura-3) < 1e-12, "UPGMA recupera una matriz ultrametrica");
        t = upgma({{0,8,7,12},{8,0,9,14},{7,9,0,11},{12,14,11,0}},nombres,pasos);
        comprobar(abs(t.back().altura-37.0/6) < 1e-12, "UPGMA pondera por cantidad de hojas");
        vector<vector<double>> aditiva = {{0,5,9,9,8},{5,0,10,10,9},{9,10,0,8,7},{9,10,8,0,3},{8,9,7,3,0}};
        t = neighborJoining(aditiva,{"A","B","C","D","E"},pasos);
        comprobar(errorRMSE(t,aditiva) < 1e-12 && abs(t[0].rama-2) < 1e-12 && abs(t[1].rama-3) < 1e-12,
                  "NJ recupera distancias y ramas de una matriz aditiva");
        t = neighborJoining({{0,3},{3,0}},{"A","B"},pasos);
        comprobar(abs(distanciaArbol(t,0,1)-3) < 1e-12, "NJ cierra la ultima arista sin duplicarla");
        t = neighborJoining({{0,1,1},{1,0,3},{1,3,0}},{"A","B","C"},pasos);
        comprobar(abs(t[0].rama+0.5) < 1e-12, "NJ conserva una rama negativa conocida");
        rechaza = false;
        try { upgma({{0,1},{2,0}},{"A","B"},pasos); } catch (const exception&) { rechaza = true; }
        comprobar(rechaza, "Rechazo de matriz asimetrica");
        comprobar(abs(levenshtein("NEYRA","NEIRA")-0.2) < 1e-12, "Levenshtein normalizado");
        comprobar(distanciaApellidos("BECERRA","RODRIGUEZ","RODRIGUEZ","BECERRA") == 0,
                  "Apellidos cruzados");
        cout << total << " pruebas correctas\n";
    } catch (const exception& error) {
        cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
