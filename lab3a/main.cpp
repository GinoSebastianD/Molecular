#include "distancias.h"
#include "salida.h"
#include <filesystem>
using namespace std;

int main(int argc, char** argv) {
    try {
        if (argc > 2) throw runtime_error("Uso: filogenia [directorio_salida]");
        string salida = argc == 2 ? argv[1] : "resultados";
        filesystem::create_directories(salida+"/adn");
        filesystem::create_directories(salida+"/apellidos");
        ostringstream registro;
        registro << "LABORATORIO 03a - UPGMA Y NEIGHBOR JOINING\n";
        vector<string> ids = {"S1","S2","S3","S4","S5"};
        // S4 tiene 8 bases, como indica el PDF de la actividad.
        vector<string> secuencias = {"ATTGCCATT","ATGGCCATT","ATCCAATTTT","ATCTTCTT","ACTGACC"};
        for (int i = 0; i < 5; i++) registro << ids[i] << " = " << secuencias[i] << '\n';
        vector<vector<double>> d(5,vector<double>(5,0)), p = d;
        ofstream alineamientos(salida+"/adn/alineamientos.txt");
        for (int i = 0; i < 5; i++) {
            for (int j = i+1; j < 5; j++) {
                Alineamiento a = alinear(secuencias[i],secuencias[j]);
                if (a.comparables == 0) throw runtime_error("No hay sitios comparables");
                p[i][j] = p[j][i] = double(a.diferencias)/a.comparables;
                d[i][j] = d[j][i] = jukesCantor(p[i][j]);
                if (i != 0 || j != 1) alineamientos << '\n';
                alineamientos << ids[i] << " - " << ids[j] << '\n' << a.a << '\n' << a.b << '\n'
                    << "cambios=" << a.diferencias << " comparables=" << a.comparables
                    << " huecos=" << a.huecos << " p=" << numero(p[i][j])
                    << " JC=" << numero(d[i][j]) << '\n';
            }
        }
        registro << "\nMATRIZ P (sin columnas con huecos)\n";
        matriz(p,ids,registro);
        registro << "\nMATRIZ JUKES-CANTOR\n";
        matriz(d,ids,registro);
        resolver(d,ids,salida+"/adn",registro);

        // Personas ficticias para la prueba de similitud de apellidos.
        vector<string> etiquetas = {"E1","E2","E3","E4","E5","E6","E7","E8"};
        vector<string> nombres = {"Ana","Luis","Carla","Diego","Elena","Pablo","Rosa","Mario"};
        vector<string> paternos = {"BECERRA","RODRIGUEZ","BRICENO","DIAZ","DIAZ","QUISPE","CRUZ","RAMIREZ"};
        vector<string> maternos = {"RODRIGUEZ","BECERRA","QUIROZ","NEYRA","NEIRA","ARROYO","LAURA","URDAY"};
        registro << "\nAPELLIDOS - MUESTRA FICTICIA\n";
        ofstream datos(salida+"/apellidos/personas.csv");
        datos << "id,nombre,paterno,materno\n";
        for (int i = 0; i < 8; i++) {
            registro << etiquetas[i] << " = " << nombres[i] << ' ' << paternos[i] << ' ' << maternos[i] << '\n';
            datos << etiquetas[i] << ',' << nombres[i] << ',' << paternos[i] << ',' << maternos[i] << '\n';
        }
        vector<vector<double>> apellidos(8,vector<double>(8,0));
        for (int i = 0; i < 8; i++)
            for (int j = i+1; j < 8; j++)
                apellidos[i][j] = apellidos[j][i] = distanciaApellidos(paternos[i],maternos[i],paternos[j],maternos[j]);
        matriz(apellidos,etiquetas,registro);
        resolver(apellidos,etiquetas,salida+"/apellidos",registro);
        registro << "\nFin. Matrices, pasos, Newick y SVG guardados en " << salida << '\n';
        cout << registro.str();
        ofstream(salida+"/ejecucion.txt") << registro.str();
    } catch (const exception& error) {
        cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
