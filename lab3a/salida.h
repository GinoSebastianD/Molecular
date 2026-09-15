#ifndef SALIDA_H
#define SALIDA_H
#include "arboles.h"
#include <fstream>
#include <sstream>

string numero(double x) {
    ostringstream s;
    s << fixed << setprecision(6) << x;
    return s.str();
}

string newick(const vector<Nodo>& t, int i) {
    if (t[i].izq == -1) return t[i].nombre;
    int a = t[i].izq, b = t[i].der;
    return "(" + newick(t,a) + ":" + numero(t[a].rama) + "," +
                 newick(t,b) + ":" + numero(t[b].rama) + ")";
}

void matriz(const vector<vector<double>>& d, vector<string> nombres, ostream& s) {
    s << fixed << setprecision(4) << setw(8) << "";
    for (string nombre : nombres) s << setw(8) << nombre;
    s << '\n';
    for (int i = 0; i < (int)d.size(); i++) {
        s << setw(8) << nombres[i];
        for (double valor : d[i]) s << setw(8) << valor;
        s << '\n';
    }
}

void posiciones(const vector<Nodo>& t, int i, int nivel, int& fila,
                vector<double>& x, vector<double>& y) {
    x[i] = 35 + 135*nivel;
    if (t[i].izq == -1) y[i] = 115 + 60*fila++;
    else {
        posiciones(t,t[i].izq,nivel+1,fila,x,y);
        posiciones(t,t[i].der,nivel+1,fila,x,y);
        y[i] = (y[t[i].izq]+y[t[i].der])/2;
    }
}

void guardarSVG(const vector<Nodo>& t, string archivo, string titulo, bool nj) {
    int n = (t.size()+1)/2, fila = 0;
    vector<double> x(t.size()), y(t.size());
    posiciones(t,t.size()-1,0,fila,x,y);
    double ancho = *max_element(x.begin(),x.end())+230;
    ofstream s(archivo);
    s << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 " << ancho << " " << 60*n+145 << "'>\n";
    s << "<rect width='100%' height='100%' fill='white'/>"
         "<g font-family='Arial' fill='#182b3a'>\n";
    s << "<text x='25' y='32' font-size='22' font-weight='bold'>" << titulo << "</text>\n";
    s << "<text x='25' y='59' font-size='14'>Cladograma sin escala; longitudes escritas sobre las ramas.</text>\n";
    s << "<text x='25' y='80' font-size='14'>" << (nj ? "NJ sin raiz evolutiva. Rojo: longitud negativa." : "UPGMA con raiz. Todas las hojas tienen igual altura.") << "</text>\n";
    for (int i = 0; i < (int)t.size()-1; i++) {
        int p = t[i].padre;
        string color = t[i].rama < -1e-12 ? "#b52e32" : "#224e69";
        s << "<path d='M " << x[p] << " " << y[p] << " V " << y[i] << " H " << x[i]
          << "' fill='none' stroke='" << color << "' stroke-width='2'/>\n";
        s << "<text x='" << x[p]+8 << "' y='" << y[i]-8 << "' font-size='14' fill='" << color << "'>" << numero(t[i].rama) << "</text>\n";
    }
    for (int i = 0; i < n; i++)
        s << "<text x='" << x[i]+9 << "' y='" << y[i]+5 << "' font-size='17'>" << t[i].nombre << "</text>\n";
    s << "</g></svg>\n";
}

void resolver(vector<vector<double>> d, vector<string> nombres, string carpeta, ostream& registro) {
    ofstream csv(carpeta+"/distancias.csv");
    csv << "id";
    for (string nombre : nombres) csv << ',' << nombre;
    csv << '\n' << setprecision(12);
    for (int i = 0; i < (int)d.size(); i++) {
        csv << nombres[i];
        for (double valor : d[i]) csv << ',' << valor;
        csv << '\n';
    }
    ofstream metricas(carpeta+"/metricas.csv");
    metricas << "metodo,rmse,ramas_negativas\n";
    for (int metodo = 0; metodo < 2; metodo++) {
        string nombre = metodo == 0 ? "upgma" : "nj";
        ofstream pasos(carpeta+"/"+nombre+"_pasos.txt");
        pasos << fixed << setprecision(6);
        vector<Nodo> t = metodo == 0 ? upgma(d,nombres,pasos) : neighborJoining(d,nombres,pasos);
        int negativas = 0;
        for (Nodo nodo : t) if (nodo.rama < -1e-12) negativas++;
        double rmse = errorRMSE(t,d);
        string arbol = newick(t,t.size()-1)+";";
        ofstream(carpeta+"/"+nombre+".nwk") << (metodo == 0 ? "[&R]" : "[&U]") << arbol << '\n';
        guardarSVG(t,carpeta+"/"+nombre+".svg",metodo == 0 ? "Arbol UPGMA" : "Arbol Neighbor Joining",metodo == 1);
        metricas << nombre << ',' << numero(rmse) << ',' << negativas << '\n';
        registro << "\n" << nombre << "  RMSE=" << numero(rmse) << "  ramas negativas=" << negativas << '\n';
        registro << arbol << '\n';
    }
}
#endif
