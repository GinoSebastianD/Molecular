#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

string leer(string archivo) {
    ifstream entrada(archivo);
    if (!entrada) throw runtime_error("No se pudo leer " + archivo);
    ostringstream texto;
    texto << entrada.rdbuf();
    return texto.str();
}
string escapar(string texto) {
    string r;
    for (char c : texto) {
        if (c == '<') r += "&lt;";
        else if (c == '>') r += "&gt;";
        else if (c == '&') r += "&amp;";
        else r += c;
    }
    return r;
}
void pagina(string archivo, string titulo, string cuerpo) {
    ofstream f("capturas/"+archivo+".html");
    f << "<!doctype html><html lang='es'><meta charset='utf-8'><title>" << titulo << "</title>"
         "<style>body{margin:0;padding:24px;font-family:Arial;background:white;color:#172b3a}"
         "main{width:940px}h1{font-size:25px;margin:0 0 16px}pre{font:18px/1.4 Consolas,monospace;"
         "padding:18px;background:#f3f5f7;white-space:pre-wrap;overflow-wrap:anywhere}"
         "p{font-size:15px;line-height:1.4}img{width:900px;display:block}</style><main><h1>"
      << titulo << "</h1>" << cuerpo << "<p>Laboratorio 03a | Resultados generados por C++ | "
         "Captura del visor de archivos de la ejecucion.</p></main></html>";
}
int main() {
    try {
        filesystem::create_directories("capturas");
        string log = leer("resultados/ejecucion.txt");
        pagina("01-adn","Secuencias y matrices de distancias",
               "<pre>"+escapar(log.substr(0,log.find("\nupgma")))+"</pre>");
        size_t inicio = log.find("APELLIDOS -");
        pagina("04-apellidos","Apellidos y matriz de similitud",
               "<pre>"+escapar(log.substr(inicio,log.find("\nupgma",inicio)-inicio))+"</pre>");
        pagina("07-pruebas","Verificacion de los algoritmos",
               "<pre>"+escapar(leer("resultados/pruebas.txt"))+"</pre>");
        pagina("02-upgma-adn","UPGMA con las cinco secuencias",
               "<img src='../resultados/adn/upgma.svg'><pre>"+escapar(leer("resultados/adn/upgma_pasos.txt"))+"</pre>");
        pagina("03-nj-adn","Neighbor Joining con las cinco secuencias",
               "<img src='../resultados/adn/nj.svg'><pre>"+escapar(leer("resultados/adn/nj_pasos.txt"))+"</pre>");
        pagina("05-upgma-apellidos","UPGMA con apellidos ficticios",
               "<img src='../resultados/apellidos/upgma.svg'>");
        pagina("06-nj-apellidos","Neighbor Joining con apellidos ficticios",
               "<img src='../resultados/apellidos/nj.svg'>");
        cout << "Se generaron siete paginas para capturar los resultados.\n";
    } catch (const exception& error) { cerr << error.what() << '\n'; return 1; }
    return 0;
}
