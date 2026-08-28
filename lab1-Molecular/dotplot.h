#ifndef DOTPLOT_H
#define DOTPLOT_H

#include <string>

void generarDotPlot(const std::string& secuenciaA, const std::string& secuenciaB,
    const std::string& archivoSalida, int ventana = 10, int umbral = 7);

#endif
