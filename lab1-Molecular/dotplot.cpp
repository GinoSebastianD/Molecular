#include "dotplot.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <vector>

void generarDotPlot(const std::string& secuenciaA, const std::string& secuenciaB,
    const std::string& archivoSalida, int ventana, int umbral) {
    if (ventana <= 0 || umbral < 0 || umbral > ventana) {
        throw std::invalid_argument("La ventana y el umbral no son validos");
    }

    const int longitudA = static_cast<int>(secuenciaA.size());
    const int longitudB = static_cast<int>(secuenciaB.size());
    const int maxLado = 500;
    const int pasoA = std::max(1, longitudA / maxLado);
    const int pasoB = std::max(1, longitudB / maxLado);
    const int filas = (longitudA + pasoA - 1) / pasoA;
    const int columnas = (longitudB + pasoB - 1) / pasoB;

    std::vector<std::vector<unsigned char>> imagen(
        filas, std::vector<unsigned char>(columnas, 255));

    for (int fi = 0; fi < filas; ++fi) {
        const int posicionA = fi * pasoA;
        if (posicionA + ventana > longitudA) continue;

        for (int fj = 0; fj < columnas; ++fj) {
            const int posicionB = fj * pasoB;
            if (posicionB + ventana > longitudB) continue;

            int coincidencias = 0;
            for (int desplazamiento = 0; desplazamiento < ventana; ++desplazamiento) {
                if (secuenciaA[posicionA + desplazamiento] == secuenciaB[posicionB + desplazamiento]) ++coincidencias;
            }

            if (coincidencias >= umbral) {
                imagen[fi][fj] = static_cast<unsigned char>(
                    255 - (coincidencias * 255 / ventana));
            }
        }
    }

    const int tamanoFila = (columnas * 3 + 3) & ~3;
    const std::uint32_t tamanoDatosPixeles = static_cast<std::uint32_t>(tamanoFila * filas);
    const std::uint32_t tamanoArchivo = 54 + tamanoDatosPixeles;
    const char fondoRojo = 18;
    const char fondoVerde = 25;
    const char fondoAzul = 38;
    std::ofstream out(archivoSalida, std::ios::binary);
    if (!out) {
        throw std::runtime_error("No se pudo crear la imagen " + archivoSalida);
    }

    auto escribir16 = [&out](std::uint16_t valor) {
        out.put(static_cast<char>(valor & 0xff));
        out.put(static_cast<char>((valor >> 8) & 0xff));
    };
    auto escribir32 = [&out](std::uint32_t valor) {
        for (int desplazamiento = 0; desplazamiento < 32; desplazamiento += 8) {
            out.put(static_cast<char>((valor >> desplazamiento) & 0xff));
        }
    };

    escribir16(0x4d42);
    escribir32(tamanoArchivo);
    escribir16(0);
    escribir16(0);
    escribir32(54);
    escribir32(40);
    escribir32(static_cast<std::uint32_t>(columnas));
    escribir32(static_cast<std::uint32_t>(filas));
    escribir16(1);
    escribir16(24);
    escribir32(0);
    escribir32(tamanoDatosPixeles);
    escribir32(2835);
    escribir32(2835);
    escribir32(0);
    escribir32(0);

    std::vector<char> fila(tamanoFila, fondoAzul);
    for (int fi = filas - 1; fi >= 0; --fi) {
        std::fill(fila.begin(), fila.end(), fondoAzul);
        for (int fj = 0; fj < columnas; ++fj) {
            const int valorPixel = imagen[fi][fj];
            if (valorPixel < 255) {
                const int fuerzaCoincidencia = 255 - valorPixel;
                const char rojo = static_cast<char>(
                    60 + fuerzaCoincidencia * 195 / 255);
                const char verde = static_cast<char>(
                    170 - fuerzaCoincidencia * 120 / 255);
                const char azul = static_cast<char>(
                    255 - fuerzaCoincidencia * 255 / 255);
                fila[fj * 3] = azul;
                fila[fj * 3 + 1] = verde;
                fila[fj * 3 + 2] = rojo;
            }
            else {
                fila[fj * 3] = fondoAzul;
                fila[fj * 3 + 1] = fondoVerde;
                fila[fj * 3 + 2] = fondoRojo;
            }
        }
        out.write(fila.data(), fila.size());
    }
}
