# Laboratorio 03a Ensamblaje de fragmentos de ADN

**Autor:** Gino Sebastian

**Curso:** Biología Molecular Computacional CB309

**Repositorio:** https://github.com/GinoSebastianD/Molecular

El trabajo contiene dos programas independientes en C++11, escritos con funciones, bucles, cadenas, vectores y recursión. Los códigos no tienen comentarios.

## Archivos principales

- [Programa 1: consenso](codigos/01_consenso.cpp).
- [Programa 2: camino hamiltoniano](codigos/02_hamiltoniano.cpp).
- [Informe en PDF con capturas](informe/Informe_Ensamblaje_Gino_Sebastian.pdf).
- `datos/`: las ocho secuencias originales y las entradas usadas para medir rendimiento.
- `capturas/`: capturas de compilación, ejecución, grafo, tiempos y pruebas.
- `resultados/`: registros reales, grafos DOT, figura SVG, mediciones y páginas locales usadas para mostrar los registros en las capturas.

## Compilación

Desde esta carpeta, en PowerShell:

```powershell
New-Item -ItemType Directory -Force build
g++ -std=c++11 -O2 -Wall -Wextra -pedantic codigos/01_consenso.cpp -o build/consenso.exe
g++ -std=c++11 -O2 -Wall -Wextra -pedantic codigos/02_hamiltoniano.cpp -o build/hamiltoniano.exe
```

También puedes abrir cada archivo `.cpp` en un proyecto de consola separado de tu editor. Cada archivo tiene su propio `main`; no deben compilarse juntos en un único ejecutable.

## Ejecución

En PowerShell:

```powershell
Get-Content datos/consenso.txt | .\build\consenso.exe
Get-Content datos/hamiltoniano_t2.txt | .\build\hamiltoniano.exe
Get-Content datos/hamiltoniano_t7.txt | .\build\hamiltoniano.exe
Get-Content datos/hamiltoniano_t8.txt | .\build\hamiltoniano.exe
```

En CMD se puede usar la redirección que aparece en las capturas:

```bat
build\consenso.exe < datos\consenso.txt
build\hamiltoniano.exe < datos\hamiltoniano_t2.txt
```

El programa hamiltoniano escribe `grafo.dot` en la carpeta desde la que se ejecuta y lo reemplaza en cada ejecución. Los resultados originales de cada umbral se conservaron por separado en `resultados/hamiltoniano_t2`, `resultados/hamiltoniano_t7` y `resultados/hamiltoniano_t8`.

Para verlo con Graphviz, si lo tienes instalado:

```text
dot -Tpng grafo.dot -o grafo.png
```

También se entrega la visualización lista en `resultados/grafo_orientaciones_elegidas.svg` y en el PDF.

## Formato de entrada

El primer programa lee `n l` y después `n` cadenas. El segundo lee `n l t` y después `n` cadenas. Se aceptan entre 1 y 10 fragmentos, longitudes objetivo positivas y las letras A, C, G y T en mayúsculas. El umbral `t` debe ser positivo. Se pueden escribir estos datos manualmente al ejecutar el programa.

El consenso explora todos los solapamientos positivos, incluidas coincidencias no máximas, y considera fragmentos contenidos. Busca una cadena conectada sin errores ni huecos. Si no existe una construcción conectada, lo informa. Primero minimiza la distancia a `l`; en empate prefiere menor longitud y luego mayor enlace mínimo. Además, muestra el consenso más corto encontrado.

El programa hamiltoniano calcula el solapamiento máximo de cada par orientado. Elige una orientación por fragmento y lo visita una sola vez. Conserva los enlaces de al menos `t` bases. Entre los caminos completos elige el resultado más cercano a `l`, con menor longitud como desempate.

## Resultados del enunciado

Con `l = 55`, el consenso y el camino con `t = 2` obtienen:

```text
GCCCGCGGCTTCAACGGATCTGTGTCGGGAGTCGCGGGCAGTACTTAACTCGAGG
```

Tiene 55 bases y contiene los ocho fragmentos en alguna de sus dos orientaciones. El camino es:

```text
f1- -> f6- -> f4+ -> f5- -> f8+ -> f3- -> f7- -> f2+
Pesos: 18, 7, 9, 2, 12, 7, 10
```

El consenso más corto y el camino con `t = 7` tienen 50 bases:

```text
CCTCGAGTTAAGTACTGCCCGCGGCTTCAACGGATCTGTGTCGGGAGTCG
```

Con `t = 8` no existe camino hamiltoniano. El grafo completo usa 16 vértices orientados; un camino válido elige 8, uno por fragmento original. Con `t = 2` existen 40 aristas y se enumeran 20 caminos completos, incluyendo caminos equivalentes por complemento reverso. El subgrafo mostrado conserva las ocho orientaciones elegidas y 14 aristas; es acíclico.

## Verificación y tiempos

Ambos programas se compilaron con g++ 10.3.0 sin errores ni advertencias. Se ejecutaron 12 pruebas de funcionamiento: datos oficiales, complementos reversos, solapamientos no máximos, fragmentos contenidos, ausencia de enlaces, un fragmento y entradas inválidas. Los registros se encuentran en `resultados/pruebas`.

Se midieron tres ejecuciones por caso. `resultados/rendimiento.json` conserva las tres mediciones, su mediana y el número de estados; `resultados/rendimiento.txt` muestra el resumen. Para repetir una medición, usa la entrada correspondiente `datos/prueba_..._consenso.txt` o `datos/prueba_..._hamiltoniano.txt` con su programa y repite la ejecución tres veces. Los tiempos pueden variar entre equipos.

Las capturas son imágenes del visor local de los registros reales. Los textos originales también se entregan para comprobar los resultados; las capturas no son imágenes generadas de una consola. Los ejecutables quedan excluidos del repositorio y se generan al compilar.
