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
- [Solución para Visual Studio 2026](visual_studio/EnsamblajeADN.sln): dos proyectos independientes.
- `capturas/visual_studio/`: capturas reales de Visual Studio y de las ventanas de CMD usadas en el informe.
- `resultados/visual_studio/`: registros de compilación y ejecución, grafos DOT, pruebas y mediciones con MSVC.
- Las capturas y resultados de la primera ejecución con g++ permanecen en las carpetas superiores como antecedentes.

## Compilación y ejecución con Visual Studio 2026

1. Abre `visual_studio/EnsamblajeADN.sln` en Visual Studio 2026 con las herramientas de desarrollo de escritorio con C++ instaladas.
2. Selecciona **Release** y **x64**.
3. Elige **Compilar > Recompilar solución**. Se generan `build/vs2026/Release/consenso.exe` y `hamiltoniano.exe`.
4. Abre `visual_studio/ejecutar_consenso.cmd` o uno de los archivos `ejecutar_hamiltoniano_t2.cmd`, `ejecutar_hamiltoniano_t7.cmd` y `ejecutar_hamiltoniano_t8.cmd`. La ventana muestra la salida y espera una tecla antes de cerrarse.

También puedes compilar abriendo `visual_studio/compilar.cmd`; utiliza la instalación de Visual Studio Community 2026 en su ruta predeterminada. La solución se verificó con Visual Studio 18.8.3, MSVC 19.51.36252 y el conjunto de herramientas v145. Los proyectos seleccionan C++14, la opción más antigua de MSVC, aunque el código es compatible con C++11.

Cada proyecto tiene su propio `main`. Para introducir datos manualmente, ejecuta el archivo `.exe` correspondiente o elige el proyecto como proyecto de inicio y usa **Ctrl+F5**.

## Compilación alternativa con g++

Desde esta carpeta, en PowerShell:

```powershell
New-Item -ItemType Directory -Force build
g++ -std=c++11 -O2 -Wall -Wextra -pedantic codigos/01_consenso.cpp -o build/consenso.exe
g++ -std=c++11 -O2 -Wall -Wextra -pedantic codigos/02_hamiltoniano.cpp -o build/hamiltoniano.exe
```

También puedes abrir cada archivo `.cpp` en un proyecto de consola separado de tu editor. Cada archivo tiene su propio `main`; no deben compilarse juntos en un único ejecutable.

## Ejecución alternativa con g++

En PowerShell:

```powershell
Get-Content datos/consenso.txt | .\build\consenso.exe
Get-Content datos/hamiltoniano_t2.txt | .\build\hamiltoniano.exe
Get-Content datos/hamiltoniano_t7.txt | .\build\hamiltoniano.exe
Get-Content datos/hamiltoniano_t8.txt | .\build\hamiltoniano.exe
```

En CMD también se puede usar la redirección de entrada:

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

Ambos programas se recompilaron desde Visual Studio 2026 y desde CMD con MSBuild: dos proyectos correctos, cero errores y cero advertencias. Se ejecutaron 12 pruebas de funcionamiento: datos oficiales, complementos reversos, solapamientos no máximos, fragmentos contenidos, ausencia de enlaces, un fragmento y entradas inválidas. Los registros actuales se encuentran en `resultados/visual_studio/pruebas`; el resumen está en `resultados/visual_studio/pruebas.txt`.

Se midieron tres ejecuciones por caso con los ejecutables Release x64 de MSVC. `resultados/visual_studio/rendimiento.json` conserva las tres mediciones, su mediana y el número de estados; `resultados/visual_studio/rendimiento.txt` muestra el resumen. Para repetir una medición, usa la entrada correspondiente `datos/prueba_..._consenso.txt` o `datos/prueba_..._hamiltoniano.txt` con su programa y repite la ejecución tres veces. Los tiempos pueden variar entre equipos.

El informe actualizado utiliza capturas reales de CMD alojado en Windows Terminal. Para mostrar la salida completa y legible se usó `type` sobre los registros guardados; solo se omitieron líneas vacías. El consenso se dividió en dos vistas de la misma ejecución. Los archivos `captura_*.cmd` vuelven a mostrar esos registros; los archivos `ejecutar_*.cmd` realizan una ejecución nueva. Las capturas de pruebas y rendimiento muestran resúmenes de ejecuciones reales. La captura del entorno Visual Studio y su registro `compilacion_ide.txt` también se entregan. Los ejecutables y las cachés del entorno quedan excluidos del repositorio y se generan al compilar.
