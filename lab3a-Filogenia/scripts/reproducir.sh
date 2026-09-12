#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
mkdir -p build resultados
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic src/main.cpp -o build/filogenia
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic tests/tests.cpp -o build/pruebas
./build/pruebas > resultados/pruebas.txt
cat resultados/pruebas.txt
./build/filogenia dna data/secuencias.fasta resultados/adn
./build/filogenia dna data/secuencias.fasta resultados/adn_gap1 -1
./build/filogenia surnames data/apellidos_demo.csv resultados/apellidos_demo
./build/filogenia surnames "${1:-data/apellidos_clase.csv}" resultados/apellidos_clase
