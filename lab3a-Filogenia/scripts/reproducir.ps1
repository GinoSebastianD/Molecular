param([string]$Apellidos = "data/apellidos_clase.csv")
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path $PSScriptRoot -Parent)
New-Item -ItemType Directory -Force build,resultados | Out-Null
function Ejecutar([string]$Programa, [string[]]$Argumentos) {
    & $Programa @Argumentos
    if ($LASTEXITCODE -ne 0) { throw "Fallo al ejecutar $Programa" }
}
Ejecutar g++ @('-std=c++17','-O2','-Wall','-Wextra','-Wpedantic','-finput-charset=UTF-8','-fexec-charset=UTF-8','src/main.cpp','-o','build/filogenia.exe')
Ejecutar g++ @('-std=c++17','-O2','-Wall','-Wextra','-Wpedantic','-finput-charset=UTF-8','-fexec-charset=UTF-8','tests/tests.cpp','-o','build/pruebas.exe')
& ./build/pruebas.exe | Tee-Object -FilePath resultados/pruebas.txt
if ($LASTEXITCODE -ne 0) { throw 'Pruebas fallidas' }
Ejecutar ./build/filogenia.exe @('dna','data/secuencias.fasta','resultados/adn')
Ejecutar ./build/filogenia.exe @('dna','data/secuencias.fasta','resultados/adn_gap1','-1')
Ejecutar ./build/filogenia.exe @('surnames','data/apellidos_demo.csv','resultados/apellidos_demo')
if ($Apellidos) {
    Ejecutar ./build/filogenia.exe @('surnames',$Apellidos,'resultados/apellidos_clase')
}
