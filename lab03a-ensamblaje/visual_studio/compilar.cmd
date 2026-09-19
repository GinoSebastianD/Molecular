@echo off
title Lab03a - Compilacion - Visual Studio 2026
mode con cols=112 lines=48
cd /d "%~dp0.."
call "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64
cls
echo VISUAL STUDIO COMMUNITY 2026 - COMPILACION C++
echo Gino Sebastian - Laboratorio 03a
echo.
echo Compilador de Microsoft instalado con Visual Studio 2026
cl
echo.
echo MSBuild visual_studio\EnsamblajeADN.sln /t:Rebuild /p:Configuration=Release /p:Platform=x64
MSBuild visual_studio\EnsamblajeADN.sln /t:Rebuild /p:Configuration=Release /p:Platform=x64 /nologo /v:minimal /clp:Summary /fl /flp:logfile=resultados\visual_studio\compilacion.log
echo.
echo Codigo de salida de MSBuild: %errorlevel%
echo.
pause
