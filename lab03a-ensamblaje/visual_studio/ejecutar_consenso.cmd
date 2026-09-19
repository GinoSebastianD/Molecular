@echo off
title Lab03a - consenso - Visual Studio 2026
mode con cols=112 lines=48
cd /d "%~dp0.."
cls
echo VISUAL STUDIO 2026 - EJECUCION EN CMD
echo Gino Sebastian - Laboratorio 03a
echo.
echo build\vs2026\Release\consenso.exe ^< datos\consenso.txt
echo.
build\vs2026\Release\consenso.exe < datos\consenso.txt
echo.
echo Codigo de salida: %errorlevel%
echo.
pause
