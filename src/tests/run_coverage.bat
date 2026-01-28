@echo off
REM Run Catch2 tests with coverage for spacc C++ code
REM Run from src/tests directory

set PATH=C:\rtools45\x86_64-w64-mingw32.static.posix\bin;%PATH%

echo Cleaning previous coverage data...
del /q *.gcda *.gcno *.gcov 2>nul

echo Compiling tests with coverage...
g++ -std=c++17 -Wall -Wextra -I../core -I. --coverage -fprofile-arcs -ftest-coverage -O0 -g -o run_tests_cov.exe ^
    test_main.cpp ^
    test_distance.cpp ^
    test_beta.cpp ^
    test_hill.cpp ^
    test_coverage.cpp ^
    test_accumulation.cpp

if %ERRORLEVEL% neq 0 (
    echo Compilation failed!
    exit /b 1
)

echo.
echo Running tests...
.\run_tests_cov.exe

echo.
echo === Coverage Summary ===
gcov -r test_distance.cpp test_beta.cpp test_hill.cpp test_coverage.cpp test_accumulation.cpp 2>nul | findstr /C:"File" /C:"Lines"
