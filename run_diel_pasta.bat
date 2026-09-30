@echo off
setlocal EnableDelayedExpansion

echo ============================================================
echo  Diel PASTA -- Launcher  (no admin required)
echo ============================================================

set "SCRIPT_DIR=%~dp0"
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"

REM Virtual environment in user home -- always writable, no admin
set "VENV_DIR=%USERPROFILE%\diel_pasta_venv"
set "VENV_PY=%VENV_DIR%\Scripts\python.exe"

REM ================================================================
REM  STEP 1: If venv already has streamlit AND pyarrow, just launch
REM ================================================================
if exist "%VENV_PY%" (
    "%VENV_PY%" -c "import streamlit, pyarrow" >nul 2>&1
    if !ERRORLEVEL! EQU 0 (
        echo [OK] Virtual env ready at %VENV_DIR%
        goto :launch
    )
    REM venv exists but packages missing -- try to fix without recreating
    echo [INFO] Venv exists but packages incomplete -- reinstalling...
    goto :install_packages
)

REM ================================================================
REM  STEP 2: Find a Python executable
REM ================================================================
set "PYTHON_EXE="

REM Try 'python' on PATH first
python --version >nul 2>&1
if %ERRORLEVEL% EQU 0 ( set "PYTHON_EXE=python" & goto :python_found )

REM Try 'python3'
python3 --version >nul 2>&1
if %ERRORLEVEL% EQU 0 ( set "PYTHON_EXE=python3" & goto :python_found )

REM Prefer x86-64 Anaconda/Miniconda Python -- it has pre-built wheels for
REM everything including pyarrow, avoiding ARM64 compilation failures.
REM ARM64 native Python entries come after as fallback.
for %%P in (
    "C:\ProgramData\Anaconda3\python.exe"
    "C:\ProgramData\miniconda3\python.exe"
    "%USERPROFILE%\Anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
    "%USERPROFILE%\AppData\Local\anaconda3\python.exe"
    "%USERPROFILE%\AppData\Local\miniconda3\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python312\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python311\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python313\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python313-arm64\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python312-arm64\python.exe"
    "C:\Python312\python.exe"
    "C:\Python311\python.exe"
    "C:\Python313\python.exe"
) do (
    if exist %%P ( set "PYTHON_EXE=%%~P" & goto :python_found )
)

echo [ERROR] Cannot find Python. Install Python 3.11+ from https://python.org
echo         and re-run this script.
pause
exit /b 1

:python_found
echo [INFO] Using Python: %PYTHON_EXE%

REM ================================================================
REM  STEP 3: Create virtual environment
REM ================================================================
if not exist "%VENV_DIR%" (
    echo [INFO] Creating virtual environment at %VENV_DIR% ...
    "%PYTHON_EXE%" -m venv "%VENV_DIR%"
    if !ERRORLEVEL! NEQ 0 (
        echo [ERROR] venv creation failed.
        pause
        exit /b 1
    )
)

REM ================================================================
REM  STEP 4: Install packages
REM  pyarrow ARM64 wheels exist for >= 14.0; force binary-only to
REM  prevent pip from falling back to a slow/failing C compilation.
REM ================================================================
:install_packages
echo [INFO] Installing packages (first run takes a few minutes)...
"%VENV_PY%" -m pip install --upgrade pip --quiet

echo [INFO] Installing packages...
"%VENV_PY%" -m pip install ^
    "pyarrow>=14.0" "streamlit>=1.50" "pandas>=1.5" "numpy>=1.23" ^
    "plotly>=5.0" ^
    --prefer-binary --progress-bar on

if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Package install failed. See messages above.
    pause
    exit /b 1
)
echo [OK] Packages installed.

:launch
echo.
echo [INFO] Starting Diel PASTA app...
echo [INFO] Browser opens automatically at http://localhost:8501
echo [INFO] Press Ctrl+C in this window to stop.
echo.
"%VENV_PY%" -m streamlit run "%SCRIPT_DIR%\diel_pasta_app.py"

endlocal
