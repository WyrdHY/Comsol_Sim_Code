@echo off
:: Set code page to UTF-8 to handle different characters
chcp 65001 > nul

:: Temporarily set system locale to English
set LC_ALL=en_US.UTF-8
set LANG=en_US.UTF-8

:: Set the repository path (modify this if needed)
set REPO_PATH=D:\Caltech\Comsol_Sim_Code

:: Navigate to the repository directory
cd /d "%REPO_PATH%"

:: Ensure git is available
where git >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Git is not installed or not found in PATH.
    echo Please install Git or add it to your PATH.
    pause
    exit /b
)

:: Pull the latest changes from the remote repository
echo Pulling latest changes from GitHub...
git pull origin main

:: Add all changes
echo Staging all changes...
git add .

:: Get English formatted date and time
for /f "tokens=2 delims==" %%I in ('wmic os get localdatetime /value') do set DATETIME=%%I
set DATE_TIME=%DATETIME:~0,4%-%DATETIME:~4,2%-%DATETIME:~6,2% %DATETIME:~8,2%:%DATETIME:~10,2%:%DATETIME:~12,2%

:: Commit with an English timestamp message
echo Committing changes with timestamp: %DATE_TIME%...
git commit -m "Auto-sync on %DATE_TIME%"

:: Push changes to GitHub
echo Pushing changes to GitHub...
git push origin main

echo Sync complete!
pause
