@echo off
:: Set the repository path (modify this if needed)
set REPO_PATH=D:\Caltech\Your_GitHub_Repo

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

:: Commit with a timestamp message
set DATE_TIME=%date% %time%
echo Committing changes with timestamp...
git commit -m "Auto-sync on %DATE_TIME%"

:: Push changes to GitHub
echo Pushing changes to GitHub...
git push origin main

echo Sync complete!
pause
