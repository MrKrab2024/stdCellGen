param(
  [string]$BuildDir = "build",
  [switch]$Clean
)

if ($Clean -and (Test-Path $BuildDir)) { Remove-Item -Recurse -Force $BuildDir }
New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null

# Configure
cmake -S . -B $BuildDir -DCMAKE_BUILD_TYPE=Release
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

# Build (handle MSVC multi-config)
cmake --build $BuildDir --config Release
exit $LASTEXITCODE
