# Build and run via -config
echo "Building..."
& 'C:\Program Files\CMake\bin\cmake.exe' -S . -B build -G 'Visual Studio 17 2022' -A x64
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& 'C:\Program Files\CMake\bin\cmake.exe' --build build --config Release
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

$exe = Join-Path 'build/Release' 'stdcellgen.exe'

# Default run: demo.config
& $exe -config data/demo.config
exit $LASTEXITCODE
