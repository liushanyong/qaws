# Building Qaws

Qaws has two build systems: CMake (cross-platform) and Sharpmake (Visual Studio).

## CMake

### Basic build

```bash
cd qaws
mkdir build && cd build
cmake ..
cmake --build .
```

### Configuration options

| Option | Default | Description |
|---|---|---|
| `QAWS_SCALAR_IS_FLOAT` | `ON` | Use `float` (f32). Set to `OFF` for `double` (f64). |
| `QAWS_ENABLE_WARNINGS` | `ON` | Enable compiler warnings |
| `QAWS_WARNINGS_AS_ERRORS` | `ON` (f64) / `OFF` (f32) | Treat warnings as errors |
| `QAWS_DISABLE_STRICT_ALIASING` | `ON` | Pass `-fno-strict-aliasing` on GCC/Clang |
| `QAWS_INSTALL` | `ON` | Enable install rules |
| `BUILD_TESTING` | `ON` | Build test suite |

### f32 vs f64 build

```bash
# Float (default)
cmake .. -DQAWS_SCALAR_IS_FLOAT=ON

# Double
cmake .. -DQAWS_SCALAR_IS_FLOAT=OFF
```

### Running tests

```bash
cmake --build . --target check
# or
ctest --output-on-failure
# or run directly
./tests/Debug/qaws_tests      # Windows
./tests/qaws_tests             # Linux/macOS
```

Tests produce SVG files in `working_dir/tests_output/` for visual verification of curve shapes.

### Installing

```bash
cmake --install . --prefix /usr/local
```

Installs:
- Library: `lib/libqaws.a` (or `.so` / `.dll`)
- Headers: `include/qaws*.h`
- CMake config: `lib/cmake/Qaws/`
- pkg-config: `lib/pkgconfig/qaws.pc`

### Using from another CMake project

```cmake
find_package(Qaws REQUIRED)
target_link_libraries(myapp PRIVATE Qaws::qaws)
```

## Sharpmake (Visual Studio)

### Prerequisites

- Visual Studio 2022
- .NET SDK 6.0 or later
- Git submodules initialized

### Bootstrap

```bash
git submodule update --init --recursive
cd buildsystem
bootstrap.bat
```

This will:
1. Build Sharpmake from source
2. Copy Sharpmake binaries to `buildsystem/sharpmake/`
3. Generate Visual Studio solution

### Opening the solution

After bootstrap, open `Qaws_vs2022_win64.sln` in the repo root.

### Regenerating projects

```bash
cd buildsystem
generate_projects.bat
```

## Compiler compatibility

| Compiler | Status |
|---|---|
| MSVC 2022 (cl) | Tested, /W4 clean |
| GCC 11+ | Supported |
| Clang 13+ | Supported |

The library requires C11 (`_Static_assert`, anonymous structs, etc.).
