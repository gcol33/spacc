# spacc C++ Refactoring Plan: Decoupling Core Logic from Rcpp

## Goal
Separate pure C++ algorithms from Rcpp bindings to enable:
1. Standalone Catch2 testing of C++ code
2. Proper C++ code coverage metrics
3. Cleaner architecture and potential reuse (Python bindings, etc.)

---

## Phase 1: Create Core Header Structure

### 1.1 Create `src/core/` directory

### 1.2 Extract `core/types.h` - Common type definitions
```cpp
// Pure C++ types used across all core modules
namespace spacc {
  using SiteSpeciesMatrix = std::vector<std::vector<int>>;
  using DistanceMatrix = std::vector<std::vector<double>>;
  using Curve = std::vector<int>;
  using CurveMatrix = std::vector<std::vector<int>>;
  // ... etc
}
```

### 1.3 Extract `core/distance_core.h` - Distance calculations
Functions to extract:
- `euclidean_distance(x1, y1, x2, y2)` → double
- `haversine_distance(lon1, lat1, lon2, lat2)` → double
- `compute_distance_matrix(xs, ys, method)` → DistanceMatrix

### 1.4 Extract `core/beta_core.h` - Beta diversity calculations
Functions to extract:
- `struct BetaComponents { total, turnover, nestedness }`
- `calc_beta_sorensen(a, b, c)` → BetaComponents
- `calc_beta_jaccard(a, b, c)` → BetaComponents
- `count_abc(set1, set2)` → tuple<a,b,c>

### 1.5 Extract `core/hill_core.h` - Hill number calculations
Functions to extract:
- `calc_hill_number(abundances, q)` → double
- `calc_shannon_entropy(abundances)` → double
- `calc_simpson_index(abundances)` → double

### 1.6 Extract `core/coverage_core.h` - Coverage calculations
Functions to extract:
- `calc_chao_coverage(abundances)` → double
- `count_singletons(abundances)` → int
- `count_doubletons(abundances)` → int

### 1.7 Extract `core/accumulation_core.h` - Core accumulation algorithms
Functions to extract:
- `knn_accumulate_single(species_pa, dist_mat, seed)` → Curve
- `find_nearest_unvisited(dist_mat, current, visited)` → int
- `accumulate_species(species_pa, site, seen)` → void

### 1.8 Move existing headers to core/
- `src/balltree.h` → `src/core/balltree.h`
- `src/kdtree_adapter.h` → `src/core/kdtree_adapter.h`

---

## Phase 2: Refactor Rcpp Files to Use Core Headers

### 2.1 Refactor `src/distance.cpp`
- Include `core/distance_core.h`
- Keep only Rcpp exports that wrap core functions
- Convert: `NumericVector` → `std::vector<double>` → call core → convert back

### 2.2 Refactor `src/beta.cpp`
- Include `core/beta_core.h`
- Remove `BetaComponents` struct (now in header)
- Remove `calc_beta_sorensen`, `calc_beta_jaccard` (now in header)
- Keep Rcpp workers that call core functions

### 2.3 Refactor `src/hill.cpp`
- Include `core/hill_core.h`
- Remove pure calculation functions (now in header)
- Keep Rcpp exports

### 2.4 Refactor `src/coverage.cpp`
- Include `core/coverage_core.h`
- Remove pure calculation functions (now in header)
- Keep Rcpp exports

### 2.5 Refactor `src/knn.cpp`
- Include `core/accumulation_core.h`
- Extract core kNN logic
- Keep Rcpp parallel workers

### 2.6 Refactor `src/kncn.cpp`
- Include `core/accumulation_core.h`
- Similar pattern

### 2.7 Refactor `src/random.cpp`
- Include `core/accumulation_core.h`
- Similar pattern

### 2.8 Refactor `src/metrics.cpp`
- Include relevant core headers
- Extract pure metric calculations

### 2.9 Refactor `src/phylo.cpp`
- Include relevant core headers
- Extract phylogenetic calculations

### 2.10 Refactor `src/methods.cpp`
- Include relevant core headers
- Extract any pure C++ logic

---

## Phase 3: Create Catch2 Test Suite

### 3.1 Setup test infrastructure
```
src/tests/
├── catch.hpp           # Already copied
├── test_main.cpp       # Already created
├── Makefile            # Build system
├── run_tests.bat       # Windows runner
└── run_coverage.bat    # Coverage runner
```

### 3.2 Create `test_distance.cpp`
- Test euclidean_distance with known values
- Test haversine_distance with known values
- Test distance matrix computation
- Edge cases: same point, antipodal points

### 3.3 Create `test_beta.cpp`
- Test calc_beta_sorensen with known a,b,c values
- Test calc_beta_jaccard with known a,b,c values
- Verify turnover + nestedness = total
- Edge cases: all shared, no shared, empty sets

### 3.4 Create `test_hill.cpp`
- Test q=0 (species richness)
- Test q=1 (Shannon exponential)
- Test q=2 (inverse Simpson)
- Edge cases: single species, uniform abundances, highly skewed

### 3.5 Create `test_coverage.cpp`
- Test Chao coverage formula
- Test singleton/doubleton counting
- Edge cases: no singletons, all singletons

### 3.6 Create `test_accumulation.cpp`
- Test kNN single curve with simple matrix
- Test species accumulation logic
- Verify monotonic increase

---

## Phase 4: Build System & Coverage

### 4.1 Create `src/tests/Makefile`
```makefile
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -I../core -I.
COV_FLAGS = --coverage -fprofile-arcs -ftest-coverage -O0 -g

TESTS = test_main.cpp test_distance.cpp test_beta.cpp test_hill.cpp \
        test_coverage.cpp test_accumulation.cpp

all: run_tests
    ./run_tests

run_tests: $(TESTS)
    $(CXX) $(CXXFLAGS) -O2 -o $@ $^

coverage: $(TESTS)
    $(CXX) $(CXXFLAGS) $(COV_FLAGS) -o run_tests_cov $^
    ./run_tests_cov
    gcov *.cpp

clean:
    rm -f run_tests run_tests_cov *.gcov *.gcda *.gcno
```

### 4.2 Create `src/tests/run_tests.bat`
```batch
@echo off
set PATH=C:\rtools45\x86_64-w64-mingw32.static.posix\bin;%PATH%
g++ -std=c++17 -Wall -I../core -I. -O2 -o run_tests.exe ^
    test_main.cpp test_distance.cpp test_beta.cpp test_hill.cpp ^
    test_coverage.cpp test_accumulation.cpp
run_tests.exe %*
```

### 4.3 Create `src/tests/run_coverage.bat`
```batch
@echo off
set PATH=C:\rtools45\x86_64-w64-mingw32.static.posix\bin;%PATH%
g++ -std=c++17 -Wall -I../core -I. --coverage -O0 -g -o run_tests_cov.exe ^
    test_main.cpp test_distance.cpp test_beta.cpp test_hill.cpp ^
    test_coverage.cpp test_accumulation.cpp
run_tests_cov.exe
gcov *.cpp 2>nul | findstr /C:"File" /C:"Lines"
```

---

## Phase 5: Verification & Cleanup

### 5.1 Verify R package still builds
```r
devtools::load_all()
devtools::check()
```

### 5.2 Verify all R tests still pass
```r
devtools::test()
```

### 5.3 Run C++ tests
```bash
cd src/tests
./run_tests.bat        # Windows
make all               # Unix
```

### 5.4 Check C++ coverage
```bash
cd src/tests
./run_coverage.bat     # Windows
make coverage          # Unix
```

### 5.5 Update .Rbuildignore
Add `src/tests/` to exclude test files from R package build

### 5.6 Commit and document
- Commit refactored code
- Update CLAUDE.md with new structure
- Document how to run C++ tests

---

## File Checklist

### New files to create:
- [ ] `src/core/types.h`
- [ ] `src/core/distance_core.h`
- [ ] `src/core/beta_core.h`
- [ ] `src/core/hill_core.h`
- [ ] `src/core/coverage_core.h`
- [ ] `src/core/accumulation_core.h`
- [ ] `src/tests/test_distance.cpp`
- [ ] `src/tests/test_beta.cpp`
- [ ] `src/tests/test_hill.cpp`
- [ ] `src/tests/test_coverage.cpp`
- [ ] `src/tests/test_accumulation.cpp`
- [ ] `src/tests/Makefile`
- [ ] `src/tests/run_tests.bat`
- [ ] `src/tests/run_coverage.bat`

### Files to modify:
- [ ] `src/distance.cpp` - use core/distance_core.h
- [ ] `src/beta.cpp` - use core/beta_core.h
- [ ] `src/hill.cpp` - use core/hill_core.h
- [ ] `src/coverage.cpp` - use core/coverage_core.h
- [ ] `src/knn.cpp` - use core/accumulation_core.h
- [ ] `src/kncn.cpp` - use core/accumulation_core.h
- [ ] `src/random.cpp` - use core/accumulation_core.h
- [ ] `src/metrics.cpp` - use core headers
- [ ] `src/phylo.cpp` - use core headers
- [ ] `src/methods.cpp` - use core headers
- [ ] `.Rbuildignore` - exclude src/tests/

### Files to move:
- [ ] `src/balltree.h` → `src/core/balltree.h`
- [ ] `src/kdtree_adapter.h` → `src/core/kdtree_adapter.h`

---

## Expected Outcome

After completion:
1. **C++ tests**: 50+ Catch2 tests covering core algorithms
2. **C++ coverage**: 80%+ on core/ headers
3. **Clean architecture**: Pure C++ separated from Rcpp bindings
4. **Fast iteration**: C++ tests run in <1 second
5. **R package**: Unchanged functionality, all tests pass
