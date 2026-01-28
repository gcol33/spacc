# Makefile for spacc C++ tests (Catch2)
# Run from src/tests directory

CXX ?= g++
ifeq ($(OS),Windows_NT)
  CXX = C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe
  TEST_EXE = run_tests.exe
  TEST_EXE_COV = run_tests_cov.exe
else
  TEST_EXE = run_tests
  TEST_EXE_COV = run_tests_cov
endif

CXXFLAGS = -std=c++17 -Wall -Wextra -I../core -I.. -I.

# Coverage flags
COV_FLAGS = --coverage -fprofile-arcs -ftest-coverage -O0 -g

# Test sources
TEST_SOURCES = test_main.cpp \
               test_distance.cpp \
               test_beta.cpp \
               test_hill.cpp \
               test_coverage.cpp \
               test_accumulation.cpp \
               test_balltree.cpp

# Default: build and run tests
all: $(TEST_EXE)
	./$(TEST_EXE)

# Build test executable (optimized)
$(TEST_EXE): $(TEST_SOURCES)
	$(CXX) $(CXXFLAGS) -O2 -o $@ $^

# Build and run with coverage
coverage: $(TEST_EXE_COV)
	./$(TEST_EXE_COV)
	@echo ""
	@echo "=== Coverage Summary ==="
	gcov -r $(TEST_SOURCES) 2>/dev/null | grep -A1 "File.*_core\.h"

$(TEST_EXE_COV): $(TEST_SOURCES)
	$(CXX) $(CXXFLAGS) $(COV_FLAGS) -o $@ $^

# Run specific test file by tag
test-distance: $(TEST_EXE)
	./$(TEST_EXE) "[distance]"

test-beta: $(TEST_EXE)
	./$(TEST_EXE) "[beta]"

test-hill: $(TEST_EXE)
	./$(TEST_EXE) "[hill]"

test-coverage: $(TEST_EXE)
	./$(TEST_EXE) "[coverage]"

test-accumulation: $(TEST_EXE)
	./$(TEST_EXE) "[accumulation]"

test-balltree: $(TEST_EXE)
	./$(TEST_EXE) "[balltree]"

# Verbose output
verbose: $(TEST_EXE)
	./$(TEST_EXE) -s

# List all tests
list: $(TEST_EXE)
	./$(TEST_EXE) --list-tests

# Generate HTML coverage report (requires lcov)
lcov: coverage
	lcov --capture --directory . --output-file coverage.info --rc lcov_branch_coverage=1
	lcov --remove coverage.info '/usr/*' '*/catch.hpp' --output-file coverage.info
	genhtml coverage.info --output-directory coverage_html --branch-coverage
	@echo "Coverage report: coverage_html/index.html"

# Clean build artifacts
clean:
	rm -f $(TEST_EXE) $(TEST_EXE_COV)
	rm -f *.gcov *.gcda *.gcno
	rm -f coverage.info
	rm -rf coverage_html

.PHONY: all coverage test-distance test-beta test-hill test-coverage test-accumulation test-balltree verbose list lcov clean
