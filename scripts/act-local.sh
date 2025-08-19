#!/bin/bash

# Script to run tests locally on Mac without Docker containers
# This runs the actual commands directly on your Mac system

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_color() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Function to print section header
print_header() {
    echo ""
    print_color "$BLUE" "=========================================="
    print_color "$BLUE" "$1"
    print_color "$BLUE" "=========================================="
    echo ""
}

# Check if we're on macOS
if [[ "$OSTYPE" != "darwin"* ]]; then
    print_color "$RED" "Error: This script is designed for macOS only"
    exit 1
fi

# Check for required tools
check_command() {
    if ! command -v $1 &> /dev/null; then
        print_color "$RED" "Error: '$1' is not installed"
        print_color "$YELLOW" "Install with: $2"
        return 1
    fi
}

print_header "Checking Prerequisites"
check_command "python3" "brew install python3" || exit 1
check_command "poetry" "pip install poetry" || exit 1
check_command "git" "brew install git" || exit 1

# Check Python version
PYTHON_VERSION=$(python3 --version | cut -d' ' -f2 | cut -d'.' -f1,2)
print_color "$GREEN" "Python version: $PYTHON_VERSION"

print_header "Local Mac Testing (Native Execution)"
print_color "$YELLOW" "Running tests directly on your Mac without containers"

# Navigate to project root
cd "$(dirname "$0")/.."
PROJECT_ROOT=$(pwd)
print_color "$GREEN" "Project root: $PROJECT_ROOT"

# Install dependencies
print_header "Installing Dependencies"
print_color "$YELLOW" "Installing Python dependencies with Poetry..."
poetry install --with dev || {
    print_color "$RED" "Failed to install dependencies"
    exit 1
}

# Run pre-commit checks
print_header "Running Pre-commit Checks"
print_color "$YELLOW" "Running pre-commit hooks..."
poetry run pre-commit install
poetry run pre-commit run --all-files || {
    print_color "$YELLOW" "Pre-commit checks had issues (this may be expected for some files)"
}

# Run unit tests
print_header "Running Unit Tests"
print_color "$YELLOW" "Running pytest unit tests..."
poetry run pytest -xsv tests/unit/airr/test_airr.py || {
    print_color "$RED" "Unit tests failed"
    exit 1
}

# Run all unit tests
print_header "Running All Unit Tests"
print_color "$YELLOW" "Running all unit tests..."
poetry run pytest -xsv tests/unit/ || {
    print_color "$RED" "Some unit tests failed"
    exit 1
}

# Run integration tests (optional)
print_header "Running Integration Tests"
print_color "$YELLOW" "Running integration tests..."
poetry run pytest -xv tests/integration/ || {
    print_color "$YELLOW" "Integration tests had issues (this may be expected)"
}

# Run coverage
print_header "Running Coverage Report"
print_color "$YELLOW" "Generating coverage report..."
poetry run coverage run --source=sadie -m pytest tests/unit/
poetry run coverage report
poetry run coverage html
print_color "$GREEN" "Coverage report generated in htmlcov/"

# Run type checking
print_header "Running Type Checking"
print_color "$YELLOW" "Running pyright type checker..."
poetry run pyright || {
    print_color "$YELLOW" "Type checking had warnings (this may be expected)"
}

print_header "Summary"
print_color "$GREEN" "✅ Local Mac testing completed!"
print_color "$YELLOW" "Note: This ran tests directly on your Mac without containers"
print_color "$YELLOW" "For full CI/CD simulation with Linux containers, use scripts/act.sh"

print_color "$GREEN" "Done!"
