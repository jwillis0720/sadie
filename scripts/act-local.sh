#!/bin/bash

# Script to run local-friendly GitHub Actions workflows using act
# These workflows are modified to work without GitHub-specific actions

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
CONTAINER_ARCH="--container-architecture linux/amd64"
PLATFORM_MAP="-P ubuntu-latest=catthehacker/ubuntu:act-latest"
WORKFLOWS_DIR=".github/workflows"

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

# Check if act is installed
if ! command -v act &> /dev/null; then
    print_color "$RED" "Error: 'act' is not installed"
    print_color "$YELLOW" "Install with: brew install act"
    exit 1
fi

print_header "GitHub Actions Local Runner (act)"
print_color "$YELLOW" "Running local-friendly workflows without GitHub Actions dependencies"

# Test pytest workflow
print_header "Running pytest-local.yml"
print_color "$GREEN" "Testing Python 3.13..."
act workflow_dispatch -W "$WORKFLOWS_DIR/pytest-local.yml" -j pytest \
    --matrix python-version:3.13 \
    $CONTAINER_ARCH $PLATFORM_MAP \
    --privileged \
    || print_color "$RED" "pytest-local.yml failed for Python 3.13"

# Test coverage workflow
print_header "Running coverage-local.yml"
act workflow_dispatch -W "$WORKFLOWS_DIR/coverage-local.yml" -j coverage \
    $CONTAINER_ARCH $PLATFORM_MAP \
    --privileged \
    || print_color "$RED" "coverage-local.yml failed"

print_header "Summary"
print_color "$GREEN" "✅ Local workflow tests completed!"
print_color "$YELLOW" "Note: These are simplified versions for local testing"
print_color "$YELLOW" "The actual GitHub workflows may have additional features"

print_color "$GREEN" "Done!"
