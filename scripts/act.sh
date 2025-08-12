#!/bin/bash

# Script to run all GitHub Actions workflows locally using act
# Requires: act (brew install act)

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
CONTAINER_ARCH="--container-architecture linux/amd64"
PLATFORM_MAP="-P ubuntu-latest=catthehacker/ubuntu:act-latest -P macos-latest=catthehacker/ubuntu:act-latest"
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

# Create event files directory
EVENTS_DIR="scripts/act-events"
mkdir -p "$EVENTS_DIR"

# Create event file for merged PR
cat > "$EVENTS_DIR/merged-pr.json" << 'EOF'
{
  "action": "closed",
  "pull_request": {
    "merged": true,
    "merge_commit_sha": "abc123",
    "base": {
      "ref": "main"
    },
    "head": {
      "ref": "feature-branch"
    }
  },
  "repository": {
    "default_branch": "main"
  }
}
EOF

# Create event file for push with tag
cat > "$EVENTS_DIR/push-tag.json" << 'EOF'
{
  "ref": "refs/tags/v1.0.0",
  "ref_type": "tag",
  "repository": {
    "default_branch": "main"
  }
}
EOF

# Create event file for push to main
cat > "$EVENTS_DIR/push-main.json" << 'EOF'
{
  "ref": "refs/heads/main",
  "repository": {
    "default_branch": "main"
  }
}
EOF

# Create event file for pull request
cat > "$EVENTS_DIR/pull-request.json" << 'EOF'
{
  "action": "opened",
  "pull_request": {
    "base": {
      "ref": "main"
    },
    "head": {
      "ref": "feature-branch"
    }
  },
  "repository": {
    "default_branch": "main"
  }
}
EOF

# Function to run a workflow
run_workflow() {
    local workflow=$1
    local event=$2
    local event_file=$3
    local job=$4

    if [ -n "$job" ]; then
        print_color "$YELLOW" "Running job: $job from $workflow with event: $event"
        if [ -n "$event_file" ]; then
            act $event -j $job -W "$WORKFLOWS_DIR/$workflow" --eventpath "$event_file" $CONTAINER_ARCH $PLATFORM_MAP $DRY_RUN
        else
            act $event -j $job -W "$WORKFLOWS_DIR/$workflow" $CONTAINER_ARCH $PLATFORM_MAP $DRY_RUN
        fi
    else
        print_color "$YELLOW" "Running all jobs from $workflow with event: $event"
        if [ -n "$event_file" ]; then
            act $event -W "$WORKFLOWS_DIR/$workflow" --eventpath "$event_file" $CONTAINER_ARCH $PLATFORM_MAP $DRY_RUN
        else
            act $event -W "$WORKFLOWS_DIR/$workflow" $CONTAINER_ARCH $PLATFORM_MAP $DRY_RUN
        fi
    fi
}

# Parse command line arguments
DRY_RUN=""
SPECIFIC_WORKFLOW=""
SPECIFIC_JOB=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN="--dryrun"
            print_color "$YELLOW" "Running in DRY RUN mode"
            shift
            ;;
        --workflow)
            SPECIFIC_WORKFLOW="$2"
            shift 2
            ;;
        --job)
            SPECIFIC_JOB="$2"
            shift 2
            ;;
        --help)
            print_header "Act Runner Script - Help"
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --dry-run           Run in dry-run mode (no actual execution)"
            echo "  --workflow NAME     Run specific workflow file"
            echo "  --job NAME          Run specific job"
            echo "  --help              Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                                    # Run all workflows"
            echo "  $0 --dry-run                          # Dry run all workflows"
            echo "  $0 --workflow pytest.yml              # Run only pytest workflow"
            echo "  $0 --workflow pypi.yml --job pypi     # Run specific job from workflow"
            exit 0
            ;;
        *)
            print_color "$RED" "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

print_header "GitHub Actions Local Runner (act)"

# If specific workflow is requested
if [ -n "$SPECIFIC_WORKFLOW" ]; then
    if [ ! -f "$WORKFLOWS_DIR/$SPECIFIC_WORKFLOW" ]; then
        print_color "$RED" "Error: Workflow file not found: $WORKFLOWS_DIR/$SPECIFIC_WORKFLOW"
        exit 1
    fi

    print_header "Running workflow: $SPECIFIC_WORKFLOW"

    # Determine the appropriate event based on workflow
    case "$SPECIFIC_WORKFLOW" in
        "pypi.yml")
            run_workflow "$SPECIFIC_WORKFLOW" "pull_request" "$EVENTS_DIR/merged-pr.json" "$SPECIFIC_JOB"
            ;;
        "pytest.yml")
            run_workflow "$SPECIFIC_WORKFLOW" "push" "$EVENTS_DIR/push-main.json" "$SPECIFIC_JOB"
            ;;
        "test-pypi-sync.yml")
            run_workflow "$SPECIFIC_WORKFLOW" "workflow_dispatch" "" "$SPECIFIC_JOB"
            ;;
        *)
            run_workflow "$SPECIFIC_WORKFLOW" "push" "$EVENTS_DIR/push-main.json" "$SPECIFIC_JOB"
            ;;
    esac
    exit 0
fi

# Run all workflows
print_header "1. pytest.yml - Build and Test"
print_color "$GREEN" "Testing with push to main..."
run_workflow "pytest.yml" "push" "$EVENTS_DIR/push-main.json" "pytest"

print_color "$GREEN" "Testing with push tag..."
run_workflow "pytest.yml" "push" "$EVENTS_DIR/push-tag.json" "pytest"

print_color "$GREEN" "Testing with pull request..."
run_workflow "pytest.yml" "pull_request" "$EVENTS_DIR/pull-request.json" "pytest"

print_header "2. pypi.yml - Publish to PyPI"
print_color "$GREEN" "Testing version-bump job..."
run_workflow "pypi.yml" "pull_request" "$EVENTS_DIR/merged-pr.json" "version-bump"

print_color "$GREEN" "Testing pypi job..."
run_workflow "pypi.yml" "pull_request" "$EVENTS_DIR/merged-pr.json" "pypi"

print_header "3. test-pypi-sync.yml - Test Version Sync"
if [ -f "$WORKFLOWS_DIR/test-pypi-sync.yml" ]; then
    print_color "$GREEN" "Testing version sync workflow..."
    run_workflow "test-pypi-sync.yml" "workflow_dispatch" "" "test-version-sync"

    print_color "$GREEN" "Testing with pull request..."
    run_workflow "test-pypi-sync.yml" "pull_request" "$EVENTS_DIR/pull-request.json" "test-version-sync"
else
    print_color "$YELLOW" "test-pypi-sync.yml not found, skipping..."
fi

print_header "Summary"
print_color "$GREEN" "✅ All workflow tests completed!"
print_color "$YELLOW" "Note: Some jobs may fail due to missing secrets or dependencies in local environment"
print_color "$YELLOW" "This is expected behavior when testing locally with act"

# Cleanup event files if not in dry-run mode
if [ -z "$DRY_RUN" ]; then
    print_color "$BLUE" "Cleaning up event files..."
    rm -rf "$EVENTS_DIR"
fi

print_color "$GREEN" "Done!"
