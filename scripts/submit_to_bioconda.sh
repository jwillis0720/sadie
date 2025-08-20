#!/bin/bash

# Script to submit sadie-antibody to Bioconda
# Usage: ./submit_to_bioconda.sh [YOUR_GITHUB_USERNAME]

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if username provided
if [ -z "$1" ]; then
    echo -e "${RED}Error: GitHub username required${NC}"
    echo "Usage: $0 <github_username>"
    exit 1
fi

GITHUB_USER=$1
RECIPE_NAME="sadie-antibody"
RECIPE_VERSION="2.0.0"

echo -e "${GREEN}Starting Bioconda submission process for ${RECIPE_NAME} v${RECIPE_VERSION}${NC}"

# Step 1: Fork reminder
echo -e "\n${YELLOW}Step 1: Fork the repository${NC}"
echo "Please fork the bioconda-recipes repository if you haven't already:"
echo "https://github.com/bioconda/bioconda-recipes/fork"
read -p "Press Enter when you've forked the repository..."

# Step 2: Clone and setup
echo -e "\n${YELLOW}Step 2: Cloning and setting up repository${NC}"
if [ ! -d "bioconda-recipes" ]; then
    git clone "https://github.com/${GITHUB_USER}/bioconda-recipes.git"
    cd bioconda-recipes
    git remote add upstream https://github.com/bioconda/bioconda-recipes.git
else
    echo "bioconda-recipes directory already exists, using existing clone"
    cd bioconda-recipes
fi

# Step 3: Update master and create branch
echo -e "\n${YELLOW}Step 3: Updating master and creating branch${NC}"
git checkout master
git pull upstream master
git push origin master
git checkout -b "add-${RECIPE_NAME}"

# Step 4: Copy recipe
echo -e "\n${YELLOW}Step 4: Adding recipe${NC}"
mkdir -p "recipes/${RECIPE_NAME}"

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SADIE_REPO_DIR="$(dirname "$SCRIPT_DIR")"

if [ -f "${SADIE_REPO_DIR}/meta.yaml" ]; then
    cp "${SADIE_REPO_DIR}/meta.yaml" "recipes/${RECIPE_NAME}/"
    echo -e "${GREEN}Recipe copied successfully${NC}"
else
    echo -e "${RED}Error: meta.yaml not found in ${SADIE_REPO_DIR}${NC}"
    exit 1
fi

# Step 5: Optional testing
echo -e "\n${YELLOW}Step 5: Testing (Optional)${NC}"
echo "To test locally, you need bioconda-utils installed."
echo "Would you like to:"
echo "  1) Lint the recipe only"
echo "  2) Build and test the recipe (requires Docker/Singularity)"
echo "  3) Skip testing"
read -p "Enter choice (1/2/3): " test_choice

case $test_choice in
    1)
        echo "Running lint check..."
        if command -v bioconda-utils &> /dev/null; then
            bioconda-utils lint "recipes/${RECIPE_NAME}"
        else
            echo -e "${YELLOW}bioconda-utils not found. Skipping lint.${NC}"
        fi
        ;;
    2)
        echo "Building and testing recipe..."
        if command -v bioconda-utils &> /dev/null; then
            bioconda-utils build "recipes/${RECIPE_NAME}"
        else
            echo -e "${YELLOW}bioconda-utils not found. Skipping build.${NC}"
        fi
        ;;
    3)
        echo "Skipping tests"
        ;;
    *)
        echo "Invalid choice, skipping tests"
        ;;
esac

# Step 6: Commit and push
echo -e "\n${YELLOW}Step 6: Committing and pushing changes${NC}"
git add "recipes/${RECIPE_NAME}"
git status
read -p "Review the changes above. Press Enter to commit and push..."
git commit -m "Add ${RECIPE_NAME} ${RECIPE_VERSION}"
git push origin "add-${RECIPE_NAME}"

# Step 7: Create PR
echo -e "\n${YELLOW}Step 7: Create Pull Request${NC}"
echo -e "${GREEN}Changes pushed successfully!${NC}"
echo ""
echo "Now go to: https://github.com/${GITHUB_USER}/bioconda-recipes"
echo "You should see a prompt to create a pull request."
echo ""
echo "PR Title: Add ${RECIPE_NAME} ${RECIPE_VERSION}"
echo "PR Description: Adding ${RECIPE_NAME} package version ${RECIPE_VERSION} to bioconda"
echo ""
echo "The CI will automatically test your recipe. Wait for all checks to pass."
echo "A bioconda maintainer will review and merge your PR."
echo ""
read -p "Press Enter after you've created the PR..."

# Step 8: Cleanup instructions
echo -e "\n${YELLOW}Step 8: Post-merge cleanup${NC}"
echo "After your PR is merged, run these commands to clean up:"
echo ""
echo "cd bioconda-recipes"
echo "git checkout master"
echo "git pull upstream master"
echo "git branch -D add-${RECIPE_NAME}"
echo "git push origin -d add-${RECIPE_NAME}"
echo ""
echo "Wait ~30 minutes after merge, then install with:"
echo "conda install -c bioconda ${RECIPE_NAME}"
echo ""
echo -e "${GREEN}Script complete!${NC}"