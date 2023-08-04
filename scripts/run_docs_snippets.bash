set -e

echo "Running docs_src/annotation"
cd docs/docs_src/annotation
for src in $(find . -name '*.py'); do
    echo "Running $src"
    python $src
done
cd ../../../

echo "Running docs_src/reference"
cd docs/docs_src/reference
for src in $(find . -name '*.py'); do
    echo "Running $src"
    python $src
done
for src in $(find . -name '*.bash'); do
    echo "Running $src"
    bash $src
done
cd ../../../
