set -e

echo "Running docs_src/object"
cd docs/docs_src/object
for src in $(find . -name '*.py'); do
    echo "Running $src"
    python $src
done
for src in $(find . -name '*.bash'); do
    echo "Running $src"
    bash $src
done
cd ../../../

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

echo "Running docs_src/renumbering"
cd docs/docs_src/renumbering
for src in $(find . -name '*.py'); do
    echo "Running $src"
    python $src
done
for src in $(find . -name '*.bash'); do
    echo "Running $src"
    bash $src
done
cd ../../../
