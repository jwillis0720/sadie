# import mkdocs.commands.build
# import mkdocs.commands.serve
# import mkdocs.config
# import mkdocs.utils
import re
from pathlib import Path

import typer

app = typer.Typer()


def take_out_termy_divs(content):
    new_content = []
    breaker = False
    for line in content.split("\n"):
        if re.match(r"<div class=\"termy\">", line):
            breaker = True
            continue
        if breaker:
            if re.match(r"</div", line):
                breaker = False
                continue
            if re.match(r"--->", line):
                continue
            if re.match(r"!!!", line):
                continue

        new_content.append(line)
    return "\n".join(new_content)


def generate_readme_content():
    index_md = Path("docs/index.md")
    content = index_md.read_text("utf-8")
    new_content = take_out_termy_divs(content)
    return new_content


@app.command()
def generate_readme():
    """
    Generate README.md content from main index.md
    """
    typer.echo("Generating README")
    readme_path = Path("README.md")
    new_content = generate_readme_content()
    readme_path.write_text(new_content, encoding="utf-8")


@app.command()
def verify_readme():
    """
    Verify README.md content from main index.md
    """
    typer.echo("Verifying README")
    readme_path = Path("README.md")
    generated_content = generate_readme_content()
    readme_content = readme_path.read_text("utf-8")
    if generated_content != readme_content:
        typer.secho("README.md outdated from the latest index.md", color=typer.colors.RED)
        raise typer.Abort()
    typer.echo("Valid README ✅")


if __name__ == "__main__":
    app()
