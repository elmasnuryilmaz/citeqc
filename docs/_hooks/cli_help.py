"""MkDocs hook: render `citeqc --help` and each subcommand's help into the CLI page.

The page contains the marker ``<!-- citeqc-cli-help -->``; this replaces it with
the argparse output at build time, so the reference cannot drift from the code.
"""
import os

MARKER = "<!-- citeqc-cli-help -->"


def _help_sections() -> str:
    os.environ["COLUMNS"] = "88"           # stable wrapping
    os.environ["NO_COLOR"] = "1"           # no ANSI colours (argparse colours help on 3.14+)
    from citeqc.cli import build_parser

    parser = build_parser()
    sub = next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")
    blocks = ["## `citeqc`\n\n```text\n" + parser.format_help() + "```\n"]
    for name, p in sub.choices.items():
        blocks.append(f"## `citeqc {name}`\n\n```text\n{p.format_help()}```\n")
    return "\n".join(blocks)


def on_page_markdown(markdown, page, config, files):
    if MARKER in markdown:
        return markdown.replace(MARKER, _help_sections())
    return markdown
