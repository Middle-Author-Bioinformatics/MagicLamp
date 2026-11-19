#!/usr/bin/env python3
"""The single supported dispatcher for the consolidated latest MagicLamp suite."""
import argparse
import sys
from latest import __version__
from latest.registry import GENIES, resolve_name


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"help", "--help", "-h", "--list"}:
        print("MagicLamp {} — latest\n".format(__version__))
        print("Usage: python MagicLamp.py <Genie> -bin_dir INPUT -bin_ext EXT -out OUTPUT [--orfs | --gbk]\n")
        for name, (_, _, description) in GENIES.items():
            print("  {:16s} {}".format(name, description))
        print("\nProtein FASTA: --orfs, no clustering. GenBank: --gbk, genomic clustering.\n"
              "Validate libraries: python MagicLamp.py validate\n"
              "AbxGenie and BTEXGenie are excluded; historical implementations are in legacy/.")
        return 0
    if argv[0] in {"--version", "version"}:
        print(__version__)
        return 0
    if argv[0] == "validate":
        from latest.models import load_models
        p = argparse.ArgumentParser(prog="MagicLamp.py validate")
        p.add_argument("-genie", help="Validate a named built-in library (default all)")
        args = p.parse_args(argv[1:])
        try:
            names = [resolve_name(args.genie)] if args.genie else [n for n, (d, _, _) in GENIES.items() if d]
            for name in names:
                if GENIES[name][0] is None:
                    raise ValueError("Choose a built-in library to validate")
                models, root = load_models(name)
                print("{}: {} validated profiles ({})".format(name, len(models), root))
            return 0
        except (ValueError, OSError) as exc:
            print("MagicLamp: {}".format(exc), file=sys.stderr)
            return 1
    try:
        name = resolve_name(argv.pop(0))
    except ValueError as exc:
        print("MagicLamp: {}".format(exc), file=sys.stderr)
        return 2
    from latest.engine import main as run
    return run(name, argv)


if __name__ == "__main__":
    raise SystemExit(main())
