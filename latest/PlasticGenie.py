"""Plastic/pollutant-associated profiles. Shared implementation: latest.engine."""
if __package__ in (None, ""):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from latest.engine import main as _run
else:
    from .engine import main as _run

def main(argv=None):
    return _run('PlasticGenie', argv)

if __name__ == "__main__":
    raise SystemExit(main())
