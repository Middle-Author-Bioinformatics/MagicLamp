#!/usr/bin/env python3
"""Compatibility entry point; installed commands use latest.cli."""
from latest.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
