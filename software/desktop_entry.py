# Source guide (2026-10-20; documentation comments only).
# Purpose: Launch the native desktop application from source
# Inputs: Python with Tk and desktop dependencies
# Generates/returns: Desktop workspace; results go to user-selected folder
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Entrypoint shared by the Python source and the windowed portable executable."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
import sys
import traceback


def worker(args):
    """Windowed executables may have no stdout; keep a UTF-8 worker log instead."""
    if args[:1] == ["--worker-log"]:
        if len(args) < 2:
            return 2
        path = Path(args[1])
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8", buffering=1) as stream:
            with redirect_stdout(stream), redirect_stderr(stream):
                return dispatch_worker(args[2:])
    return dispatch_worker(args)


def dispatch_worker(args):
    try:
        from zyme_designer.cli import main
        main(args)
        return 0
    except SystemExit as exc:
        return exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
    except Exception:
        traceback.print_exc()
        return 1


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args[:1] == ["--worker"]:
        return worker(args[1:])
    from zyme_designer.desktop import main as desktop_main
    return desktop_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
