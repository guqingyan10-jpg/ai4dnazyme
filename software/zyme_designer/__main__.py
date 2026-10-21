# Source guide (2026-10-20; documentation comments only).
# Purpose: python -m zyme_designer entry
# Inputs: Command-line arguments
# Generates/returns: Dispatch to CLI
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
from .cli import main

main()
