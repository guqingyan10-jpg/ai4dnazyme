# Source guide (2026-10-20; documentation comments only).
# Purpose: Native desktop UI and worker orchestration
# Inputs: User-selected FASTA/settings/output folder
# Generates/returns: Design jobs, progress, reports and selected-design display
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Local native desktop UI. No HTTP server, cloud upload or shell commands."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
import json
import math
import os
from pathlib import Path
import queue
import re
import subprocess
import sys
import threading
import time
import uuid
import webbrowser

from .core import normalize_rna

STRATEGY_817 = "8–17 / AG sites"
STRATEGY_GENERIC = "Custom dinucleotide motif"
FOLD_PYTHON = "Calculate with ViennaRNA"
FOLD_NONE = "Sequence analysis only"
FOLD_SITE = "Import site probabilities"
FOLD_LUNP = "Import RNAplfold output"
SCAFFOLD_ARMS = "Binding arms only"
SCAFFOLD_BUILTIN = "BIRC5 / miR-21 scaffold"
SCAFFOLD_CUSTOM = "Custom scaffold JSON"
PHASE_NAMES = {
    "SEQUENCE_ENUMERATION": "Enumerating target sites",
    "RNA_ACCESSIBILITY": "Evaluating RNA accessibility",
    "TRANSCRIPTOME_SEQUENCE_SCAN": "Screening the reference transcriptome",
    "OPTIONAL_3D_DIAGNOSTIC": "Analysing imported coordinates",
    "ALPHAFOLD_INPUT_EXPORT": "Preparing AlphaFold 3 inputs",
    "COMPLETE": "Analysis complete",
}

# A restrained native theme shared with the software's Wiki-facing results.
# These tokens affect presentation only; no calculation or validation uses them.
UI_COLORS = {
    "paper": "#FFF9F1",
    "card": "#FFFFFF",
    "ink": "#163C52",
    "muted": "#627681",
    "teal": "#147D78",
    "teal_dark": "#0C625F",
    "teal_wash": "#EAF5F1",
    "peach": "#F1C9B3",
    "field": "#FAFCFB",
    "line": "#DBE5E2",
}


def set_windows_dpi_awareness() -> str:
    """Request native pixel rendering before Tk creates any Windows HWND."""
    if sys.platform != "win32":
        return "platform-default"
    import ctypes
    try:
        if ctypes.windll.user32.SetProcessDpiAwarenessContext(ctypes.c_void_p(-4)):
            return "per-monitor-v2"
    except (AttributeError, OSError):
        pass
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(2)
        return "per-monitor-or-manifest"
    except (AttributeError, OSError):
        try:
            ctypes.windll.user32.SetProcessDPIAware()
        except (AttributeError, OSError):
            pass
        return "manifest-or-system"


def resource_root() -> Path:
    """PyInstaller data tree or the source-project root."""
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS)
    return Path(__file__).resolve().parents[1]


def default_results_parent() -> Path:
    return Path.home() / "Documents" / "ZymeSiteDesigner" / "results"


@dataclass
class DesktopConfig:
    input_mode: str = "file"
    target_fasta: str = ""
    pasted_text: str = ""
    gene: str = ""
    transcriptome: str = ""
    reference_label: str = "No reference supplied"
    strategy: str = STRATEGY_817
    motif: str = "AG"
    cut_positions: str = ""
    left_lengths: str = "7,8,9"
    right_lengths: str = "7,8,9"
    preferred_min: int = 7
    preferred_max: int = 9
    temperature: float = 37.0
    window: int = 200
    span: int = 150
    rna_salt: float = 1.021
    fold_mode: str = FOLD_PYTHON
    accessibility_path: str = ""
    accessibility_manifest: str = ""
    screen_top: int = 20
    mismatches: int = 2
    scaffold_mode: str = SCAFFOLD_ARMS
    scaffold_path: str = ""
    output_parent: str = ""
    demonstration: bool = False


def parse_lengths(value: str) -> tuple[int, ...]:
    try:
        values = tuple(sorted(set(int(x.strip()) for x in value.split(",") if x.strip())))
    except ValueError as exc:
        raise ValueError("Enter comma-separated integer arm lengths, for example 7,8,9.") from exc
    if not values or any(x < 4 or x > 20 for x in values):
        raise ValueError("Arm lengths must be 4–20 nt.")
    return values


def parse_cut_positions(value: str) -> tuple[int, ...]:
    if not value.strip():
        return ()
    try:
        positions = tuple(sorted(set(int(x.strip()) for x in value.split(","))))
    except ValueError as exc:
        raise ValueError("Enter comma-separated cut positions, for example 408,420.") from exc
    if any(x < 1 for x in positions):
        raise ValueError("Cut positions must be positive, 1-based coordinates.")
    return positions


def pasted_fasta_text(text: str, gene: str = "") -> str:
    """Accept one FASTA record or raw unambiguous RNA/DNA; never invent a gene."""
    stripped = text.strip()
    if not stripped:
        raise ValueError("Enter a target transcript sequence.")
    lines = [line.strip() for line in stripped.splitlines() if line.strip()]
    headers = [line for line in lines if line.startswith(">")]
    if headers:
        if len(headers) != 1 or not lines[0].startswith(">"):
            raise ValueError("Provide one transcript, with a single FASTA header on the first line.")
        header = lines[0]
        if len(header) < 2:
            raise ValueError("The FASTA header must not be empty.")
        sequence = normalize_rna("".join(lines[1:]))
    else:
        header = ">USER_INPUT" + (f" gene={gene}" if gene else "")
        sequence = normalize_rna("".join(lines))
    return header + "\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"


def validate_config(config: DesktopConfig) -> DesktopConfig:
    if config.input_mode not in {"file", "paste"}:
        raise ValueError("Choose a FASTA file or paste a sequence.")
    if config.input_mode == "file" and not Path(config.target_fasta).is_file():
        raise ValueError("Select an existing target FASTA file.")
    if config.input_mode == "paste":
        pasted_fasta_text(config.pasted_text, config.gene)
    if config.gene and not re.fullmatch(r"[A-Za-z0-9_.-]+", config.gene):
        raise ValueError("Gene labels may contain letters, numbers, underscores, dots and hyphens.")
    left, right = parse_lengths(config.left_lengths), parse_lengths(config.right_lengths)
    motif = config.motif.strip().upper().replace("T", "U")
    if len(motif) != 2 or set(motif) - set("ACGURYSWKMBDHVN"):
        raise ValueError("The motif must contain two IUPAC RNA characters, for example AG, RY or NN.")
    if config.strategy == STRATEGY_817 and motif != "AG":
        raise ValueError("The 8–17 preset uses AG sites. Choose a custom strategy for other motifs.")
    if config.scaffold_mode == SCAFFOLD_BUILTIN and motif != "AG":
        raise ValueError("The BIRC5 / miR-21 scaffold requires the AG motif.")
    parse_cut_positions(config.cut_positions)
    if not 6 <= config.preferred_min <= config.preferred_max <= 12:
        raise ValueError("Preferred arm lengths must satisfy 6 <= minimum <= maximum <= 12.")
    if not math.isfinite(config.temperature) or not 0 <= config.temperature <= 100:
        raise ValueError("Temperature must be a finite value between 0 and 100 °C.")
    if config.window < 1 or not 1 <= config.span <= config.window:
        raise ValueError("RNA folding settings must satisfy 1 <= span <= window.")
    if config.fold_mode == FOLD_PYTHON and config.window < max(left) + max(right) + 2:
        raise ValueError("The folding window must cover the longest candidate target sequence.")
    if not math.isfinite(config.rna_salt) or config.rna_salt <= 0:
        raise ValueError("The RNA folding model requires a positive monovalent salt concentration in M.")
    if config.screen_top < 0 or config.mismatches not in {0, 1, 2}:
        raise ValueError("Screening count must be non-negative; mismatch allowance must be 0, 1 or 2.")
    if config.transcriptome and not Path(config.transcriptome).is_file():
        raise ValueError("The selected reference transcriptome file does not exist.")
    if config.transcriptome and not config.reference_label.strip():
        raise ValueError("Provide a reference label including its organism, version and coverage.")
    if config.fold_mode in {FOLD_SITE, FOLD_LUNP}:
        if not Path(config.accessibility_path).is_file() or not Path(config.accessibility_manifest).is_file():
            raise ValueError("Select both an accessibility file and its matching manifest JSON.")
    if config.scaffold_mode == SCAFFOLD_CUSTOM and not Path(config.scaffold_path).is_file():
        raise ValueError("Select an existing custom scaffold JSON.")
    if not config.output_parent.strip():
        raise ValueError("Choose a results folder.")
    return config


def unique_output_path(parent: Path, prefix: str = "design") -> Path:
    """Only return a fresh path; CLI owns creation of the leaf directory."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return parent / f"{prefix}_{stamp}_{uuid.uuid4().hex[:8]}"


def build_cli_args(config: DesktopConfig, target: Path, output: Path) -> list[str]:
    config = validate_config(config)
    args = ["design", "--target-fasta", str(target), "--out", str(output),
            "--left-lengths", ",".join(map(str, parse_lengths(config.left_lengths))),
            "--right-lengths", ",".join(map(str, parse_lengths(config.right_lengths))),
            "--motif", config.motif.strip().upper().replace("T", "U"),
            "--preferred-arm-min", str(config.preferred_min), "--preferred-arm-max", str(config.preferred_max),
            "--temperature", str(config.temperature), "--window", str(config.window), "--span", str(config.span),
            "--screen-top", str(config.screen_top), "--mismatches", str(config.mismatches)]
    if config.gene:
        args += ["--gene", config.gene]
    for cut in parse_cut_positions(config.cut_positions):
        args += ["--cut-after", str(cut)]
    if config.transcriptome:
        args += ["--transcriptome", config.transcriptome, "--reference-label", config.reference_label]
    if config.fold_mode == FOLD_PYTHON:
        args += ["--fold", "--fold-backend", "python", "--rna-salt", str(config.rna_salt)]
    elif config.fold_mode in {FOLD_SITE, FOLD_LUNP}:
        args += ["--site-table" if config.fold_mode == FOLD_SITE else "--lunp", config.accessibility_path,
                 "--accessibility-manifest", config.accessibility_manifest]
    if config.scaffold_mode == SCAFFOLD_BUILTIN:
        args += ["--scaffold", "builtin-birc5-mir21"]
    elif config.scaffold_mode == SCAFFOLD_CUSTOM:
        args += ["--scaffold", config.scaffold_path]
    if config.demonstration:
        for position in (401, 416, 419):
            args += ["--include-start", str(position)]
    return args


def worker_command(args: list[str], log_path: Path) -> list[str]:
    prefix = [sys.executable] if getattr(sys, "frozen", False) else [sys.executable, str(resource_root() / "desktop_entry.py")]
    return prefix + ["--worker", "--worker-log", str(log_path)] + args


def demo_config(parent: Path) -> DesktopConfig:
    examples = resource_root() / "examples"
    return DesktopConfig(target_fasta=str(examples / "BIRC5_NM_001168_3.fasta"), gene="BIRC5",
                         transcriptome=str(examples / "tiny_reference_NOT_HUMAN_TRANSCRIPTOME.fasta"),
                         reference_label="TINY DEMONSTRATION: BIRC5 + short hit windows; NOT full human transcriptome",
                         left_lengths="7", right_lengths="7", fold_mode=FOLD_SITE,
                         accessibility_path=str(examples / "BIRC5_accessibility_W200_selected.csv"),
                         accessibility_manifest=str(examples / "BIRC5_accessibility_W200_selected.manifest.json"),
                         scaffold_mode=SCAFFOLD_BUILTIN, output_parent=str(parent), demonstration=True)


class DesignerApp:
    def __init__(self, root, results_parent=None, open_browser=True, on_finished=None):
        import tkinter as tk
        from tkinter import ttk
        self.tk, self.ttk, self.root = tk, ttk, root
        self.open_browser, self.on_finished = open_browser, on_finished
        self.process = None
        self.output = None
        self.worker_log = None
        self.log_offset = 0
        self.events = queue.Queue()
        self.current_token = 0
        self.cancelled = False
        self.closing = False
        self.demo_mode = False
        self.result_data = None
        self.candidate_lookup = {}
        from . import __version__
        self.root.title("Zyme Site Designer")
        # Tk fonts use points, while geometry/padding use physical pixels.
        # Scale the preferred window as well as fonts on high-DPI displays.
        screen_scale = max(1.0, float(root.tk.call("tk", "scaling")) / (4 / 3))
        width = min(round(1120 * screen_scale), max(840, root.winfo_screenwidth() - round(100 * screen_scale)))
        height = min(round(800 * screen_scale), max(560, root.winfo_screenheight() - round(110 * screen_scale)))
        self.root.geometry(f"{width}x{height}")
        self.root.minsize(840, 560)
        self.root.configure(bg=UI_COLORS["paper"])
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TFrame", background=UI_COLORS["paper"])
        style.configure("Card.TFrame", background=UI_COLORS["card"], relief="solid", borderwidth=1,
                        bordercolor=UI_COLORS["line"])
        style.configure("TLabel", background=UI_COLORS["paper"], foreground=UI_COLORS["ink"], font=("Segoe UI", 10))
        style.configure("Muted.TLabel", foreground=UI_COLORS["muted"])
        style.configure("Card.TLabel", background=UI_COLORS["card"])
        style.configure("CardMuted.TLabel", background=UI_COLORS["card"], foreground=UI_COLORS["muted"])
        style.configure("Section.TLabel", font=("Segoe UI", 12, "bold"), foreground=UI_COLORS["ink"])
        style.configure("CardSection.TLabel", font=("Segoe UI", 12, "bold"), background=UI_COLORS["card"])
        style.configure("TButton", padding=(13, 8), font=("Segoe UI", 10), background=UI_COLORS["teal_wash"],
                        foreground=UI_COLORS["ink"], borderwidth=0, focusthickness=2, focuscolor=UI_COLORS["teal"])
        style.map("TButton", background=[("active", "#D9ECE5")], foreground=[("disabled", "#8C9D9F")])
        style.configure("Accent.TButton", background=UI_COLORS["teal"], foreground="white", font=("Segoe UI", 10, "bold"), padding=(22, 11))
        style.map("Accent.TButton", background=[("disabled", "#A3BCB7"), ("active", UI_COLORS["teal_dark"])],
                  foreground=[("disabled", "#F5F8F6"), ("active", "white")])
        style.configure("TEntry", padding=8, font=("Segoe UI", 10), fieldbackground=UI_COLORS["field"],
                        foreground=UI_COLORS["ink"], bordercolor=UI_COLORS["line"], lightcolor=UI_COLORS["line"],
                        darkcolor=UI_COLORS["line"], insertcolor=UI_COLORS["ink"])
        style.map("TEntry", bordercolor=[("focus", UI_COLORS["teal"])])
        style.configure("TCombobox", padding=8, font=("Segoe UI", 10), foreground=UI_COLORS["ink"],
                        fieldbackground=UI_COLORS["field"], background=UI_COLORS["teal_wash"],
                        bordercolor=UI_COLORS["line"], arrowcolor=UI_COLORS["teal"])
        style.map("TCombobox", fieldbackground=[("readonly", UI_COLORS["field"])],
                  foreground=[("readonly", UI_COLORS["ink"])], bordercolor=[("focus", UI_COLORS["teal"])])
        style.configure("TRadiobutton", background=UI_COLORS["card"], foreground=UI_COLORS["ink"], font=("Segoe UI", 10),
                        indicatorcolor=UI_COLORS["card"])
        style.map("TRadiobutton", background=[("active", UI_COLORS["card"])],
                  indicatorcolor=[("selected", UI_COLORS["teal"])])
        style.configure("TNotebook", background=UI_COLORS["paper"], borderwidth=0, tabmargins=(0, 0, 0, 0))
        style.configure("TNotebook.Tab", padding=(20, 11), font=("Segoe UI", 10, "bold"),
                        background="#EFECE3", foreground=UI_COLORS["muted"], borderwidth=0)
        style.map("TNotebook.Tab", background=[("selected", UI_COLORS["teal_wash"]), ("active", "#E6EFE9")],
                  foreground=[("selected", UI_COLORS["teal"])])
        style.configure("Horizontal.TProgressbar", background=UI_COLORS["teal"], troughcolor=UI_COLORS["teal_wash"],
                        borderwidth=0, lightcolor=UI_COLORS["teal"], darkcolor=UI_COLORS["teal"])
        rowheight = max(28, int(22 * float(root.tk.call("tk", "scaling"))))
        style.configure("Treeview", rowheight=rowheight, font=("Segoe UI", 9), background=UI_COLORS["card"],
                        fieldbackground=UI_COLORS["card"], foreground=UI_COLORS["ink"], borderwidth=0)
        style.map("Treeview", background=[("selected", UI_COLORS["teal"])], foreground=[("selected", "white")])
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"), background=UI_COLORS["teal_wash"],
                        foreground=UI_COLORS["ink"], padding=(10, 9), relief="flat")
        style.configure("TMenubutton", padding=(12, 8), font=("Segoe UI", 10), background=UI_COLORS["teal_wash"],
                        foreground=UI_COLORS["ink"], borderwidth=0, arrowcolor=UI_COLORS["teal"])
        header = tk.Frame(root, bg=UI_COLORS["paper"], padx=28, pady=18)
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(0, weight=1)
        tk.Label(header, text="RNA TARGETS / DNA BINDING ARMS", font=("Segoe UI", 8, "bold"),
                 fg=UI_COLORS["teal"], bg=UI_COLORS["paper"]).grid(row=0, column=0, sticky="w")
        tk.Label(header, text="Zyme Site Designer", font=("Segoe UI", 22, "bold"),
                 fg=UI_COLORS["ink"], bg=UI_COLORS["paper"]).grid(row=1, column=0, sticky="w", pady=(2, 3))
        tk.Label(header, text="Find accessible RNA sites. Design complementary DNA arms.", font=("Segoe UI", 10),
                 fg=UI_COLORS["muted"], bg=UI_COLORS["paper"]).grid(row=2, column=0, sticky="w")
        tk.Label(header, text="LOCAL\nWORKSPACE", font=("Segoe UI", 8, "bold"), justify="left",
                 fg=UI_COLORS["teal"], bg=UI_COLORS["teal_wash"], padx=14, pady=9).grid(row=1, column=1, rowspan=2, padx=(18, 0))
        menu = tk.Menu(root)
        file_menu = tk.Menu(menu, tearoff=False)
        file_menu.add_command(label="Load BIRC5 example", command=self.load_example)
        file_menu.add_command(label="Open results folder", command=self.open_output)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.close)
        menu.add_cascade(label="File", menu=file_menu)
        help_menu = tk.Menu(menu, tearoff=False)
        help_menu.add_command(label="Quick start", command=self.show_help)
        help_menu.add_command(label=f"About · version {__version__}", command=self.show_about)
        menu.add_cascade(label="Help", menu=help_menu)
        root.configure(menu=menu)
        self.vars = {}
        config = DesktopConfig(output_parent=str(results_parent or default_results_parent()))
        for name, value in vars(config).items():
            if name not in {"pasted_text", "demonstration"}:
                self.vars[name] = tk.StringVar(value=str(value))
                self.vars[name].trace_add("write", self.clear_demo_marker)
        self.status = tk.StringVar(value="Ready")
        self.notebook = ttk.Notebook(root)
        self.notebook.grid(row=1, column=0, sticky="nsew", padx=22, pady=(0, 0))
        self.target_tab, self.param_tab, self.log_tab, self.result_tab = [ttk.Frame(self.notebook) for _ in range(4)]
        for tab, title in zip((self.target_tab, self.param_tab, self.result_tab, self.log_tab), ("01  Transcript", "02  Design", "03  Results", "Activity log")):
            self.notebook.add(tab, text=title)
        self.build_target_tab()
        self.build_param_tab()
        self.build_log_tab()
        self.build_result_tab()
        self.footer = ttk.Frame(root, padding=(18, 12), style="Card.TFrame")
        self.footer.grid(row=2, column=0, sticky="ew", padx=22, pady=(10, 18))
        self.footer.columnconfigure(0, weight=1)
        self.footer_status = ttk.Label(self.footer, textvariable=self.status, style="CardMuted.TLabel", wraplength=330)
        self.footer_status.grid(row=0, column=0, sticky="w")
        self.cancel_button = ttk.Button(self.footer, text="Cancel", command=self.cancel, state="disabled")
        self.cancel_button.grid(row=0, column=1, padx=8)
        self.cancel_button.grid_remove()
        self.report_button = ttk.Button(self.footer, text="Open report", command=self.open_result, state="disabled")
        self.report_button.grid(row=0, column=2, padx=(0, 10))
        self.start_button = ttk.Button(self.footer, text="Run analysis", style="Accent.TButton", command=self.start)
        self.start_button.grid(row=0, column=3)
        # Kept as an API handle for the smoke-demo worker tests, not another
        # default action in the footer. Loading an example is a File-menu task.
        self.demo_button = ttk.Button(self.footer, text="Example", command=self.run_demo)
        self.footer.bind("<Configure>", lambda e: self.footer_status.configure(wraplength=max(120, e.width - 380)))
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.after(300, self.poll)

    def scroll_page(self, parent):
        """Scrollable content stays inside the expandable middle grid row."""
        tk, ttk = self.tk, self.ttk
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        canvas = tk.Canvas(parent, bg=UI_COLORS["paper"], highlightthickness=0, width=1, height=1)
        bar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=bar.set)
        canvas.grid(row=0, column=0, sticky="nsew")
        bar.grid(row=0, column=1, sticky="ns")
        body = ttk.Frame(canvas, padding=(4, 20, 4, 12))
        item = canvas.create_window((0, 0), window=body, anchor="nw")
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(item, width=event.width))
        body.bind("<Configure>", lambda event: canvas.configure(scrollregion=canvas.bbox("all")))
        def wheel(event):
            widget = event.widget
            # Scope the wheel to this page. Never steal a Text or table scroll.
            if widget.winfo_class() in {"Text", "Treeview", "TCombobox"}:
                return
            while widget is not None:
                if widget == parent:
                    canvas.yview_scroll(-1 if event.delta > 0 else 1, "units")
                    break
                widget = getattr(widget, "master", None)
        parent.winfo_toplevel().bind("<MouseWheel>", wheel, add="+")
        return body

    def section(self, parent, title, subtitle=None):
        frame = self.ttk.Frame(parent, padding=18, style="Card.TFrame")
        frame.pack(fill="x", pady=(0, 18))
        self.ttk.Label(frame, text=title, style="CardSection.TLabel").pack(anchor="w")
        if subtitle:
            self.ttk.Label(frame, text=subtitle, style="CardMuted.TLabel", wraplength=640).pack(anchor="w", pady=(5, 10))
        fields = self.ttk.Frame(frame, style="Card.TFrame", padding=0, borderwidth=0)
        fields.pack(fill="x", pady=(8, 0))
        return fields

    def entry_row(self, parent, row, label, key, browse=False, combo=None):
        ttk = self.ttk
        ttk.Label(parent, text=label, style="Card.TLabel").grid(row=row, column=0, sticky="w", padx=(0, 16), pady=6)
        if combo:
            widget = ttk.Combobox(parent, textvariable=self.vars[key], values=combo, state="readonly")
        else:
            widget = ttk.Entry(parent, textvariable=self.vars[key])
        widget.grid(row=row, column=1, sticky="ew", pady=6)
        parent.columnconfigure(1, weight=1)
        if browse:
            ttk.Button(parent, text="Browse", command=lambda: self.browse(key, directory=key == "output_parent")).grid(row=row, column=2, padx=(8, 0))
        return widget

    def build_target_tab(self):
        ttk, tk = self.ttk, self.tk
        frame = self.scroll_page(self.target_tab)
        target = ttk.Frame(frame, style="Card.TFrame", padding=18)
        target.pack(fill="x", pady=(0, 18))
        ttk.Label(target, text="Your target transcript", style="CardSection.TLabel").pack(anchor="w", pady=(0, 6))
        ttk.Label(target, text="Choose one complete transcript, or paste its sequence.", style="CardMuted.TLabel").pack(anchor="w", pady=(0, 16))
        choice = ttk.Frame(target, style="Card.TFrame", borderwidth=0)
        choice.pack(fill="x")
        self.input_choice = choice
        ttk.Radiobutton(choice, text="FASTA file", variable=self.vars["input_mode"], value="file", command=self.input_mode_changed).pack(side="left")
        ttk.Radiobutton(choice, text="Paste sequence", variable=self.vars["input_mode"], value="paste", command=self.input_mode_changed).pack(side="left", padx=24)
        self.file_fields = ttk.Frame(target, style="Card.TFrame", borderwidth=0)
        self.file_fields.pack(fill="x", pady=12)
        self.entry_row(self.file_fields, 0, "Target FASTA", "target_fasta", browse=True)
        self.paste_fields = ttk.Frame(target, style="Card.TFrame", borderwidth=0)
        self.paste_fields.pack(fill="x", pady=12)
        self.paste = tk.Text(self.paste_fields, height=6, wrap="word", font=("Consolas", 10), relief="flat", borderwidth=0,
                             bg=UI_COLORS["field"], fg=UI_COLORS["ink"], insertbackground=UI_COLORS["ink"],
                             highlightthickness=1, highlightbackground=UI_COLORS["line"], highlightcolor=UI_COLORS["teal"], padx=12, pady=10)
        self.paste.pack(fill="x")
        self.paste_fields.pack_forget()
        identity = ttk.Frame(target, style="Card.TFrame", borderwidth=0)
        identity.pack(fill="x")
        self.entry_row(identity, 0, "Gene label (optional)", "gene")
        output = self.section(frame, "Results destination", "Every run gets a new folder. Existing results stay unchanged.")
        self.entry_row(output, 0, "Results folder", "output_parent", browse=True)
        reference = self.section(frame, "Reference transcriptome", "Optional. Supply a transcriptome FASTA to screen for similar binding sites.")
        self.entry_row(reference, 0, "Reference FASTA", "transcriptome", browse=True)
        self.entry_row(reference, 1, "Reference label", "reference_label")
        ttk.Label(reference, text="No reference? Off-target screening remains not assessed.", style="CardMuted.TLabel", wraplength=640).grid(row=2, column=0, columnspan=3, sticky="w", pady=(8, 0))

    def input_mode_changed(self):
        if self.vars["input_mode"].get() == "paste":
            self.file_fields.pack_forget()
            self.paste_fields.pack(fill="x", pady=12, after=self.input_choice)
        else:
            self.paste_fields.pack_forget()
            self.file_fields.pack(fill="x", pady=12, after=self.input_choice)

    def build_param_tab(self):
        ttk = self.ttk
        frame = self.scroll_page(self.param_tab)
        basic = self.section(frame, "Build your design", "Start with AG sites and 7–9 nt arms, or choose your own search.")
        self.entry_row(basic, 0, "Site strategy", "strategy", combo=(STRATEGY_817, STRATEGY_GENERIC)).bind("<<ComboboxSelected>>", self.strategy_changed)
        self.entry_row(basic, 1, "Motif", "motif")
        self.entry_row(basic, 2, "Left lengths (nt)", "left_lengths")
        self.entry_row(basic, 3, "Right lengths (nt)", "right_lengths")
        self.entry_row(basic, 4, "Temperature (°C)", "temperature")
        self.entry_row(basic, 5, "RNA accessibility", "fold_mode", combo=(FOLD_PYTHON, FOLD_NONE, FOLD_SITE, FOLD_LUNP)).bind("<<ComboboxSelected>>", self.conditional_fields)
        self.entry_row(basic, 6, "Output sequences", "scaffold_mode", combo=(SCAFFOLD_ARMS, SCAFFOLD_BUILTIN, SCAFFOLD_CUSTOM)).bind("<<ComboboxSelected>>", self.conditional_fields)
        self.import_fields = ttk.Frame(frame, style="Card.TFrame", padding=18)
        self.entry_row(self.import_fields, 0, "Accessibility file", "accessibility_path", browse=True)
        self.entry_row(self.import_fields, 1, "Matching manifest", "accessibility_manifest", browse=True)
        self.scaffold_fields = ttk.Frame(frame, style="Card.TFrame", padding=18)
        self.entry_row(self.scaffold_fields, 0, "Scaffold JSON", "scaffold_path", browse=True)
        self.advanced_button = ttk.Button(frame, text="Advanced settings ▸", command=self.toggle_advanced)
        self.advanced_button.pack(anchor="w", pady=(4, 10))
        self.advanced_frame = ttk.Frame(frame, style="Card.TFrame", padding=18)
        self.entry_row(self.advanced_frame, 0, "Cut positions (optional)", "cut_positions")
        self.entry_row(self.advanced_frame, 1, "Preferred minimum (nt)", "preferred_min")
        self.entry_row(self.advanced_frame, 2, "Preferred maximum (nt)", "preferred_max")
        self.entry_row(self.advanced_frame, 3, "Screen top sites (0 = all)", "screen_top")
        self.entry_row(self.advanced_frame, 4, "Total arm mismatches", "mismatches", combo=("0", "1", "2"))
        self.entry_row(self.advanced_frame, 5, "Folding window W (nt)", "window")
        self.entry_row(self.advanced_frame, 6, "Pairing span L (nt)", "span")
        self.entry_row(self.advanced_frame, 7, "RNA model salt (M)", "rna_salt")
        self.advanced_open = False

    def toggle_advanced(self):
        self.advanced_open = not self.advanced_open
        if self.advanced_open:
            self.advanced_frame.pack(fill="x", pady=(0, 12))
        else:
            self.advanced_frame.pack_forget()
        self.advanced_button.configure(text="Advanced settings ▾" if self.advanced_open else "Advanced settings ▸")

    def conditional_fields(self, event=None):
        self.import_fields.pack_forget()
        self.scaffold_fields.pack_forget()
        if self.vars["fold_mode"].get() in {FOLD_SITE, FOLD_LUNP}:
            self.import_fields.pack(fill="x", pady=(0, 16), before=self.advanced_button)
        if self.vars["scaffold_mode"].get() == SCAFFOLD_CUSTOM:
            self.scaffold_fields.pack(fill="x", pady=(0, 16), before=self.advanced_button)

    def build_log_tab(self):
        ttk, tk = self.ttk, self.tk
        frame = ttk.Frame(self.log_tab, padding=20, style="Card.TFrame")
        frame.pack(fill="both", expand=True, pady=(20, 0), padx=4)
        ttk.Label(frame, text="Analysis activity", style="CardSection.TLabel").pack(anchor="w", pady=(0, 12))
        self.progress = ttk.Progressbar(frame, mode="indeterminate")
        self.progress.pack(fill="x", pady=(0, 10))
        self.progress_detail = tk.StringVar(value="Run an analysis to view its progress and log.")
        ttk.Label(frame, textvariable=self.progress_detail, style="CardMuted.TLabel", wraplength=640).pack(anchor="w", pady=(0, 12))
        from tkinter.scrolledtext import ScrolledText
        self.log_text = ScrolledText(frame, height=6, wrap="word", bg=UI_COLORS["field"], fg=UI_COLORS["ink"], font=("Consolas", 9), relief="flat", padx=12, pady=12, state="disabled")
        self.log_text.pack(fill="both", expand=True)

    def build_result_tab(self):
        ttk, tk = self.ttk, self.tk
        frame = self.scroll_page(self.result_tab)
        controls = ttk.Frame(frame)
        controls.pack(fill="x")
        ttk.Label(controls, text="Explore candidate sites", style="Section.TLabel").pack(side="left")
        exports = tk.Menu(controls, tearoff=False)
        exports.add_command(label="Results folder", command=self.open_output)
        exports.add_command(label="AlphaFold 3 input files", command=self.open_af3)
        exports.add_command(label="Structure viewer", command=self.open_structure)
        export_button = ttk.Menubutton(controls, text="Exports", menu=exports)
        export_button.pack(side="right")
        self.result_scope = tk.StringVar(value="Run an analysis to inspect and compare candidate sites.")
        ttk.Label(frame, textvariable=self.result_scope, wraplength=640, style="Muted.TLabel").pack(anchor="w", pady=8)
        self.canvas = tk.Canvas(frame, height=155, bg=UI_COLORS["card"], highlightthickness=1,
                                highlightbackground=UI_COLORS["line"])
        self.canvas.pack(fill="x", pady=(0, 10))
        self.canvas.bind("<Configure>", lambda event: self.draw_profile())
        self.canvas.bind("<Button-1>", self.select_nearest)
        columns = ("cut", "arms", "p", "exact", "score", "stage")
        table = ttk.Frame(frame)
        table.pack(fill="x")
        self.tree = ttk.Treeview(table, columns=columns, show="headings", height=7, selectmode="browse")
        for column, title, width in zip(columns, ("Cut site", "Arms L / R", "Accessibility", "Other-gene hits", "Score", "Evidence"), (90, 90, 100, 130, 80, 190)):
            self.tree.heading(column, text=title)
            self.tree.column(column, width=width, minwidth=70, stretch=True)
        self.tree.grid(row=0, column=0, sticky="nsew")
        table.columnconfigure(0, weight=1)
        ttk.Scrollbar(table, orient="vertical", command=self.tree.yview).grid(row=0, column=1, sticky="ns")
        xbar = ttk.Scrollbar(table, orient="horizontal", command=self.tree.xview)
        xbar.grid(row=1, column=0, sticky="ew")
        self.tree.configure(xscrollcommand=xbar.set)
        ybar = table.grid_slaves(row=0, column=1)[0]
        self.tree.configure(yscrollcommand=ybar.set)
        self.tree.bind("<<TreeviewSelect>>", self.select_candidate)
        self.selected_detail = tk.StringVar(value="Select a site to view its target sequence and binding arms.")
        detail = ttk.Frame(frame, style="Card.TFrame", padding=18)
        detail.pack(fill="x", pady=(18, 0))
        ttk.Label(detail, text="Selected design", style="CardSection.TLabel").pack(anchor="w", pady=(0, 10))
        ttk.Label(detail, textvariable=self.selected_detail, justify="left", wraplength=640,
                  style="Card.TLabel", font=("Consolas", 10)).pack(anchor="w")

    def browse(self, key, directory=False):
        from tkinter import filedialog
        value = filedialog.askdirectory(parent=self.root) if directory else filedialog.askopenfilename(parent=self.root, filetypes=[("All supported files", "*.*")])
        if value:
            self.vars[key].set(value)
            self.demo_mode = False

    def strategy_changed(self, event=None):
        if self.vars["strategy"].get() == STRATEGY_817:
            self.vars["motif"].set("AG")
        else:
            self.vars["scaffold_mode"].set(SCAFFOLD_ARMS)
        self.demo_mode = False

    def clear_demo_marker(self, *args):
        self.demo_mode = False

    def read_config(self):
        values = {name: variable.get().strip() for name, variable in self.vars.items()}
        for field in ("preferred_min", "preferred_max", "window", "span", "screen_top", "mismatches"):
            try:
                values[field] = int(values[field])
            except ValueError as exc:
                raise ValueError(f"{field.replace('_', ' ').capitalize()} must be an integer.") from exc
        for field in ("temperature", "rna_salt"):
            try:
                values[field] = float(values[field])
            except ValueError as exc:
                raise ValueError(f"{field.replace('_', ' ').capitalize()} must be a number.") from exc
        values["pasted_text"] = self.paste.get("1.0", "end")
        values["demonstration"] = self.demo_mode
        return validate_config(DesktopConfig(**values))

    def set_config(self, config):
        for name, value in vars(config).items():
            if name in self.vars:
                self.vars[name].set(str(value))
        self.paste.delete("1.0", "end")
        self.paste.insert("1.0", config.pasted_text)
        self.demo_mode = config.demonstration
        self.input_mode_changed()
        self.conditional_fields()

    def load_example(self):
        if self.process is not None:
            return
        parent = Path(self.vars["output_parent"].get()).expanduser()
        self.set_config(demo_config(parent))
        self.notebook.select(self.target_tab)
        self.status.set("BIRC5 example loaded")

    def show_help(self):
        from tkinter import messagebox
        messagebox.showinfo("Quick start", "1. Select a FASTA file or paste one transcript.\n"
                            "2. Review Design settings, or keep the defaults.\n"
                            "3. Choose Run analysis, then inspect Results.\n\n"
                            "Open report provides the full candidate table and exports.\n"
                            "File > Load BIRC5 example loads a small reference example.", parent=self.root)

    def show_about(self):
        from tkinter import messagebox
        from . import __version__
        messagebox.showinfo("About Zyme Site Designer", f"Zyme Site Designer {__version__}\n\n"
                            "RNA target selection and DNA binding-arm design.\n"
                            "FASTA input · CSV / JSON output · ViennaRNA accessibility\n\n"
                            "Candidate scores prioritize sequence and structure features;\n"
                            "they are not cleavage-rate predictions.\n"
                            "RNA folding salt and duplex reference energies use separate models.\n"
                            "AlphaFold exports are inputs, not predictions.\n\n"
                            "Source: MIT license. Runtime dependencies retain their own licenses.\n"
                            "See docs/METHODS.md and docs/LIMITATIONS.md.", parent=self.root)

    def append_log(self, message):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", message)
        if int(self.log_text.index("end-1c").split(".")[0]) > 10000:
            self.log_text.delete("1.0", "3000.0")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def run_demo(self):
        if self.process is not None:
            return
        parent = Path(self.vars["output_parent"].get()).expanduser()
        self.set_config(demo_config(parent))
        self.start()

    def start(self):
        from tkinter import messagebox
        if self.process is not None:
            return
        try:
            config = self.read_config()
            parent = Path(config.output_parent).expanduser().resolve()
            parent.mkdir(parents=True, exist_ok=True)
            output = unique_output_path(parent, "BIRC5_demo" if config.demonstration else "design")
            session = parent / "desktop_logs" / output.name
            session.mkdir(parents=True, exist_ok=False)
            if config.input_mode == "paste":
                target = session / "input_transcript.fasta"
                target.write_text(pasted_fasta_text(config.pasted_text, config.gene), encoding="utf-8")
            else:
                target = Path(config.target_fasta).resolve()
            log_path = session / "worker.log"
            args = build_cli_args(config, target, output)
            command = worker_command(args, log_path)
            (session / "desktop_request.json").write_text(json.dumps({"config": vars(config), "argv": command, "output": str(output), "scope": "local worker; no cloud submission"}, ensure_ascii=False, indent=2), encoding="utf-8")
            environment = os.environ.copy()
            environment.update({"PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"})
            flags = (subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP) if os.name == "nt" else 0
            process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       text=True, encoding="utf-8", errors="replace", shell=False,
                                       env=environment, creationflags=flags, cwd=str(resource_root()))
        except (ValueError, OSError) as exc:
            self.status.set("Check input")
            if self.on_finished:
                self.on_finished(False, {"error": str(exc)})
            else:
                messagebox.showerror("Check input", str(exc), parent=self.root)
            return
        self.process, self.output, self.worker_log = process, output, log_path
        self.log_offset = 0
        self.current_token += 1
        self.cancelled = False
        self.progress.start(12)
        self.start_button.configure(state="disabled")
        self.demo_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")
        self.cancel_button.grid()
        self.notebook.select(self.log_tab)
        self.status.set("Starting analysis")
        self.report_button.configure(state="disabled")
        self.append_log("\nResults: " + str(output) + "\nCommand:\n" + json.dumps(command, ensure_ascii=False) + "\n")
        token = self.current_token
        threading.Thread(target=self.read_worker, args=(process, token), daemon=True).start()

    def read_worker(self, process, token):
        try:
            for line in process.stdout:
                self.events.put(("log", token, line))
            code = process.wait()
            self.events.put(("done", token, code))
        except Exception as exc:
            self.events.put(("log", token, f"Worker reader: {exc}\n"))
            self.events.put(("done", token, process.poll() if process.poll() is not None else 1))

    def poll(self):
        if self.worker_log and self.worker_log.exists():
            try:
                with self.worker_log.open("r", encoding="utf-8", errors="replace") as handle:
                    handle.seek(self.log_offset)
                    chunk = handle.read()
                    self.log_offset = handle.tell()
                if chunk:
                    self.append_log(chunk)
            except OSError:
                pass
        if self.process is not None and self.output:
            try:
                data = json.loads((self.output / "progress.json").read_text(encoding="utf-8"))
                phase = PHASE_NAMES.get(data.get("phase"), data.get("phase", "Starting"))
                self.status.set("Cancelling" if self.cancelled else phase)
                parts = [f"Elapsed {float(data.get('elapsed_seconds', 0)):.1f} s"]
                for key, label in (("candidate_count", "Candidates"), ("records_scanned", "Transcripts scanned"), ("nucleotides_scanned", "Nucleotides scanned"), ("binding_hits_counted", "Binding hits")):
                    if data.get(key) is not None:
                        parts.append(f"{label} {data[key]:,}")
                self.progress_detail.set("  ·  ".join(parts))
            except (OSError, ValueError, TypeError):
                pass
        while not self.events.empty():
            kind, token, value = self.events.get_nowait()
            if token != self.current_token:
                continue
            if kind == "log":
                self.append_log(value)
            else:
                self.finish(value)
        if not self.closing:
            self.root.after(300, self.poll)

    def finish(self, exit_code):
        self.process = None
        self.progress.stop()
        self.start_button.configure(state="normal")
        self.demo_button.configure(state="normal")
        self.cancel_button.configure(state="disabled")
        self.cancel_button.grid_remove()
        success = exit_code == 0 and not self.cancelled and (self.output / "result.json").is_file()
        outcome = {"exit_code": exit_code, "cancelled": self.cancelled, "output": str(self.output)}
        if success:
            try:
                self.result_data = json.loads((self.output / "result.json").read_text(encoding="utf-8"))
                self.populate_result()
                outcome["candidate_count"] = len(self.result_data.get("candidates", []))
                self.status.set(f"Complete · {outcome['candidate_count']:,} candidates")
                self.report_button.configure(state="normal")
                self.notebook.select(self.result_tab)
                if self.open_browser:
                    self.open_result()
            except (OSError, ValueError, TypeError) as exc:
                success = False
                outcome["error"] = str(exc)
                self.status.set("Unable to read results; files retained")
        else:
            self.status.set("Cancelled · files retained" if self.cancelled else "Analysis stopped · see Activity")
        if self.worker_log:
            (self.worker_log.parent / "desktop_outcome.json").write_text(json.dumps(outcome, ensure_ascii=False, indent=2), encoding="utf-8")
        if self.on_finished:
            self.on_finished(success, outcome)
        if self.closing:
            self.root.destroy()

    def cancel(self):
        process = self.process
        if process is None or process.poll() is not None:
            return
        self.cancelled = True
        self.status.set("Cancelling · files retained")
        process.terminate()
        self.root.after(3000, lambda: process.kill() if process.poll() is None else None)

    def close(self):
        from tkinter import messagebox
        if self.process is not None:
            if not messagebox.askyesno("Exit", "Stop the current analysis and exit? Saved files will be retained.", parent=self.root):
                return
            self.closing = True
            self.cancel()
            # Keep the event poll alive until the worker exit is observed.
            self.root.after(3500, self.poll)
        else:
            self.root.destroy()

    def populate_result(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        candidates = self.result_data.get("candidates", [])
        self.candidate_lookup = {}
        for index, candidate in enumerate(candidates[:300]):
            iid = str(index)
            self.candidate_lookup[iid] = candidate
            p = candidate.get("p_unpaired_target")
            exact = candidate.get("other_gene_exact_compatible", candidate.get("other_gene_exact_ag"))
            score = candidate.get("ranking_score")
            self.tree.insert("", "end", iid=iid, values=(f"{candidate.get('cut_after_1based', '?')} | {candidate.get('cut_after_1based', 0) + 1}",
                              f"{candidate.get('left_length', '?')} / {candidate.get('right_length', '?')}",
                              "Not assessed" if p is None else f"{p:.2%}", "Not assessed" if exact is None else exact,
                              "—" if score is None else f"{score:.1f}", candidate.get("ranking_stage", "unknown").replace("_", " ").title()))
        metadata = self.result_data.get("metadata", {})
        reference = metadata.get("offtarget", {}).get("reference_label") or metadata.get("scope_label") or "See report"
        self.result_scope.set(f"{len(candidates):,} candidates · showing {min(300, len(candidates)):,}\nReference: {reference}")
        if self.tree.get_children():
            original = next((key for key, row in self.candidate_lookup.items() if row.get("start_1based") == 401 and row.get("end_1based") == 416), self.tree.get_children()[0])
            self.tree.selection_set(original)
            self.tree.see(original)
        self.draw_profile()

    def draw_profile(self):
        self.canvas.delete("all")
        if not self.result_data:
            self.canvas.create_text(20, 55, anchor="w", fill=UI_COLORS["muted"], text="RNA accessibility along your transcript", font=("Segoe UI", 10))
            return
        width = max(self.canvas.winfo_width(), 300)
        length = len(self.result_data.get("target", {}).get("sequence", "")) or 1
        left, right, top, bottom = 48, width - 22, 22, 114
        self.canvas.create_line(left, bottom, right, bottom, fill=UI_COLORS["line"])
        self.canvas.create_line(left, top, left, bottom, fill=UI_COLORS["line"])
        for value in (0, .5, 1):
            y = bottom - value * (bottom - top)
            self.canvas.create_text(left - 8, y, anchor="e", text=f"{value:g}", fill=UI_COLORS["muted"])
        selected = self.tree.selection()
        selected_id = self.candidate_lookup.get(selected[0], {}).get("candidate_id") if selected else None
        rows = self.result_data.get("candidates", [])
        stride = max(1, math.ceil(len(rows) / 2500))
        for candidate in rows[::stride]:
            p = candidate.get("p_unpaired_target")
            x = left + (right - left) * candidate.get("cut_after_1based", 0) / length
            y = bottom + 9 if p is None else bottom - min(1, max(0, p)) * (bottom - top)
            active = candidate.get("candidate_id") == selected_id
            radius, color = (5, "#CA7855") if active else (2.5, UI_COLORS["teal"] if p is not None else "#A9B6B6")
            self.canvas.create_oval(x - radius, y - radius, x + radius, y + radius, fill=color, outline="")
        self.canvas.create_text(left, 144, anchor="w", text="1", fill=UI_COLORS["muted"])
        self.canvas.create_text(right, 144, anchor="e", text=f"{length} nt", fill=UI_COLORS["muted"])
        self.canvas.create_text((left + right) / 2, 144, text="Transcript position (1-based)" + (" · sampled display" if stride > 1 else ""), fill=UI_COLORS["muted"])

    def select_nearest(self, event):
        if not self.candidate_lookup or not self.result_data:
            return
        width = max(self.canvas.winfo_width(), 300)
        length = len(self.result_data.get("target", {}).get("sequence", "")) or 1
        position = (event.x - 48) / max(1, width - 70) * length
        key = min(self.candidate_lookup, key=lambda iid: abs(self.candidate_lookup[iid].get("cut_after_1based", 0) - position))
        self.tree.selection_set(key)
        self.tree.see(key)

    def select_candidate(self, event=None):
        selected = self.tree.selection()
        if not selected:
            return
        row = self.candidate_lookup[selected[0]]
        description = [f"{row.get('candidate_id')} · target {row.get('start_1based')}–{row.get('end_1based')}",
                       f"Target RNA 5′→3′: {row.get('target_rna')}",
                       f"Left DNA arm:  {row.get('left_arm_dna')}\nRight DNA arm: {row.get('right_arm_dna')}"]
        if row.get("left_partzyme") or row.get("right_partzyme"):
            description += [f"Left partzyme:  {row.get('left_partzyme')}", f"Right partzyme: {row.get('right_partzyme')}"]
        else:
            description.append("Output: binding arms only")
        if row.get("length_tradeoff_score") is not None:
            reasons = row.get("length_tradeoff_reasons") or []
            reason_labels = {
                "arm_length_outside_user_preferred_range": "outside preferred arm range",
                "arm_length_outside_user_practical_range": "outside practical arm range",
                "very_strong_1M_reference_duplex_stability_penalized_not_a_release_rate": "strong reference duplex penalty",
                "strong_binding_penalty_not_assessed_missing_reference_energy": "reference duplex penalty not assessed",
                "within_user_preferences_no_excess_reference_binding_penalty": "within arm preferences",
            }
            note = "; ".join(reason_labels.get(reason, str(reason)) for reason in reasons) if reasons else "within preferences"
            description.append(f"Arm tradeoff factor: {row['length_tradeoff_score']:.3f} · {note}")
        self.selected_detail.set("\n".join(description))
        self.draw_profile()

    def open_path(self, path):
        from tkinter import messagebox
        if not path or not path.exists():
            messagebox.showinfo("No result available", "Run an analysis first. The selected output is not available yet.", parent=self.root)
            return
        if os.name == "nt" and path.is_dir():
            os.startfile(str(path))
        else:
            webbrowser.open(path.resolve().as_uri())

    def open_result(self):
        self.open_path(self.output / "viewer" / "index.html" if self.output else None)

    def open_output(self):
        self.open_path(self.output)

    def open_af3(self):
        self.open_path(self.output / "alphafold3_inputs" if self.output else None)

    def open_structure(self):
        local = self.output / "viewer" / "structure.html" if self.output else None
        self.open_path(local if local and local.exists() else resource_root() / "ui" / "structure.html")


def main(argv=None):
    import tkinter as tk
    p = argparse.ArgumentParser(description="Local Zyme Site Designer desktop")
    p.add_argument("--results-dir", type=Path, help="Optional persistent output-root override")
    p.add_argument("--smoke-test", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--smoke-demo", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--smoke-layout", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--smoke-report", type=Path, help=argparse.SUPPRESS)
    p.add_argument("--smoke-max-seconds", type=float, default=90, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    dpi = set_windows_dpi_awareness()
    root = tk.Tk()
    if args.smoke_test or args.smoke_demo or args.smoke_layout:
        root.withdraw()
    status = {"tk_version": root.tk.call("info", "patchlevel"), "dpi_awareness": dpi, "resource_root": str(resource_root()), "frozen": bool(getattr(sys, "frozen", False)), "completed": False}

    def save_and_exit(success, outcome):
        status.update({"completed": True, "success": success, **outcome})
        if args.smoke_report:
            args.smoke_report.parent.mkdir(parents=True, exist_ok=True)
            args.smoke_report.write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
        root.after(100, root.destroy)

    app = DesignerApp(root, args.results_dir, open_browser=not (args.smoke_test or args.smoke_demo or args.smoke_layout), on_finished=save_and_exit if args.smoke_demo else None)
    root.update_idletasks()
    if args.smoke_layout:
        root.overrideredirect(True)
        cases = []
        for scaling in (4 / 3, 5 / 3, 2.0):
            root.tk.call("tk", "scaling", scaling)
            for width, height in ((900, 600), (1024, 680), (1280, 800)):
                root.geometry(f"{width}x{height}+32000+32000")
                root.deiconify()
                for tab in (app.target_tab, app.param_tab, app.log_tab, app.result_tab):
                    app.notebook.select(tab)
                    root.update_idletasks()
                    root.update()
                    button = app.start_button
                    x = button.winfo_rootx() - root.winfo_rootx()
                    y = button.winfo_rooty() - root.winfo_rooty()
                    visible = bool(button.winfo_ismapped() and x >= 0 and y >= 0
                                   and x + button.winfo_width() <= root.winfo_width()
                                   and y + button.winfo_height() <= root.winfo_height())
                    cases.append({"scaling": scaling, "size": [width, height],
                                  "tab": app.notebook.tab(tab, "text"), "run_button_visible": visible})
        root.withdraw()
        status["layout_cases"] = cases
        save_and_exit(all(case["run_button_visible"] for case in cases), {"primary_button_label": app.start_button.cget("text")})
    elif args.smoke_demo:
        root.after(100, app.run_demo)

        def timeout():
            if not status["completed"]:
                app.cancel()
                save_and_exit(False, {"error": "smoke demo timeout"})
        root.after(int(args.smoke_max_seconds * 1000), timeout)
    elif args.smoke_test:
        status["tabs"] = len(app.notebook.tabs())
        status["worker_running"] = app.process is not None
        root.after(500, lambda: save_and_exit(True, {"widgets_created": True}))
    root.mainloop()
    return 0 if not (args.smoke_test or args.smoke_demo or args.smoke_layout) or status.get("success") else 1
