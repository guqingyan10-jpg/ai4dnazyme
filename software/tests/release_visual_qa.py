"""Capture this application's native widgets without displaying a visible UI.

A fully transparent, disabled window must occupy a drawable viewport because
off-screen Tk children are clipped by Windows. PrintWindow renders only this
application's client into a private bitmap; desktop/other-app capture is never
used. A predominantly black image is rejected rather than shipped as evidence.
"""
from __future__ import annotations

import argparse
import ctypes
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from zyme_designer.desktop import DesignerApp, DesktopConfig, SCAFFOLD_BUILTIN, resource_root, set_windows_dpi_awareness


def native_capture(widget, flags=3, method="printwindow"):
    """Ask this window to paint its own client; never capture the desktop."""
    from ctypes import wintypes
    from PIL import Image
    user, gdi = ctypes.windll.user32, ctypes.windll.gdi32
    user.GetDC.argtypes = [wintypes.HWND]
    user.GetDC.restype = wintypes.HDC
    user.ReleaseDC.argtypes = [wintypes.HWND, wintypes.HDC]
    user.PrintWindow.argtypes = [wintypes.HWND, wintypes.HDC, wintypes.UINT]
    user.PrintWindow.restype = wintypes.BOOL
    user.SendMessageW.argtypes = [wintypes.HWND, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM]
    user.SendMessageW.restype = wintypes.LPARAM
    gdi.CreateCompatibleDC.argtypes = [wintypes.HDC]
    gdi.CreateCompatibleDC.restype = wintypes.HDC
    gdi.CreateCompatibleBitmap.argtypes = [wintypes.HDC, ctypes.c_int, ctypes.c_int]
    gdi.CreateCompatibleBitmap.restype = wintypes.HBITMAP
    gdi.SelectObject.argtypes = [wintypes.HDC, wintypes.HGDIOBJ]
    gdi.SelectObject.restype = wintypes.HGDIOBJ
    gdi.DeleteDC.argtypes = [wintypes.HDC]
    gdi.DeleteObject.argtypes = [wintypes.HGDIOBJ]

    class Header(ctypes.Structure):
        _fields_ = [("size", wintypes.DWORD), ("width", wintypes.LONG), ("height", wintypes.LONG),
                    ("planes", wintypes.WORD), ("bitcount", wintypes.WORD), ("compression", wintypes.DWORD),
                    ("image_size", wintypes.DWORD), ("x", wintypes.LONG), ("y", wintypes.LONG),
                    ("colors", wintypes.DWORD), ("important", wintypes.DWORD)]

    class Info(ctypes.Structure):
        _fields_ = [("header", Header), ("palette", wintypes.DWORD * 3)]

    hwnd = widget.winfo_id()
    width, height = widget.winfo_width(), widget.winfo_height()
    source = user.GetDC(hwnd)
    dc = gdi.CreateCompatibleDC(source)
    bitmap = gdi.CreateCompatibleBitmap(source, width, height)
    previous = gdi.SelectObject(dc, bitmap)
    try:
        if method == "printwindow":
            if not user.PrintWindow(hwnd, dc, flags):
                raise OSError("PrintWindow did not render the native Tk client")
        else:
            # WM_PRINT / PRF_CLIENT | PRF_ERASEBKGND | PRF_CHILDREN.
            user.SendMessageW(hwnd, 0x0317, dc, 0x04 | 0x08 | 0x10)
        info = Info()
        info.header = Header(ctypes.sizeof(Header), width, -height, 1, 32, 0, width * height * 4, 0, 0, 0, 0)
        buffer = ctypes.create_string_buffer(width * height * 4)
        gdi.GetDIBits.argtypes = [wintypes.HDC, wintypes.HBITMAP, wintypes.UINT, wintypes.UINT,
                                ctypes.c_void_p, ctypes.POINTER(Info), wintypes.UINT]
        if gdi.GetDIBits(dc, bitmap, 0, height, buffer, ctypes.byref(info), 0) != height:
            raise OSError("Could not read native client bitmap")
        return Image.frombuffer("RGB", (width, height), buffer.raw, "raw", "BGRX", 0, 1)
    finally:
        gdi.SelectObject(dc, previous)
        gdi.DeleteObject(bitmap)
        gdi.DeleteDC(dc)
        user.ReleaseDC(hwnd, source)


def black_fraction(image):
    from PIL import ImageChops
    red, green, blue = image.convert("RGB").split()
    maximum = ImageChops.lighter(ImageChops.lighter(red, green), blue)
    return sum(maximum.histogram()[:3]) / (image.width * image.height)


def main():
    import tkinter as tk
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    set_windows_dpi_awareness()
    root = tk.Tk()
    root.withdraw()
    # Apply transparency and disable input before the window is mapped. The
    # application is never visibly displayed and cannot intercept interaction.
    root.attributes("-alpha", 0.0)
    root.attributes("-disabled", True)
    root.overrideredirect(True)
    app = DesignerApp(root, results_parent=args.out / "unused_results", open_browser=False)
    examples = resource_root() / "examples"
    app.set_config(DesktopConfig(target_fasta=str(examples / "BIRC5_NM_001168_3.fasta"), gene="BIRC5",
                                 transcriptome=str(examples / "tiny_reference_NOT_HUMAN_TRANSCRIPTOME.fasta"),
                                 reference_label="QA tiny reference, NOT full human transcriptome",
                                 scaffold_mode=SCAFFOLD_BUILTIN, output_parent=str(args.out / "unused_results")))
    app.output = args.result.parent
    app.result_data = json.loads(args.result.read_text(encoding="utf-8"))
    app.populate_result()
    app.report_button.configure(state="normal")
    app.status.set(f"Complete · {len(app.result_data['candidates']):,} candidates")
    root.geometry("1600x1080+0+0")
    root.deiconify()
    report = []
    try:
        for tab, name in ((app.target_tab, "Desktop_Transcript_EN.png"),
                          (app.param_tab, "Desktop_Design_EN.png"),
                          (app.result_tab, "Desktop_Results_EN.png")):
            app.notebook.select(tab)
            root.update_idletasks()
            root.update()
            image = native_capture(root)
            low, high = image.convert("L").getextrema()
            if high - low < 50 or black_fraction(image) > .70:
                raise RuntimeError("The native window did not render; no verified screenshot produced")
            image.save(args.out / name)
            button = app.start_button
            report.append({"file": name, "pixels": list(image.size), "run_button_mapped": bool(button.winfo_ismapped()),
                           "tab": app.notebook.tab(tab, "text"),
                           "black_fraction": black_fraction(image),
                           "capture": "PrintWindow PW_CLIENTONLY | PW_RENDERFULLCONTENT; native application client only",
                           "window_alpha": float(root.attributes("-alpha")),
                           "window_input_disabled": bool(root.attributes("-disabled"))})
    finally:
        for callback in root.tk.splitlist(root.tk.call("after", "info")):
            root.after_cancel(callback)
        root.destroy()
    (args.out / "visual_capture.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
