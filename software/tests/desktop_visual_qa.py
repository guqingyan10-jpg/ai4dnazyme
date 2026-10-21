"""Optional native-window visual QA. Requires Pillow only in the QA environment.

python tests/desktop_visual_qa.py --result PATH/result.json --out PATH/qa
The controlled window closes after capturing three tabs. No worker is started.
"""
from __future__ import annotations

import argparse
import ctypes
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zyme_designer.desktop import DesignerApp, demo_config


def native_capture(widget):
    """PrintWindow fallback for Windows sessions without a capturable desktop."""
    from ctypes import wintypes
    from PIL import Image
    user, gdi = ctypes.windll.user32, ctypes.windll.gdi32
    user.GetDC.restype = wintypes.HDC
    user.GetDC.argtypes = [wintypes.HWND]
    gdi.CreateCompatibleDC.restype = wintypes.HDC
    gdi.CreateCompatibleDC.argtypes = [wintypes.HDC]
    gdi.CreateCompatibleBitmap.restype = wintypes.HBITMAP
    gdi.CreateCompatibleBitmap.argtypes = [wintypes.HDC, ctypes.c_int, ctypes.c_int]
    gdi.SelectObject.restype = wintypes.HGDIOBJ
    gdi.SelectObject.argtypes = [wintypes.HDC, wintypes.HGDIOBJ]
    user.PrintWindow.argtypes = [wintypes.HWND, wintypes.HDC, wintypes.UINT]
    user.ReleaseDC.argtypes = [wintypes.HWND, wintypes.HDC]
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
        if not user.PrintWindow(hwnd, dc, 1):
            raise OSError("PrintWindow did not render the native Tk client")
        info = Info()
        info.header = Header(ctypes.sizeof(Header), width, -height, 1, 32, 0, width * height * 4, 0, 0, 0, 0)
        buffer = ctypes.create_string_buffer(width * height * 4)
        gdi.GetDIBits.argtypes = [wintypes.HDC, wintypes.HBITMAP, wintypes.UINT, wintypes.UINT, ctypes.c_void_p, ctypes.POINTER(Info), wintypes.UINT]
        if gdi.GetDIBits(dc, bitmap, 0, height, buffer, ctypes.byref(info), 0) != height:
            raise OSError("Could not read native client bitmap")
        return Image.frombuffer("RGB", (width, height), buffer.raw, "raw", "BGRX", 0, 1)
    finally:
        gdi.SelectObject(dc, previous)
        gdi.DeleteObject(bitmap)
        gdi.DeleteDC(dc)
        user.ReleaseDC(hwnd, source)


def main():
    import tkinter as tk
    from PIL import ImageGrab
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if sys.platform == "win32":
        ctypes.windll.user32.SetProcessDPIAware()
    root = tk.Tk()
    root.geometry("1120x860+15+15")
    app = DesignerApp(root, results_parent=args.out / "results_not_created", open_browser=False)
    app.set_config(demo_config(args.out / "results_not_created"))
    app.output = args.result.parent
    app.result_data = json.loads(args.result.read_text(encoding="utf-8"))
    app.populate_result()
    root.attributes("-topmost", True)
    args.out.mkdir(parents=True, exist_ok=True)
    try:
        for tab, name in ((app.target_tab, "Desktop_Target.png"), (app.param_tab, "Desktop_Parameters.png"), (app.result_tab, "Desktop_Results.png")):
            app.notebook.select(tab)
            root.update_idletasks()
            root.update()
            time.sleep(.35)
            x, y = root.winfo_rootx(), root.winfo_rooty()
            try:
                image = ImageGrab.grab(bbox=(x, y, x + root.winfo_width(), y + root.winfo_height()))
            except OSError:
                image = native_capture(root)
            if sum(1 for pixel in image.convert("RGB").getdata() if max(pixel) < 3) > image.width * image.height * .7:
                raise OSError("Native screenshot is predominantly unrendered black; use an unlocked interactive Windows session. No verified screenshot saved.")
            image.save(args.out / name)
            print(args.out / name)
    finally:
        root.destroy()


if __name__ == "__main__":
    main()
