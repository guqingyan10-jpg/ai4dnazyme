"""Native layout regressions for the English desktop release.

These checks create off-screen Tk windows only. No worker is launched, browser
opened, sequence analysed, or persistent application output created. On systems
without a usable Tk display the native checks are explicitly skipped.
"""
from __future__ import annotations

from pathlib import Path
import re
import unittest

from zyme_designer import desktop


CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")


class DesktopLayoutReleaseTests(unittest.TestCase):
    def create_app(self, scaling=4 / 3):
        try:
            import tkinter as tk
            root = tk.Tk()
        except (ImportError, RuntimeError) as exc:
            self.skipTest(f"Native Tk is unavailable: {exc}")
        except Exception as exc:
            # TclError is intentionally handled without importing tkinter at
            # module scope, so backend-only test discovery remains portable.
            if type(exc).__name__ == "TclError":
                self.skipTest(f"No native Tk display: {exc}")
            raise
        root.withdraw()
        root.overrideredirect(True)
        root.tk.call("tk", "scaling", scaling)
        self.addCleanup(self.destroy_root, root)
        app = desktop.DesignerApp(
            root,
            results_parent=Path(__file__).resolve().parents[1] / "results_not_created_by_layout_test",
            open_browser=False,
        )
        return root, app

    @staticmethod
    def destroy_root(root):
        # Tk shares its event queue across interpreter instances. Cancel the
        # application's scheduled poll before destroying an off-screen root,
        # avoiding stale "invalid command name ...poll" diagnostic noise.
        for callback in root.tk.splitlist(root.tk.call("after", "info")):
            root.after_cancel(callback)
        root.destroy()

    def assert_primary_action_in_client(self, root, app):
        button = app.start_button
        self.assertTrue(button.winfo_ismapped(), "The primary action was removed from the visible layout")
        x = button.winfo_rootx() - root.winfo_rootx()
        y = button.winfo_rooty() - root.winfo_rooty()
        width, height = button.winfo_width(), button.winfo_height()
        self.assertGreater(width, 40, "The primary action has no usable width")
        self.assertGreater(height, 20, "The primary action has no usable height")
        self.assertGreaterEqual(x, 0)
        self.assertGreaterEqual(y, 0)
        self.assertLessEqual(x + width, root.winfo_width(), "The primary action overflows the client horizontally")
        self.assertLessEqual(y + height, root.winfo_height(), "The primary action is hidden below the client")
        self.assertIsNone(app.process, "A layout check must never start an analysis worker")

    def test_primary_action_at_small_window_and_high_scaling(self):
        # Tk points-to-pixels factors approximately correspond to 100%, 125%
        # and 150% Windows scaling. This tests geometry, not monitor hardware.
        for scaling in (4 / 3, 5 / 3, 2.0):
            root, app = self.create_app(scaling)
            for width, height in ((900, 600), (1024, 680), (1280, 800)):
                with self.subTest(scaling=scaling, size=f"{width}x{height}"):
                    root.geometry(f"{width}x{height}+32000+32000")
                    root.deiconify()
                    root.update_idletasks()
                    root.update()
                    self.assertEqual((root.winfo_width(), root.winfo_height()), (width, height),
                                     "The minimum window size exceeds the tested desktop viewport")
                    for tab in (app.target_tab, app.param_tab, app.log_tab, app.result_tab):
                        app.notebook.select(tab)
                        root.update_idletasks()
                        root.update()
                        self.assert_primary_action_in_client(root, app)
            root.withdraw()

    def test_advanced_settings_cannot_hide_primary_action(self):
        root, app = self.create_app(2.0)
        self.assertTrue(callable(getattr(app, "toggle_advanced", None)),
                        "Advanced settings should be available without crowding the main workflow")
        root.geometry("900x600+32000+32000")
        root.deiconify()
        app.notebook.select(app.param_tab)
        app.toggle_advanced()
        root.update_idletasks()
        root.update()
        self.assert_primary_action_in_client(root, app)
        app.toggle_advanced()
        root.update_idletasks()
        root.update()
        self.assert_primary_action_in_client(root, app)

    def test_public_labels_and_status_are_english(self):
        root, app = self.create_app()
        labels = [root.title(), app.status.get(), app.progress_detail.get(),
                  app.result_scope.get(), app.selected_detail.get()]
        labels += [str(app.notebook.tab(tab, "text")) for tab in app.notebook.tabs()]

        def collect_text(widget):
            if "text" in widget.keys():
                labels.append(str(widget.cget("text")))
            for child in widget.winfo_children():
                collect_text(child)

        collect_text(root)
        untranslated = [text for text in labels if CJK.search(text)]
        self.assertFalse(untranslated,
                         f"English-release UI labels are still in Chinese: {untranslated[:8]!r}")
        self.assertRegex(str(app.start_button.cget("text")).lower(), r"\b(run|design)\b")


if __name__ == "__main__":
    unittest.main()
