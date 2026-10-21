# Third-party notices

The project's own source code is licensed under the [MIT License](LICENSE). Dependencies and bundled runtime components retain their original licenses; they are not relicensed by that statement. Preserve this file and the full texts in `third_party_licenses/` with redistributed builds.

This document accompanies the **1.0.0 research-software release**. The component versions and source audit below describe the archived **2026-10-06 v0.2 r3 Windows build**. Use the current build's component and acceptance records for release-specific inventory verification. The historical inventory is not a claim that a new build has been independently audited.

## Components and original terms

| Component | Observed version | License and distribution notes | Full text |
|---|---|---|---|
| ViennaRNA / RNAlib Python bindings | 2.7.2 | ViennaRNA terms allow research, educational and commercial use/modification with attribution and no redistribution fee beyond media costs. Contact its authors before inclusion in a commercial product. | [COPYING](third_party_licenses/ViennaRNA_2.7.2_COPYING.txt) |
| NAVIEW within ViennaRNA | Included with 2.7.2 source | Separate Robert E. Bruccoleri terms require noncommercial copying, credit and retained notices; other copying requires specific permission. This notice is retained conservatively with RNAlib. | [NAVIEW notice](third_party_licenses/ViennaRNA_2.7.2_NAVIEW_notice.txt) |
| CPython Windows interpreter | 3.13.9 | PSF License Version 2 and historical license chain; the Windows distribution also includes Microsoft distributable-code conditions and component notices. | [Windows license](third_party_licenses/CPython_3.13.9_Windows_distribution_LICENSE.txt), [source license chain](third_party_licenses/CPython_3.13.9_LICENSE.txt) |
| Tcl / Tk | 8.6.15 / 8.6.15 | Retain the original copyright and license notices with distribution. | [Tcl](third_party_licenses/Tcl_8.6.15_license.terms.txt), [Tk](third_party_licenses/Tk_8.6.15_license.terms.txt) |
| PyInstaller | 6.22.3 | GPL-2.0-or-later with the bootloader exception; some hooks/support code use Apache-2.0. The exception does not remove the applicable terms when distributing PyInstaller itself or modifications to it. | [License and exception](third_party_licenses/PyInstaller_6.22.3_COPYING.txt) |
| pyahocorasick | 2.3.1; optional | BSD-3-Clause; retain copyright, conditions and disclaimer. The standard-library fallback preserves the screening policy. | [LICENSE](third_party_licenses/pyahocorasick_2.3.1_LICENSE.txt) |
| OpenSSL | 3.0.18 | Apache-2.0 and retained OpenSSL Project / original author notices. | [LICENSE](third_party_licenses/OpenSSL_3.0.18_LICENSE.txt), [copyright](third_party_licenses/OpenSSL_3.0.18_COPYRIGHT_notice.txt) |
| libffi | 3.4.4 | MIT; retain the original author/copyright notices. | [LICENSE](third_party_licenses/libffi_3.4.4_LICENSE.txt) |
| zlib | 1.3.1 | zlib license; retain notices and do not misrepresent authorship. | [LICENSE](third_party_licenses/zlib_1.3.1_LICENSE.txt) |
| bzip2 used by `_bz2.pyd` | 1.0.8 | Retain Julian R. Seward's copyright, conditions and disclaimer. | [Windows distribution license](third_party_licenses/CPython_3.13.9_Windows_distribution_LICENSE.txt) |
| liblzma used by `_lzma.pyd` | XZ 5.2.5 | That version's liblzma is public domain; other XZ programs/scripts have separate terms. | [COPYING](third_party_licenses/xz_5.2.5_COPYING.txt) |
| libmpdec used by `_decimal.pyd` | 4.0.0 | BSD-2-Clause; retain Stefan Krah's notices. | [COPYRIGHT](third_party_licenses/mpdecimal_4.0.0_COPYRIGHT.txt) |
| HACL* code and support headers in CPython | CPython 3.13.9 distribution | Generated C code uses MIT; Kremlin support headers use Apache-2.0. Keep both sets of notices. | [C-code notice](third_party_licenses/HACL_cpython_3.13.9_generated_C_notice.txt), [support notice](third_party_licenses/HACL_cpython_3.13.9_Kremlin_runtime_notice.txt), [Apache text](third_party_licenses/HACL_cpython_3.13.9_LICENSE.txt) |
| BLAKE2 reference implementation in CPython | CPython 3.13.9 distribution | CC0-1.0; retain Samuel Neves' source attribution. | [Notice](third_party_licenses/BLAKE2_cpython_3.13.9_notice.txt), [CC0 text](third_party_licenses/libb2_0.98.1_COPYING.txt) |
| Microsoft C/C++ runtime | MSVCP140 14.44.35211.0; VCRUNTIME140 / `_1` 14.42.34438.0 | Proprietary Microsoft runtime terms, not MIT or PSF. The runtime end-user EULA is not evidence of independent redistribution authorization. | [Python Windows conditions](third_party_licenses/CPython_3.13.9_Windows_distribution_LICENSE.txt), [runtime terms](third_party_licenses/Microsoft_Visual_C_2015_2022_runtime_terms.txt) |

We acknowledge the ViennaRNA authors and the Institute for Theoretical Chemistry of the University of Vienna. NAVIEW copying is by permission of Robert E. Bruccoleri under the retained notice. Third-party authors and institutions do not endorse this application.

See the [ViennaRNA license](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/license.html) and [PyInstaller license explanation](https://pyinstaller.org/en/v6.22.0/license.html). The complete versioned license texts take precedence over this summary.

## Historical runtime provenance

The archived r3 build used python.org Python 3.13.9. Tcl/Tk, Python extensions, OpenSSL, libffi, zlib and VCRUNTIME came from that official Windows installer; ViennaRNA and pyahocorasick came from versioned wheels. `MSVCP140.dll` came from the local Windows System32 directory, not the Python installer. Signatures and source hashes were checked as provenance evidence.

The dated [component inventory](third_party_licenses/WINDOWS_BUNDLE_COMPONENTS.json) records 25 observed DLL/extension files and their hashes. It is an observed binary inventory, not a complete source-level SPDX SBOM or a legal opinion. The [Python installer SBOM](https://www.python.org/ftp/python/3.13.9/python-3.13.9-amd64.exe.spdx.json) describes its installer, not this application.

That build relied on Windows system UCRT instead of copying `ucrtbase.dll` or `api-ms-win-*`. It did not include an Anaconda runtime, libheif, Pillow, SQLite or Expat extension modules. Some additional historic license texts remain for provenance; a notice's presence does not establish that its binary is bundled.

## Redistribution checklist

- Preserve all original license texts, copyrights, disclaimers and the PyInstaller exception.
- Distinguish project-source MIT licensing from dependency licenses and scientific-data provenance.
- Check the actual frozen inventory for the build being distributed.
- Confirm the team's Microsoft runtime redistribution qualification and the applicable [Visual Studio REDIST terms](https://learn.microsoft.com/en-us/visualstudio/releases/2022/redistribution) before public redistribution. Source location or a valid signature does not establish authorization; see [Microsoft distribution guidance](https://learn.microsoft.com/en-us/cpp/windows/determining-which-dlls-to-redistribute?view=msvc-170).
- Observe the ViennaRNA/NAVIEW distribution conditions. Do not infer permission for paid redistribution or commercial-product integration from the project's MIT license.
- Record any public repository or download publication only after it has actually occurred.

The earlier same-machine sanitized-PATH tests are documented separately from independent clean-machine testing; no independent VM certification is asserted. Source URLs are recorded in [SOURCES.md](third_party_licenses/SOURCES.md).

This is an engineering attribution and distribution checklist, not legal advice. The distributing team is responsible for compliance.
