# Sources of Third-Party License Texts

The original inventory was checked on 2026-10-06. Upstream license texts are retained, with line endings normalized where recorded; their terms have not been translated or rewritten. Microsoft runtime text was extracted from the linked DOCX documents, whose original numbering and formatting remain authoritative. This inventory documents provenance, not legal advice or independent redistribution authorization.

This distribution contains source code only. References to an inspected Windows build explain the provenance of retained notices; they do not mean that its executable, DLLs or Python runtime are included here.

| Local file | Version and recorded source |
|---|---|
| `ViennaRNA_2.7.2_COPYING.txt` | Full COPYING text from the ViennaRNA 2.7.2 wheel; checked against the [versioned upstream file](https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/COPYING). |
| `ViennaRNA_2.7.2_NAVIEW_notice.txt` | Complete opening copyright/license comment from the [versioned naview.c](https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/src/ViennaRNA/plotting/naview/naview.c), without the algorithm implementation. |
| `PyInstaller_6.22.3_COPYING.txt` | Full license from the versioned wheel, including GPL, the bootloader exception and Apache-2.0 material. |
| `pyahocorasick_2.3.1_LICENSE.txt` | Full LICENSE from the versioned wheel; its installation metadata also identifies BSD-3-Clause and Public-Domain material. |
| `CPython_3.13.9_LICENSE.txt` | Local Python 3.13.9 license checked against the [upstream version](https://raw.githubusercontent.com/python/cpython/v3.13.9/LICENSE), including its license chain. |
| `CPython_3.13.9_Windows_distribution_LICENSE.txt` | Complete license from the inspected python.org Windows installation, including Microsoft distributable-code conditions and additional component notices. Recorded sources: [installer](https://www.python.org/ftp/python/3.13.9/python-3.13.9-amd64.exe) and [installer SBOM](https://www.python.org/ftp/python/3.13.9/python-3.13.9-amd64.exe.spdx.json). |
| `Tcl_8.6.15_license.terms.txt` | Full Tcl 8.6.15 distribution license, checked against the [versioned upstream file](https://raw.githubusercontent.com/tcltk/tcl/core-8-6-15/license.terms). |
| `Tk_8.6.15_license.terms.txt` | Full local Tk license, checked against [Tk 8.6.15](https://raw.githubusercontent.com/tcltk/tk/core-8-6-15/license.terms). |
| `OpenSSL_3.0.18_LICENSE.txt` | Full text from the [Python dependency tag](https://raw.githubusercontent.com/python/cpython-source-deps/openssl-3.0.18/LICENSE.txt). The inspected runtime reported OpenSSL 3.0.18. |
| `OpenSSL_3.0.18_COPYRIGHT_notice.txt` | Complete Copyright section of the [OpenSSL 3.0.18 README](https://raw.githubusercontent.com/openssl/openssl/openssl-3.0.18/README.md), without product-description text. |
| `libffi_3.4.4_LICENSE.txt` | Full [versioned Python dependency license](https://raw.githubusercontent.com/python/cpython-source-deps/libffi-3.4.4/LICENSE); version recorded from the installer SBOM. |
| `zlib_1.3.1_LICENSE.txt` | Full [versioned dependency license](https://raw.githubusercontent.com/python/cpython-source-deps/zlib-1.3.1/LICENSE). The inspected runtime reported zlib 1.3.1. |
| `xz_5.2.5_COPYING.txt` | Full [COPYING](https://raw.githubusercontent.com/python/cpython-source-deps/xz-5.2.5/COPYING). The inspected build embedded that version of liblzma, not the xz CLI or shell scripts. |
| `mpdecimal_4.0.0_COPYRIGHT.txt` | Full [versioned COPYRIGHT](https://raw.githubusercontent.com/python/cpython-source-deps/mpdecimal-4.0.0/COPYRIGHT.txt). The inspected runtime reported libmpdec 4.0.0. |
| `HACL_cpython_3.13.9_generated_C_notice.txt` | Complete opening MIT notice from [Hacl_Hash_SHA2.c](https://raw.githubusercontent.com/python/cpython/v3.13.9/Modules/_hacl/Hacl_Hash_SHA2.c), without implementation code. |
| `HACL_cpython_3.13.9_Kremlin_runtime_notice.txt` | Complete opening notices from [target.h](https://raw.githubusercontent.com/python/cpython/v3.13.9/Modules/_hacl/include/krml/internal/target.h) and [FStar_UInt128_Verified.h](https://raw.githubusercontent.com/python/cpython/v3.13.9/Modules/_hacl/include/krml/FStar_UInt128_Verified.h). |
| `HACL_cpython_3.13.9_LICENSE.txt` | Full [Apache license at the recorded HACL* commit](https://raw.githubusercontent.com/hacl-star/hacl-star/bb3d0dc8d9d15a5cd51094d5b69e70aa09005ff0/LICENSE). Generated C-code MIT notices are listed separately. |
| `libb2_0.98.1_COPYING.txt` | Full [libb2 0.98.1 CC0 text](https://raw.githubusercontent.com/BLAKE2/libb2/v0.98.1/COPYING). The [CPython reference implementation](https://raw.githubusercontent.com/python/cpython/v3.13.9/Modules/_blake2/impl/blake2.h) also records CC0. This does not assert a separately loaded libb2 DLL. |
| `BLAKE2_cpython_3.13.9_notice.txt` | Complete opening reference-implementation, Samuel Neves and CC0 notices from the [CPython source](https://raw.githubusercontent.com/python/cpython/v3.13.9/Modules/_blake2/impl/blake2.h), without algorithm code. |
| `Microsoft_Visual_C_2015_2022_runtime_terms.txt` | Text extracted from the [original DOCX](https://visualstudio.microsoft.com/wp-content/uploads/2021/09/Visual-C-Runtime-2015-2022-License-1.docx) linked by the [Microsoft license entry](https://visualstudio.microsoft.com/license-terms/vs2022-cruntime/). The inspected DLL versions were 14.42/14.44. End-user terms alone are not independent redistribution authorization. |
| `Expat_2.7.3_COPYING.txt` | Full [upstream versioned COPYING](https://raw.githubusercontent.com/libexpat/libexpat/R_2_7_3/expat/COPYING), retained as supplementary provenance; the inspected final runtime did not include the relevant Expat extension modules. |
| `Microsoft_Visual_C_V14_2026_runtime_terms.txt` | Text from the [original DOCX](https://visualstudio.microsoft.com/wp-content/uploads/2025/10/Visual-C-V14-License-Redistributable_and_Runtime_ENU.docx) linked by the [Microsoft license entry](https://visualstudio.microsoft.com/license-terms/vs2026-ga-visualcpp-v14-redist-runtime/), retained as supplementary provenance. It does not assert that the corresponding 14.51 DLL is included. |

The recorded build used python.org Python 3.13.9 and Tcl/Tk 8.6.15. The corresponding copyright/license material is also present in the Windows distribution license. The retained ViennaRNA, PyInstaller and pyahocorasick texts were compared against the installed versioned files after line-ending normalization. No original third-party component was modified by this project.

The runtime provenance record identifies the sources and signatures of the inspected Microsoft runtime files. Those engineering checks do not replace the distributing team's license obligations; consult `../THIRD_PARTY_NOTICES.md` and the original terms.

`WINDOWS_BUNDLE_COMPONENTS.json` records 25 DLL/extension files from that inspected build, with hashes and source categories. It is not an inventory of this source-only ZIP, does not cover all source/documentation files and is not a complete SPDX SBOM.
