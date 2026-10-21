"""Restore the losslessly compressed browser example before testing/building."""
from hashlib import sha256
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
target = ROOT / "software/ui/demo-data.js"
expected = "4ce98792f88352b9d4a05da987c6ef913a4cc810a7f11dcc45b744754f74503c"
payload = gzip.decompress(target.with_suffix(".js.gz").read_bytes())
assert len(payload) == 5898693 and sha256(payload).hexdigest() == expected
if target.exists():
    assert target.read_bytes() == payload, "Existing example differs; no overwrite performed"
else:
    target.write_bytes(payload)
print("Complete browser example restored and SHA256 verified.")
