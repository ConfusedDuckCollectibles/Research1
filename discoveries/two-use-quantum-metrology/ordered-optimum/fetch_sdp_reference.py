#!/usr/bin/env python3
"""Download only the pinned author dataset, checking it before writing."""
from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen

URL = ("https://raw.githubusercontent.com/qiushi-liu/strategies_in_metrology/"
       "07748a4e94d7b29dc888d2cbf41d0b0f421c58dc/data/QFI_seq_rangep_0_1_N_2.txt")
SHA = "a73a6686f3e305b4500119c96e861f89504de0ab7d432348b8737fd7f0e35002"
if __name__ == "__main__":
    with urlopen(URL, timeout=20) as response:
        data = response.read(5228)
    if len(data) != 5227 or sha256(data).hexdigest() != SHA:
        raise SystemExit("Unexpected reference dataset")
    dst = Path(__file__).resolve().parent / "source" / "QFI_seq_rangep_0_1_N_2.txt"
    dst.parent.mkdir(exist_ok=True)
    dst.write_bytes(data)
    print(f"Verified reference SHA-256: {SHA}")
