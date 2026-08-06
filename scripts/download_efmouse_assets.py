#!/usr/bin/env python3
"""Download the official EFMouse tetrahedral mesh without installing AllenSDK."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


FILES = {
    "mesh": (
        "https://raw.githubusercontent.com/klabhub/EFMouse/main/aux_files/"
        "Alekseichuk_mesh_clean.mat",
        "data/efmouse/Alekseichuk_mesh_clean.mat",
    ),
}


def download(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": "ti-vpl-fem/3.0"}
    if existing:
        headers["Range"] = f"bytes={existing}-"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        accepts_resume = response.status == 206 and existing > 0
        mode = "ab" if accepts_resume else "wb"
        if not accepts_resume:
            existing = 0
        total_header = response.headers.get("Content-Length")
        total = existing + int(total_header) if total_header else None
        written = existing
        with partial.open(mode) as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                written += len(chunk)
                if total:
                    print(f"\r{destination.name}: {written / total:6.1%}", end="", flush=True)
    print()
    if partial.stat().st_size < 10_000_000:
        raise RuntimeError("Downloaded file is unexpectedly small; retained as .part for inspection")
    partial.replace(destination)
    sha256 = hashlib.sha256()
    with destination.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha256.update(chunk)
    return {"url": url, "path": str(destination), "bytes": destination.stat().st_size, "sha256": sha256.hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root = Path(args.project_root)
    manifest = []
    for url, relative in FILES.values():
        destination = root / relative
        if destination.exists() and destination.stat().st_size > 10_000_000 and not args.force:
            print(f"Already present: {destination}")
            continue
        manifest.append(download(url, destination))
    manifest_path = root / "data/efmouse/download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
