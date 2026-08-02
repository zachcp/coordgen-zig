#!/usr/bin/env python3
"""Compute a Git tree object ID from a single-root source tar archive."""

from __future__ import annotations

import hashlib
import io
import sys
import tarfile


def object_id(kind: bytes, payload: bytes) -> bytes:
    return hashlib.sha1(kind + b" " + str(len(payload)).encode() + b"\0" + payload).digest()


def tree_id(node: dict[bytes, object]) -> bytes:
    entries: list[tuple[bytes, bytes]] = []
    for name, value in node.items():
        if isinstance(value, dict):
            mode = b"40000"
            digest = tree_id(value)
            sort_name = name + b"/"
        else:
            mode, payload = value
            digest = object_id(b"blob", payload)
            sort_name = name
        encoded = mode + b" " + name + b"\0" + digest
        entries.append((sort_name, encoded))
    body = b"".join(encoded for _, encoded in sorted(entries, key=lambda entry: entry[0]))
    return object_id(b"tree", body)


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} ARCHIVE", file=sys.stderr)
        return 2

    root: dict[bytes, object] = {}
    with tarfile.open(sys.argv[1], "r:*") as archive:
        members = [member for member in archive.getmembers() if member.name.rstrip("/")]
        roots = {member.name.rstrip("/").split("/", 1)[0] for member in members}
        if len(roots) != 1:
            raise SystemExit("archive must contain exactly one top-level directory")
        prefix = next(iter(roots)) + "/"

        for member in members:
            if member.isdir() or member.name.rstrip("/") == prefix.rstrip("/"):
                continue
            if not member.name.startswith(prefix):
                raise SystemExit(f"entry escapes archive root: {member.name}")
            relative = member.name[len(prefix) :].rstrip("/")
            parts = [part.encode() for part in relative.split("/") if part]
            if not parts:
                continue
            parent = root
            for part in parts[:-1]:
                child = parent.setdefault(part, {})
                if not isinstance(child, dict):
                    raise SystemExit(f"file/directory collision at {relative}")
                parent = child

            if member.issym():
                value = (b"120000", member.linkname.encode())
            elif member.isreg():
                source = archive.extractfile(member)
                if source is None:
                    raise SystemExit(f"cannot read {member.name}")
                mode = b"100755" if member.mode & 0o111 else b"100644"
                value = (mode, source.read())
            else:
                raise SystemExit(f"unsupported tar member type: {member.name}")
            parent[parts[-1]] = value

    print(tree_id(root).hex())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
