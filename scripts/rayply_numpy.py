# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Reads and writes ray cloud .ply files (binary little-endian, float or double precision) as numpy record arrays,
with the fields named in the file's header (e.g. x, y, z, time, nx, ny, nz or rayx, rayy, rayz, red, green, blue,
alpha). Shared by the study and plotting scripts, which import it from their own directory."""
import numpy as np

PLY_TYPES = {"float": "<f4", "float32": "<f4", "double": "<f8", "float64": "<f8", "uchar": "u1", "uint8": "u1",
             "char": "i1", "int8": "i1", "short": "<i2", "int16": "<i2", "ushort": "<u2", "uint16": "<u2",
             "int": "<i4", "int32": "<i4", "uint": "<u4", "uint32": "<u4"}


def read_header(f):
    """Reads a .ply header from the binary file f, returning the header bytes, the vertex record type and the number
    of vertices. Raises ValueError if the header is truncated or isn't a binary little-endian vertex list"""
    lines = []
    while True:
        line = f.readline()
        if not line:
            raise ValueError(f"{f.name}: truncated .ply header (no end_header)")
        lines.append(line)
        if line.strip() == b"end_header":
            break
    text = [line.decode("ascii", "replace").split() for line in lines]
    if not text[0] or text[0][0] != "ply":
        raise ValueError(f"{f.name}: not a .ply file")
    formats = [t for t in text if t and t[0] == "format"]
    if len(formats) != 1 or len(formats[0]) < 2 or formats[0][1] != "binary_little_endian":
        raise ValueError(f"{f.name}: only binary little-endian .ply files are supported")
    elements = [t for t in text if t and t[0] == "element"]
    if len(elements) != 1 or len(elements[0]) != 3 or elements[0][1] != "vertex":
        raise ValueError(f"{f.name}: expected a single vertex element, as in a ray cloud")
    fields = []
    for t in text:
        if t and t[0] == "property":
            if len(t) != 3 or t[1] not in PLY_TYPES:
                raise ValueError(f"{f.name}: unsupported property: {' '.join(t)}")
            fields.append((t[2], PLY_TYPES[t[1]]))
    return b"".join(lines), np.dtype(fields), int(elements[0][2])


def read_ply(name):
    """Reads a ray cloud .ply file, returning its header (bytes) and its vertices as a numpy record array"""
    with open(name, "rb") as f:
        header, dtype, count = read_header(f)
        data = np.fromfile(f, dtype=dtype, count=count)
    if len(data) != count:
        raise ValueError(f"{name}: expected {count} vertices, found {len(data)}")
    return header, data


def write_ply(name, header, data):
    """Writes vertices (a record array of the header's vertex type) to a .ply file, with the given header from
    read_ply, updating its vertex count (zero padded to the same width, as raycloudtools writes it)"""
    lines = header.split(b"\n")
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 3 and parts[0] == b"element" and parts[1] == b"vertex":
            lines[i] = b"element vertex %0*d" % (len(parts[2]), len(data))
    with open(name, "wb") as f:
        f.write(b"\n".join(lines))
        data.tofile(f)


def read_rays(name):
    """Reads a ray cloud .ply file, returning its ray starts, ends (as float64 N x 3 arrays), times and colours
    (N x 4: red, green, blue, alpha)"""
    _, data = read_ply(name)
    ray = ("nx", "ny", "nz") if "nx" in data.dtype.names else ("rayx", "rayy", "rayz")
    ends = np.stack([data["x"], data["y"], data["z"]], axis=1).astype(np.float64)
    starts = ends + np.stack([data[r] for r in ray], axis=1)  # the ray vector is start - end
    colours = np.stack([data["red"], data["green"], data["blue"], data["alpha"]], axis=1)
    return starts, ends, data["time"], colours
