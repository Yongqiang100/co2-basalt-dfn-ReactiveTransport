# -*- coding: utf-8 -*-
"""
generate_xdmf.py -- Create ParaView-compatible XDMF + clean HDF5

PFLOTRAN HDF5 output has special characters in group/variable names
(spaces, brackets, plus signs) that ParaView's XDMF reader cannot handle.
This script creates a clean HDF5 with sanitized names + matching XDMF.

Usage:
    python generate_xdmf.py pflotran_results/p32_100_s42
    python generate_xdmf.py pflotran_results/p32_100_s42 --inp
    python generate_xdmf.py --all
    python generate_xdmf.py --all --inp
"""

import os
import sys
import re
import glob
import argparse
import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("ERROR: h5py required.  pip install h5py")

SEC_PER_YEAR = 3.156e7


def sanitize_name(name):
    """Clean variable name for ParaView compatibility.
    'Al(OH)2+ [M]' -> 'Al_OH_2p_M'
    'Calcite VF [m^3 mnrl_m^3 bulk]' -> 'Calcite_VF'
    """
    # Extract just the variable name (before the unit bracket)
    match = re.match(r'^(.*?)(?:\s*\[.*\])?\s*$', name)
    if match:
        clean = match.group(1).strip()
    else:
        clean = name.strip()
    # Replace special characters
    clean = clean.replace("+++", "3p")
    clean = clean.replace("++", "2p")
    clean = clean.replace("+", "p")
    clean = clean.replace("--", "2m")
    clean = clean.replace("-", "m")
    clean = clean.replace("(", "_")
    clean = clean.replace(")", "_")
    clean = clean.replace(" ", "_")
    clean = clean.replace("^", "")
    clean = clean.replace("/", "_per_")
    # Remove double underscores
    while "__" in clean:
        clean = clean.replace("__", "_")
    clean = clean.strip("_")
    return clean


def find_time_groups(h5):
    """Find and sort time groups in PFLOTRAN HDF5."""
    results = []
    for k in h5.keys():
        if "Time" not in k:
            continue
        if "failure" in k.lower() or "cut" in k.lower():
            continue
        parts = k.strip().split()
        try:
            time_idx = parts.index("Time")
            t_val = float(parts[time_idx + 1])
            if len(parts) > time_idx + 2 and parts[time_idx + 2] == "y":
                t_yr = t_val
            else:
                t_yr = t_val / SEC_PER_YEAR
            results.append((k, t_yr))
        except (ValueError, IndexError):
            continue
    results.sort(key=lambda x: x[1])
    return results


def read_uge(uge_path):
    """Return cell centres and volumes from .uge."""
    xs, ys, zs, vols = [], [], [], []
    with open(uge_path) as f:
        header = f.readline().strip().split()
        n_cells = int(header[1])
        for _ in range(n_cells):
            parts = f.readline().strip().split()
            xs.append(float(parts[1]))
            ys.append(float(parts[2]))
            zs.append(float(parts[3]))
            vols.append(float(parts[4]))
    coords = np.column_stack([xs, ys, zs])
    return coords, np.array(vols), n_cells


def read_inp(inp_path):
    """Read LaGriT AVS/UCD .inp file."""
    with open(inp_path) as f:
        header = f.readline().strip().split()
        n_nodes = int(header[0])
        n_cells = int(header[1])
        coords = np.zeros((n_nodes, 3))
        for i in range(n_nodes):
            parts = f.readline().strip().split()
            coords[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
        cells = []
        for _ in range(n_cells):
            parts = f.readline().strip().split()
            cell_type = parts[2].lower()
            node_ids = [int(x) - 1 for x in parts[3:]]
            cells.append((cell_type, node_ids))
    return coords, cells, n_nodes, n_cells


def inp_cell_xdmf_type(cell_type):
    mapping = {
        "tet": ("Tetrahedron", 4),
        "hex": ("Hexahedron", 8),
        "pri": ("Wedge", 6),
        "pyr": ("Pyramid", 5),
        "tri": ("Triangle", 3),
        "quad": ("Quadrilateral", 4),
    }
    return mapping.get(cell_type, (None, None))


def create_clean_h5(src_path, dst_path, n_expected):
    """Copy PFLOTRAN HDF5 data into a clean file with sanitized names.
    Returns list of (clean_group_name, t_yr, {clean_var: orig_var}) or None."""

    with h5py.File(src_path, "r") as src:
        time_groups = find_time_groups(src)
        if not time_groups:
            return None

        # Build variable name mapping from first time group
        first_grp = src[time_groups[0][0]]
        var_map = {}  # clean_name -> original_name
        for orig_name in first_grp.keys():
            ds = first_grp[orig_name]
            if ds.ndim == 1 and ds.shape[0] == n_expected:
                clean = sanitize_name(orig_name)
                var_map[clean] = orig_name

        print(f"  Variables: {len(var_map)} (of {len(first_grp.keys())} total)")

        with h5py.File(dst_path, "w") as dst:
            step_info = []
            for i, (tg_name, t_yr) in enumerate(time_groups):
                clean_grp_name = f"Step_{i:04d}"
                grp_src = src[tg_name]
                grp_dst = dst.create_group(clean_grp_name)
                grp_dst.attrs["time_yr"] = t_yr

                for clean_name, orig_name in var_map.items():
                    if orig_name in grp_src:
                        data = grp_src[orig_name][:]
                        grp_dst.create_dataset(clean_name, data=data)

                step_info.append((clean_grp_name, t_yr, var_map))

        return step_info


def write_xdmf_points(xdmf_path, h5_basename, coords, n_cells, step_info):
    """Write XDMF for point-cloud with clean HDF5 references."""
    lines = []
    lines.append('<?xml version="1.0" ?>')
    lines.append('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    lines.append('<Xdmf Version="3.0">')
    lines.append('  <Domain>')
    lines.append('    <Grid Name="TimeSeries" GridType="Collection" '
                  'CollectionType="Temporal">')

    for grp_name, t_yr, var_map in step_info:
        lines.append(f'      <Grid Name="mesh" GridType="Uniform">')
        lines.append(f'        <Time Value="{t_yr:.6e}" />')
        lines.append(f'        <Topology TopologyType="Polyvertex" '
                     f'NumberOfElements="{n_cells}" />')
        lines.append(f'        <Geometry GeometryType="XYZ">')
        lines.append(f'          <DataItem Dimensions="{n_cells} 3" '
                     f'NumberType="Float" Precision="8" Format="XML">')
        for i in range(n_cells):
            lines.append(f'            {coords[i, 0]:.6e} '
                         f'{coords[i, 1]:.6e} {coords[i, 2]:.6e}')
        lines.append(f'          </DataItem>')
        lines.append(f'        </Geometry>')

        for clean_name in sorted(var_map.keys()):
            lines.append(
                f'        <Attribute Name="{clean_name}" '
                f'AttributeType="Scalar" Center="Node">')
            lines.append(
                f'          <DataItem Dimensions="{n_cells}" '
                f'NumberType="Float" Precision="8" '
                f'Format="HDF">{h5_basename}:/{grp_name}/{clean_name}'
                f'</DataItem>')
            lines.append(f'        </Attribute>')
        lines.append(f'      </Grid>')

    lines.append('    </Grid>')
    lines.append('  </Domain>')
    lines.append('</Xdmf>')

    with open(xdmf_path, "w") as f:
        f.write("\n".join(lines))


def write_xdmf_inp(xdmf_path, h5_basename, coords, cells, n_nodes,
                    n_data, step_info):
    """Write XDMF for tetrahedral mesh with clean HDF5 references."""
    type_counts = {}
    for ctype, _ in cells:
        type_counts[ctype] = type_counts.get(ctype, 0) + 1
    dominant_type = max(type_counts, key=type_counts.get)
    xdmf_type, nodes_per_cell = inp_cell_xdmf_type(dominant_type)

    if xdmf_type is None:
        return False

    filtered = [c for c in cells if c[0] == dominant_type]
    n_cells_used = len(filtered)

    if n_cells_used != len(cells):
        print(f"  NOTE: Using {n_cells_used}/{len(cells)} cells "
              f"(type={dominant_type})")

    conn = np.zeros((n_cells_used, nodes_per_cell), dtype=np.int64)
    for i, (_, nids) in enumerate(filtered):
        conn[i, :] = nids[:nodes_per_cell]

    if n_data == n_cells_used:
        centre = "Cell"
    elif n_data == n_nodes:
        centre = "Node"
    else:
        centre = "Cell"

    lines = []
    lines.append('<?xml version="1.0" ?>')
    lines.append('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    lines.append('<Xdmf Version="3.0">')
    lines.append('  <Domain>')
    lines.append('    <Grid Name="TimeSeries" GridType="Collection" '
                  'CollectionType="Temporal">')

    for grp_name, t_yr, var_map in step_info:
        lines.append(f'      <Grid Name="mesh" GridType="Uniform">')
        lines.append(f'        <Time Value="{t_yr:.6e}" />')

        lines.append(
            f'        <Topology TopologyType="{xdmf_type}" '
            f'NumberOfElements="{n_cells_used}">')
        lines.append(
            f'          <DataItem Dimensions="{n_cells_used} '
            f'{nodes_per_cell}" NumberType="Int" Format="XML">')
        for i in range(n_cells_used):
            lines.append("            " +
                         " ".join(str(x) for x in conn[i]))
        lines.append(f'          </DataItem>')
        lines.append(f'        </Topology>')

        lines.append(f'        <Geometry GeometryType="XYZ">')
        lines.append(
            f'          <DataItem Dimensions="{n_nodes} 3" '
            f'NumberType="Float" Precision="8" Format="XML">')
        for i in range(n_nodes):
            lines.append(f'            {coords[i, 0]:.6e} '
                         f'{coords[i, 1]:.6e} {coords[i, 2]:.6e}')
        lines.append(f'          </DataItem>')
        lines.append(f'        </Geometry>')

        for clean_name in sorted(var_map.keys()):
            lines.append(
                f'        <Attribute Name="{clean_name}" '
                f'AttributeType="Scalar" Center="{centre}">')
            lines.append(
                f'          <DataItem Dimensions="{n_data}" '
                f'NumberType="Float" Precision="8" '
                f'Format="HDF">{h5_basename}:/{grp_name}/{clean_name}'
                f'</DataItem>')
            lines.append(f'        </Attribute>')
        lines.append(f'      </Grid>')

    lines.append('    </Grid>')
    lines.append('  </Domain>')
    lines.append('</Xdmf>')

    with open(xdmf_path, "w") as f:
        f.write("\n".join(lines))
    return True


def process_dfn(dfn_dir, use_inp=False):
    """Generate clean HDF5 + XDMF for a single DFN."""
    dfn_name = os.path.basename(dfn_dir.rstrip("/"))

    h5_files = sorted(glob.glob(os.path.join(dfn_dir, "pflotran_co2.h5")))
    if not h5_files:
        print(f"  [{dfn_name}] No pflotran_co2.h5 found, skipping")
        return False

    src_h5 = h5_files[0]
    uge_path = os.path.join(dfn_dir, "full_mesh.uge")
    inp_path = os.path.join(dfn_dir, "full_mesh.inp")

    if not os.path.exists(uge_path):
        print(f"  [{dfn_name}] No full_mesh.uge found, skipping")
        return False

    # Read mesh for cell count
    coords_uge, vols, n_cells_uge = read_uge(uge_path)

    # Output files
    clean_h5_name = "pflotran_clean.h5"
    clean_h5_path = os.path.join(dfn_dir, clean_h5_name)
    xdmf_path = os.path.join(dfn_dir, "pflotran_co2.xmf")

    # Create clean HDF5
    print(f"  [{dfn_name}] Creating clean HDF5...")
    step_info = create_clean_h5(src_h5, clean_h5_path, n_cells_uge)
    if step_info is None:
        print(f"  [{dfn_name}] No time groups found")
        return False

    print(f"  [{dfn_name}] Time steps: {len(step_info)}")
    print(f"  [{dfn_name}] Wrote: {clean_h5_path}")

    # Try INP mesh
    if use_inp and os.path.exists(inp_path):
        print(f"  [{dfn_name}] Reading INP mesh...")
        try:
            coords_inp, cells, n_nodes, n_cells_inp = read_inp(inp_path)
            print(f"  [{dfn_name}] Mesh: {n_nodes} nodes, "
                  f"{n_cells_inp} cells")
            ok = write_xdmf_inp(xdmf_path, clean_h5_name, coords_inp,
                                cells, n_nodes, n_cells_uge, step_info)
            if ok:
                print(f"  [{dfn_name}] Wrote: {xdmf_path} (tetrahedral)")
                return True
            else:
                print(f"  [{dfn_name}] INP failed, falling back to UGE")
        except Exception as e:
            print(f"  [{dfn_name}] INP error: {e}, falling back to UGE")

    # Point cloud fallback
    print(f"  [{dfn_name}] Writing point cloud XDMF...")
    write_xdmf_points(xdmf_path, clean_h5_name, coords_uge,
                       n_cells_uge, step_info)
    print(f"  [{dfn_name}] Wrote: {xdmf_path} (point cloud)")
    print(f"  [{dfn_name}] TIP: In ParaView use 'Point Gaussian' or "
          f"'Delaunay 3D' filter")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate ParaView-compatible XDMF + clean HDF5 "
                    "from PFLOTRAN output")
    parser.add_argument("dfn_dir", nargs="?", default=None,
                        help="Path to a single DFN result directory")
    parser.add_argument("--all", action="store_true",
                        help="Process all DFNs in pflotran_results/")
    parser.add_argument("--results", type=str,
                        default=os.path.join(os.getcwd(), "pflotran_results"),
                        help="Root results directory")
    parser.add_argument("--inp", action="store_true",
                        help="Use full_mesh.inp for tetrahedral mesh")
    args = parser.parse_args()

    if args.all:
        dfn_dirs = sorted(glob.glob(os.path.join(args.results, "*")))
        dfn_dirs = [d for d in dfn_dirs if os.path.isdir(d)]
    elif args.dfn_dir:
        dfn_dirs = [args.dfn_dir]
    else:
        parser.print_help()
        sys.exit(1)

    print(f"Processing {len(dfn_dirs)} DFN(s)...\n")

    success = 0
    for dfn_dir in dfn_dirs:
        if process_dfn(dfn_dir, use_inp=args.inp):
            success += 1
        print()

    print(f"Done. Generated {success}/{len(dfn_dirs)} XDMF files.")
    print(f"\nOpen in ParaView:")
    print(f"  File -> Open -> pflotran_co2.xmf")
    print(f"\nNote: XDMF references pflotran_clean.h5 (not pflotran_co2.h5)")
    print(f"      Both files must be in the same directory.")


if __name__ == "__main__":
    main()