import os
import re

class DLGParser:
    def __init__(self, dlg_file_path):
        self.dlg_file_path = dlg_file_path
        self.poses = []
        self.top_poses = []

    def parse_ranking_table(self):
        ranking_pattern = r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\-\d.]+)\s+([\-\d.]+).*RANKING"
        try:
            with open(self.dlg_file_path, "r") as file:
                lines = file.readlines()
            for line in lines:
                match = re.match(ranking_pattern, line)
                if match:
                    rank, sub_rank, run, binding_energy, cluster_rmsd = match.groups()
                    self.poses.append({
                        "rank": int(rank),
                        "sub_rank": int(sub_rank),
                        "run": int(run),
                        "binding_energy": float(binding_energy),
                        "cluster_rmsd": float(cluster_rmsd),
                    })
            print(f"Parsed {len(self.poses)} poses from the ranking table.")
        except Exception as e:
            print(f"Error parsing DLG file: {e}")

    def get_top_scoring_pose(self, top_n=1):
        try:
            self.top_poses = sorted(self.poses, key=lambda x: (x["rank"], x["sub_rank"]))[:top_n]
            if not self.top_poses:
                print("No poses found. Check your DLG file.")
            else:
                print(f"Retrieved top {len(self.top_poses)} pose(s).")
        except Exception as e:
            print(f"Error retrieving top pose(s): {e}")

    def extract_pose_coordinates(self, run_number, include_flexible=False):
        pose_coordinates = []
        recording = False
        try:
            with open(self.dlg_file_path, "r") as file:
                for line in file:
                    if recording:
                        if line.startswith("DOCKED: HETATM"):
                            # Process ligand (non-flexible) lines.
                            pdb_line = self.format_pdb_line(line[8:].strip(), flexible=False)
                            pose_coordinates.append(pdb_line)
                        elif include_flexible and line.startswith("DOCKED: ATOM"):
                            # Process flexible residue lines.
                            pdb_line = self.format_pdb_line(line[8:].strip(), flexible=True)
                            pose_coordinates.append(pdb_line)
                        elif line.startswith("DOCKED: ENDMDL"):
                            print(f"Finished reading coordinates for run {run_number}.")
                            break
                    elif re.match(f"DOCKED:\\s*MODEL\\s*{run_number}\\s*", line):
                        recording = True
                        print(f"Started reading coordinates for run {run_number}.")
            if not pose_coordinates:
                print(f"No coordinates found for run {run_number}.")
        except Exception as e:
            print(f"Error extracting coordinates for run {run_number}: {e}")
        return pose_coordinates

    @staticmethod
    def format_pdb_line(docked_line, flexible=False):
        """
        Process a dlg line by preserving its original fields (including spacing) and 
        then modifying only a few positions:
        - Force the record (columns 1-6) to "ATOM  "
        - For flexible atoms (DOCKED: ATOM lines), use the chain from the dlg, but if blank, force "X".
            For non‑flexible (ligand) lines, force the chain to "L".
        - Force the occupancy (columns 55–60) to 1.00 and temperature factor (columns 61–66) to 0.00.
        
        This method assumes the dlg line is fixed‐width. It pads the line to 80 characters if needed.
        """
        # Remove any ending newline and pad the line to 80 characters.
        orig_line = docked_line.rstrip("\n").ljust(80)
        
        # (1) Force the record label (columns 1–6) to "ATOM  ":
        #    PDB columns are 1-indexed; in Python 0-index, columns [0:6] are the record.
        new_line = "ATOM  " + orig_line[6:]
        
        # (2) Modify the chain field, which is in column 22 (index 21).
        if flexible:
            # Use the chain from the dlg if available; if blank, force "X" (or any default matching your receptor).
            chain = new_line[21:22]
            if chain.strip() == "":
                new_line = new_line[:21] + "X" + new_line[22:]
        else:
            # For ligand (non-flexible) atoms, force chain to "L".
            new_line = new_line[:21] + "L" + new_line[22:]
        
        # (3) Overwrite occupancy (columns 55–60, Python indices [54:60]) and
        #     temperature factor (columns 61–66, Python indices [60:66]) with fixed values.
        occ = f"{1.00:6.2f}"   # creates a 6-character string (e.g., " 1.00")
        temp = f"{0.00:6.2f}"  # creates a 6-character string (e.g., " 0.00")
        new_line = new_line[:54] + occ + temp + new_line[66:]
        
        return new_line


    def write_complex_pdb(self, pose, receptor_file_path, output_filename, include_flexible=True):
        """
        Merge the dlg flexible atoms into the rigid receptor pdbqt file by inserting them
        into the corresponding residue block, based on the residue number and chain.
        """
        try:
            # Read receptor pdbqt file lines.
            with open(receptor_file_path, "r") as rec_file:
                receptor_lines = [line.rstrip("\n") for line in rec_file]

            # Extract dlg coordinates from the dlg file.
            pose_coordinates = self.extract_pose_coordinates(pose["run"], include_flexible=include_flexible)
            if not pose_coordinates:
                print("No dlg coordinates extracted; aborting merge.")
                return

            # Separate flexible (ATOM) and ligand (HETATM) lines.
            dlg_flexible_lines = [line for line in pose_coordinates if line.startswith("ATOM")]
            ligand_lines = [line for line in pose_coordinates if line.startswith("HETATM")]

            # Group flexible dlg atoms by residue key. We use (chain, residue number) as the key.
            flexible_groups = {}
            for line in dlg_flexible_lines:
                if len(line) >= 26:
                    key = (line[21], line[22:26].strip())
                    flexible_groups.setdefault(key, []).append(line)
            
            # Process the receptor file by grouping contiguous ATOM/HETATM lines into residue blocks.
            new_lines = []
            current_block = []
            current_key = None

            def flush_block(block, key):
                """Return the block, inserting flexible atoms if the residue key is found."""
                result = block[:]
                if key in flexible_groups:
                    # Insert flexible atoms immediately after this residue block.
                    result.extend(flexible_groups[key])
                    del flexible_groups[key]
                return result

            for line in receptor_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Extract the key from the receptor line.
                    if len(line) >= 26:
                        key = (line[21], line[22:26].strip())
                    else:
                        key = None
                    if current_key is None:
                        current_key = key
                        current_block.append(line)
                    elif key == current_key:
                        current_block.append(line)
                    else:
                        new_lines.extend(flush_block(current_block, current_key))
                        current_block = [line]
                        current_key = key
                else:
                    if current_block:
                        new_lines.extend(flush_block(current_block, current_key))
                        current_block = []
                        current_key = None
                    new_lines.append(line)
            if current_block:
                new_lines.extend(flush_block(current_block, current_key))
            
            # Any remaining flexible groups that did not match will be appended at the end with a warning.
            if flexible_groups:
                for key, lines_group in flexible_groups.items():
                    print(f"Warning: Flexible residue {key} not found in receptor; appending its atoms at the end.")
                    new_lines.extend(lines_group)
            
            # Append ligand atoms (typically placed at the end).
            if ligand_lines:
                new_lines.extend(ligand_lines)
            
            # Reassign atom serial numbers sequentially.
            final_lines = []
            atom_serial = 1
            for line in new_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    new_line = f"{line[:6]}{atom_serial:5d}" + line[11:]
                    final_lines.append(new_line)
                    atom_serial += 1
                else:
                    final_lines.append(line)
            
            # Write the merged pdbqt file.
            with open(output_filename, "w") as outfile:
                outfile.write("\n".join(final_lines) + "\n")
            print(f"Complex pdbqt file written: {output_filename}")
        except Exception as e:
            print(f"Error writing complex pdbqt file: {e}")

    def write_pdb(self, pose, output_dir=".", pose_number=1, include_flexible=False):
        try:
            coordinates = self.extract_pose_coordinates(pose["run"], include_flexible=include_flexible)
            if not coordinates:
                print(f"No coordinates to write for pose {pose['rank']}.")
                return
            base_name = os.path.splitext(os.path.basename(self.dlg_file_path))[0]
            out_filename = os.path.join(output_dir, f"{base_name}_pose_{pose_number}.pdb")
            with open(out_filename, "w") as pdb_file:
                pdb_file.write("\n".join(coordinates) + "\n")
            print(f"PDB file written: {out_filename}")
        except Exception as e:
            print(f"Error writing PDB file: {e}")

    def write_top_poses(self, output_dir=".", top_n=1, include_flexible=False, receptor_file=None):
        try:
            self.get_top_scoring_pose(top_n=top_n)
            for idx, pose in enumerate(self.top_poses, start=1):
                if receptor_file:
                    base_name = os.path.splitext(os.path.basename(self.dlg_file_path))[0]
                    # For receptor merging, write out as pdbqt first.
                    out_filename = os.path.join(output_dir, f"{base_name}_pose_{idx}_complex.pdbqt")
                    self.write_complex_pdb(pose, receptor_file, out_filename, include_flexible=True)
                else:
                    self.write_pdb(pose, output_dir=output_dir, pose_number=idx, include_flexible=include_flexible)
        except Exception as e:
            print(f"Error writing top poses: {e}")

    def print_summary(self, log_file=None):
        summary = ""
        if not self.top_poses:
            summary = "No top poses available."
        else:
            summary_lines = ["Top Scoring Poses:"]
            for pose in self.top_poses:
                summary_lines.append(
                    f"Rank: {pose['rank']}, Sub-Rank: {pose['sub_rank']}, Run: {pose['run']}, Binding Energy: {pose['binding_energy']} kcal/mol"
                )
            summary = "\n".join(summary_lines)
        print(summary)
        if log_file:
            try:
                with open(log_file, "w") as lf:
                    lf.write(summary + "\n")
                print(f"Summary written to {log_file}")
            except Exception as e:
                print(f"Error writing summary log: {e}")


# -------------------------------------------------------------------
# New function: Use MDAnalysis to guess connectivity and write final PDB.
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_bonds

def pdbqt_to_pdb_with_conect(input_pdbqt, output_pdb):
    """
    Read a pdbqt file (merged output), guess bonds using MDAnalysis,
    and write out a final PDB file with CONECT records.
    """
    try:
        u = mda.Universe(input_pdbqt)
    except Exception as e:
        print(f"Error reading {input_pdbqt}: {e}")
        return

    temp_pdb = "temp_no_connect.pdb"
    try:
        u.atoms.write(temp_pdb)
    except Exception as e:
        print(f"Error writing temporary PDB: {e}")
        return

    try:
        bonds = guess_bonds(u.atoms)
    except Exception as e:
        print(f"Error guessing bonds: {e}")
        bonds = []

    bonds_dict = {}
    for bond in bonds:
        i = bond[0] + 1  # Convert 0-index to 1-index.
        j = bond[1] + 1
        bonds_dict.setdefault(i, []).append(j)
        bonds_dict.setdefault(j, []).append(i)

    try:
        with open(temp_pdb, "r") as f:
            pdb_lines = f.readlines()
    except Exception as e:
        print(f"Error reading temporary PDB: {e}")
        return

    conect_lines = []
    for i in sorted(bonds_dict):
        partners = sorted(bonds_dict[i])
        # Write multiple CONECT lines if more than four bonds.
        for idx in range(0, len(partners), 4):
            chunk = partners[idx:idx+4]
            line = "CONECT%5d" % i
            for partner in chunk:
                line += "%5d" % partner
            conect_lines.append(line + "\n")

    new_pdb_lines = []
    inserted = False
    for line in pdb_lines:
        if line.startswith("END"):
            new_pdb_lines.extend(conect_lines)
            new_pdb_lines.append(line)
            inserted = True
        else:
            new_pdb_lines.append(line)
    if not inserted:
        new_pdb_lines.extend(conect_lines)
        new_pdb_lines.append("END\n")

    try:
        with open(output_pdb, "w") as f:
            f.write("".join(new_pdb_lines))
        print(f"Final PDB with connectivity written to: {output_pdb}")
    except Exception as e:
        print(f"Error writing final PDB: {e}")
    finally:
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)


# -------------------------------------------------------------------
# Main function and argument parsing.
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Process AutoDock DLG files, output PDB/pdbqt files, merge with receptor if provided, and optionally guess connectivity."
    )
    parser.add_argument("dlg_file", type=str, help="Path to the input DLG file.")
    parser.add_argument("--top_n", type=int, default=1, help="Number of top poses to output.")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory for output files.")
    parser.add_argument("--include_flexible", action="store_true", help="Include flexible residues in the output.")
    parser.add_argument("--log_file", type=str, help="Optional summary log file.")
    parser.add_argument("--receptor", type=str, help="Path to a rigid receptor pdbqt file to merge with dlg coordinates.")
    parser.add_argument("--guess_connectivity", action="store_true", help="If set, guess connectivity using MDAnalysis to write final PDB with CONECT records.")

    args = parser.parse_args()

    dlg_parser = DLGParser(args.dlg_file)
    dlg_parser.parse_ranking_table()
    dlg_parser.write_top_poses(
        output_dir=args.output_dir,
        top_n=args.top_n,
        include_flexible=args.include_flexible,
        receptor_file=args.receptor
    )
    dlg_parser.print_summary(log_file=args.log_file)

    # If connectivity guessing is requested and a receptor was provided (thus pdbqt output), process each file.
    if args.guess_connectivity and args.receptor:
        base_name = os.path.splitext(os.path.basename(args.dlg_file))[0]
        for i in range(1, args.top_n + 1):
            input_file = os.path.join(args.output_dir, f"{base_name}_pose_{i}_complex.pdbqt")
            output_file = os.path.join(args.output_dir, f"{base_name}_pose_{i}_final.pdb")
            pdbqt_to_pdb_with_conect(input_file, output_file)

if __name__ == "__main__":
    main()
