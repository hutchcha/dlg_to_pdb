import os
import re

class DLGParser:
    def __init__(self, dlg_file_path):
        """
        Initialize the DLGParser with the path to a DLG file.

        Parameters:
            dlg_file_path (str): Path to the AutoDock DLG file.
        """
        self.dlg_file_path = dlg_file_path
        self.poses = []  # List to store pose data (rank, run, binding energy, etc.)
        self.top_poses = []

    def parse_ranking_table(self):
        """
        Parse the ranking table near the end of the DLG file to extract pose rankings.
        """
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
        except FileNotFoundError:
            print(f"Error: File not found - {self.dlg_file_path}")
        except Exception as e:
            print(f"An error occurred while parsing the DLG file: {e}")

    def get_top_scoring_pose(self, top_n=1):
        """
        Get the top scoring poses based on the ranking.

        Parameters:
            top_n (int): Number of top poses to retrieve.
        """
        try:
            self.top_poses = sorted(self.poses, key=lambda x: (x["rank"], x["sub_rank"]))[:top_n]
            if not self.top_poses:
                print("No poses found. Ensure the DLG file was parsed correctly.")
            else:
                print(f"Retrieved top {len(self.top_poses)} scoring poses.")
        except Exception as e:
            print(f"An error occurred while retrieving top scoring poses: {e}")

    def extract_pose_coordinates(self, run_number, include_flexible=False):
        """
        Extract ligand atom coordinates (and optionally flexible residues) for a specific run.

        Parameters:
            run_number (int): The run number to extract coordinates for.
            include_flexible (bool): Whether to include flexible residue atoms.

        Returns:
            list: A list of strings representing the atom coordinates in PDB format.
        """
        pose_coordinates = []
        recording = False

        try:
            with open(self.dlg_file_path, "r") as file:
                for line in file:
                    if recording:
                        if line.startswith("DOCKED: HETATM"):
                            pdb_line = self.format_pdb_line(line[8:].strip())
                            pose_coordinates.append(pdb_line)
                        elif include_flexible and line.startswith("DOCKED: ATOM"):
                            pdb_line = self.format_pdb_line(line[8:].strip(), flexible=True)
                            pose_coordinates.append(pdb_line)
                        elif line.startswith("DOCKED: ENDMDL"):
                            print(f"Finished reading coordinates for run {run_number}.")
                            break
                    elif re.match(f"DOCKED:\s*MODEL\s*{run_number}\s*", line):
                        recording = True
                        print(f"Started reading coordinates for run {run_number}.")

            if not pose_coordinates:
                print(f"No coordinates found for run {run_number}.")
        except Exception as e:
            print(f"An error occurred while extracting pose coordinates: {e}")

        return pose_coordinates

    @staticmethod
    def format_pdb_line(docked_line, flexible=False):
        """
        Format a DOCKED: HETATM or ATOM line into standard PDB format.

        Parameters:
            docked_line (str): Line from the DLG file starting with HETATM or ATOM.
            flexible (bool): Whether the line represents a flexible residue.

        Returns:
            str: Reformatted line in standard PDB format.
        """
        parts = docked_line.split()
        atom_name = parts[2]
        residue_name = parts[3]
        chain_id = parts[4] if flexible else "A"  # Use chain ID from file if flexible
        residue_id = int(parts[5]) if flexible else int(parts[4])
        x, y, z = map(float, parts[6:9]) if flexible else map(float, parts[5:8])
        element = parts[-1]

        return f"HETATM{int(parts[1]):>5} {atom_name:<4} {residue_name:<3} {chain_id}{residue_id:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}           {element:>2}"

    def write_pdb(self, pose, output_dir=".", pose_number=1, include_flexible=False):
        """
        Write a PDB file for a given pose.

        Parameters:
            pose (dict): A pose dictionary containing rank, sub_rank, run, etc.
            output_dir (str): Directory to save the PDB file.
            pose_number (int): Pose number for unique file naming.
            include_flexible (bool): Whether to include flexible residue atoms.
        """
        try:
            coordinates = self.extract_pose_coordinates(pose["run"], include_flexible=include_flexible)
            if not coordinates:
                print(f"No coordinates to write for pose with rank {pose['rank']}.")
                return

            base_name = os.path.splitext(os.path.basename(self.dlg_file_path))[0]
            pdb_filename = os.path.join(output_dir, f"{base_name}_pose_{pose_number}.pdb")
            with open(pdb_filename, "w") as pdb_file:
                pdb_file.write("\n".join(coordinates) + "\n")

            print(f"PDB file written: {pdb_filename}")
        except Exception as e:
            print(f"An error occurred while writing PDB file: {e}")

    def write_top_poses(self, output_dir=".", top_n=1, include_flexible=False):
        """
        Write PDB files for the top scoring poses.

        Parameters:
            output_dir (str): Directory to save the PDB files.
            top_n (int): Number of top poses to write.
            include_flexible (bool): Whether to include flexible residue atoms.
        """
        try:
            self.get_top_scoring_pose(top_n=top_n)
            for idx, pose in enumerate(self.top_poses, start=1):
                self.write_pdb(pose, output_dir=output_dir, pose_number=idx, include_flexible=include_flexible)
        except Exception as e:
            print(f"An error occurred while writing top PDB poses: {e}")

    def print_summary(self, log_file=None):
        """
        Print a summary of the top poses and optionally write to a log file.

        Parameters:
            log_file (str): Path to the log file (optional).
        """
        summary = ""
        if not self.top_poses:
            summary = "No top poses available to print."
        else:
            summary_lines = ["Top Scoring Poses:"]
            for pose in self.top_poses:
                summary_lines.append(f"Rank: {pose['rank']}, Sub-Rank: {pose['sub_rank']}, Run: {pose['run']}, Binding Energy: {pose['binding_energy']} kcal/mol")
            summary = "\n".join(summary_lines)

        print(summary)
        if log_file:
            try:
                with open(log_file, "w") as file:
                    file.write(summary + "\n")
                print(f"Summary written to {log_file}")
            except Exception as e:
                print(f"An error occurred while writing the log file: {e}")

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Process AutoDock DLG files and output PDB files.")
    parser.add_argument("dlg_file", type=str, help="Path to the input DLG file.")
    parser.add_argument("--top_n", type=int, default=1, help="Number of top poses to output.")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to save the PDB files.")
    parser.add_argument("--include_flexible", action="store_true", help="Include flexible residues in the output PDB files.")
    parser.add_argument("--log_file", type=str, help="Optional log file to save the summary.")

    args = parser.parse_args()

    dlg_parser = DLGParser(args.dlg_file)
    dlg_parser.parse_ranking_table()
    dlg_parser.write_top_poses(output_dir=args.output_dir, top_n=args.top_n, include_flexible=args.include_flexible)
    dlg_parser.print_summary(log_file=args.log_file)

if __name__ == "__main__":
    main()
