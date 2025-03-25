import os
import subprocess
import sys
from pathlib import Path

def run_autodock_vina(receptor_path, ligand_path):
    """
    Run AutoDock Vina with given receptor and ligand files
    
    Args:
        receptor_path (str): Path to receptor PDBQT file
        ligand_path (str): Path to ligand PDBQT file
    """
    # Create results directory if it doesn't exist
    Path("results").mkdir(exist_ok=True)
    
    # Generate output filename based on ligand name
    ligand_name = Path(ligand_path).stem
    output_file = f"results/output_{ligand_name}.pdbqt"
    
    # Default configuration parameters
    config = {
        "center_x": 0,
        "center_y": 0,
        "center_z": 0,
        "size_x": 20,
        "size_y": 20,
        "size_z": 20,
        "exhaustiveness": 8,
        "num_modes": 9,
        "energy_range": 3
    }
    
    # Construct Vina command
    command = [
        "vina",
        "--receptor", receptor_path,
        "--ligand", ligand_path,
        "--center_x", str(config["center_x"]),
        "--center_y", str(config["center_y"]),
        "--center_z", str(config["center_z"]),
        "--size_x", str(config["size_x"]),
        "--size_y", str(config["size_y"]),
        "--size_z", str(config["size_z"]),
        "--exhaustiveness", str(config["exhaustiveness"]),
        "--num_modes", str(config["num_modes"]),
        "--energy_range", str(config["energy_range"]),
        "--out", output_file
    ]
    
    try:
        print(f"Starting AutoDock Vina...")
        print(f"Receptor: {receptor_path}")
        print(f"Ligand: {ligand_path}")
        print(f"Output will be saved to: {output_file}")
        
        # Run Vina command
        result = subprocess.run(command, capture_output=True, text=True)
        
        # Print output
        print("\nAutoDock Vina Output:")
        print(result.stdout)
        
        if result.returncode == 0:
            print(f"\nDocking completed successfully!")
            print(f"Results saved to: {output_file}")
        else:
            print("\nError running AutoDock Vina:")
            print(result.stderr)
            
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

def main():
    # Check if correct number of arguments provided
    if len(sys.argv) != 3:
        print("Usage: python3 dock.py <receptor.pdbqt> <ligand.pdbqt>")
        print("Example: python3 dock.py receptor/pocket.pdbqt ligand/ligand-b.pdbqt")
        sys.exit(1)
    
    # Get file paths from command line arguments
    receptor_path = sys.argv[1]
    ligand_path = sys.argv[2]
    
    # Check if files exist
    if not os.path.exists(receptor_path):
        print(f"Error: Receptor file '{receptor_path}' not found")
        sys.exit(1)
    
    if not os.path.exists(ligand_path):
        print(f"Error: Ligand file '{ligand_path}' not found")
        sys.exit(1)
    
    # Run docking
    run_autodock_vina(receptor_path, ligand_path)

if __name__ == "__main__":
    main()

