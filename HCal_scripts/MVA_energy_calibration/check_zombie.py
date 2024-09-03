import os
import ROOT
import shutil


def is_zombie(file_path):
    """Check if a ROOT file is a zombie."""
    try:
        root_file = ROOT.TFile.Open(file_path)
        if not root_file or root_file.IsZombie():
            return True
        root_file.Close()
    except OSError as e:
        print(f"Error opening file {file_path}: {e}")
        return True
    except Exception as e:
        print(f"Unexpected error with file {file_path}: {e}")
        return True
    return False

def move_zombie_files(directory, zombie_dir):
    """Find and move all zombie ROOT files to a subdirectory named 'zombies'."""
    if not os.path.exists(zombie_dir):
        os.makedirs(zombie_dir)
    
    zombie_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.root'):
                file_path = os.path.join(root, file)
                if is_zombie(file_path):
                    try:
                        shutil.move(file_path, os.path.join(zombie_dir, file))
                        zombie_files.append(file_path)
                        print(f"Moved zombie file: {file_path}")
                    except Exception as e:
                        print(f"Error moving file {file_path}: {e}")
    return zombie_files


def find_zombie_files(directory):
    """Find and list all zombie ROOT files in the given directory."""
    zombie_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.root'):
                file_path = os.path.join(root, file)
                if is_zombie(file_path):
                    zombie_files.append(file_path)
    return zombie_files

def main():
    #directory = '/eos/user/m/mmlynari/FCC_fellow/FCC_rootfile_storage/MVA_training_v28Aug24_FSR/240902_energies_1mil_SWandTopo_noNoise_ECalHCal_oldBenchmark_BRT_training'  # Replace with your directory
    
    directory = '/eos/user/m/mmlynari/FCC_fellow/FCC_rootfile_storage/MVA_training_v28Aug24_FSR/240902_energies_10kevt_cells_SW_noNoise_ECalHCal_oldBenchmark_BRT_validation'
    zombie_dir = os.path.join(directory, 'zombies')
    
    zombies = move_zombie_files(directory, zombie_dir)
    #zombies = find_zombie_files(directory)
    
    if zombies:
        print("Zombie files found:")
        for zombie in zombies:
            print(zombie)
    else:
        print("No zombie files found.")

if __name__ == "__main__":
    main()
