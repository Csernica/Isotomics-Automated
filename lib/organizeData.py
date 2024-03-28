import os
import pandas as pd

def create_sample_standard_subfolders(parent_folder):
    '''
    helper function for create_nested_folders. Doesn't return anything, just makes subdirectories. 
    '''
    for subfolder in ["Smp", "Std"]:
        subfolder_path = os.path.join(parent_folder, subfolder)
        os.makedirs(subfolder_path, exist_ok=True)

def create_folder_with_subfolders(folder_path):
    """
    Attempts to create the folder and its subfolders, handling exceptions appropriately.
    """
    try:
        os.makedirs(folder_path, exist_ok=True)  # Adjusted for idempotency with exist_ok=True
        print(f"Folder created or already exists: {folder_path}")
        # Add "sample" and "standard" subfolders
        create_sample_standard_subfolders(folder_path)

    except Exception as e:
        print(f"Error creating folder {folder_path}: {e}")

def create_nested_folders(fragment_dict, parent_folder=".", oneFragFile=False):
    """
    Create nested folders based on a dictionary.

    Parameters:
    - fragment_dict: A dictionary where keys represent folder names.
    - parent_folder: The parent folder where the structure will be created. Default is the current directory.
    - oneFragFile: Indicates whether a single .isox file with data from all fragments is being provided. 
    """ 
    file_paths_dict = {}
    
    # Define special subfolders for specific fragment names
    special_subfolders = {
        'full': ['full_relative_abundance', 'full_molecular_average']
    }
    
    for fragment_name in fragment_dict.keys():
        # Handle special subfolders for the 'full' fragment
        if fragment_name in special_subfolders:
            for subfolder in special_subfolders[fragment_name]:
                special_folder_path = os.path.join(parent_folder, subfolder)
                create_folder_with_subfolders(special_folder_path)
                file_paths_dict.setdefault(fragment_name, []).append(special_folder_path)
            
        elif oneFragFile:
            # Use a common folder for all fragments except when handling 'full'
            folder_path = os.path.join(parent_folder, 'all_fragments')
            create_folder_with_subfolders(folder_path)
            file_paths_dict[fragment_name] = [folder_path]
        else:
            folder_path = os.path.join(parent_folder, fragment_name)
            create_folder_with_subfolders(folder_path)
            file_paths_dict[fragment_name] = [folder_path]

    return file_paths_dict, parent_folder

def get_file_paths_in_subfolders(folder_path, file_extensions='.isox'):
    """
    Recursively collects paths of files with specified extensions within subfolders of the specified folder.

    Parameters:
    - folder_path: A string representing the path to the main folder.
    - file_extensions: A list of strings representing the desired file extensions (e.g., ['.txt', '.csv']). 
                       If None, all files will be included.

    Returns:
    - A dictionary where keys are subfolder names and values are lists of file paths.
    """
    subfolder_files = []
    subfolder_file_order = []

    # Walk through the folder and its subfolders
    for root, dirs, files in os.walk(folder_path):
        # Iterate through subfolders
        for subfolder in dirs:
            subfolder_path = os.path.join(root, subfolder)


            # Collect paths of files within the subfolder
            for file in os.listdir(subfolder_path):
                file_path = os.path.join(subfolder_path, file)
                if os.path.isfile(file_path):
                    # Check if the file has the desired extension
                    if file_extensions is None or any(file.endswith(ext) for ext in file_extensions):
                        subfolder_files.append(file_path)
                        subfolder_file_order.append(os.path.basename(os.path.dirname(file_path)))

    return subfolder_files, subfolder_file_order

def get_subfolder_paths(folder_path):
    """
    Get a list of paths of subfolders within a folder.

    Parameters:
    - folder_path: A string representing the path to the main folder.

    Returns:
    - A list of strings representing the paths of subfolders within the main folder.
    """
    subfolder_paths = []

    # Get the list of entries (files and subfolders) within the main folder
    entries = os.listdir(folder_path)

    # Iterate through the entries
    for entry in entries:
        # Construct the full path of the entry
        entry_path = os.path.join(folder_path, entry)
        # Check if the entry is a directory (subfolder)
        if os.path.isdir(entry_path):
            # Add the subfolder path to the list
            subfolder_paths.append(entry_path)

    return subfolder_paths

