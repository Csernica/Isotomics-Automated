import os
import pandas as pd

def create_nested_folders(fragment_dict, parent_folder="."):
    """
    Create nested folders based on a dictionary.

    Parameters:
    - fragment_dict: A dictionary where keys represent folder names.
    - parent_folder: The parent folder where the structure will be created. Default is the current directory.
    """

    file_paths_dict = {}
    rtnParentFolder = parent_folder

    for fragment_name in fragment_dict:
        folder_path = os.path.join(parent_folder, fragment_name)

        # Add the folder path to the dictionary
        if fragment_name in file_paths_dict:
            file_paths_dict[fragment_name].append(folder_path)
        else:
            file_paths_dict[fragment_name] = [folder_path]

        try:
            if fragment_name == 'full':
                # Create separate folders for full_relative_abundance and full_molecular_average
                relative_abundance_folder = os.path.join(parent_folder, 'full_relative_abundance')
                molecular_average_folder = os.path.join(parent_folder, 'full_molecular_average')
                os.makedirs(relative_abundance_folder)
                print(f"Folder created: {relative_abundance_folder}")
                os.makedirs(molecular_average_folder)
                print(f"Folder created: {molecular_average_folder}")

                # Add "sample" and "standard" subfolders for the new folders
                sample_folder = os.path.join(relative_abundance_folder, "Smp")
                standard_folder = os.path.join(relative_abundance_folder, "Std")
                os.makedirs(sample_folder)
                os.makedirs(standard_folder)

                sample_folder = os.path.join(molecular_average_folder, "Smp")
                standard_folder = os.path.join(molecular_average_folder, "Std")
                os.makedirs(sample_folder)
                os.makedirs(standard_folder)
            else:
                os.makedirs(folder_path)
                print(f"Folder created: {folder_path}")
                # Add "sample" and "standard" subfolders
                sample_folder = os.path.join(folder_path, "Smp")
                standard_folder = os.path.join(folder_path, "Std")
                os.makedirs(sample_folder)
                os.makedirs(standard_folder)

        except FileExistsError:
            print(f"Folder already exists: {folder_path}")
        except Exception as e:
            print(f"Error creating folder {folder_path}: {e}")
            
    return file_paths_dict, rtnParentFolder


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

