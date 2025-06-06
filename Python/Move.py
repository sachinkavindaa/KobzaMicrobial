import os
import shutil

# Set source and destination directories
source_dir = "/Users/echandrasekara2/Documents/Kobza/Kobza Data Files/all_fastq_gz"
destination_dir = "/Users/echandrasekara2/Documents/Kobza"

# Create destination directory if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)

# Loop from 57 to 84 and move matching files
for i in range(57, 85):  # includes 84
    filename = f"{i}_S{i}_L001_R1_001.fastq.gz"
    source_path = os.path.join(source_dir, filename)
    destination_path = os.path.join(destination_dir, filename)
    
    if os.path.exists(source_path):
        shutil.move(source_path, destination_path)
        print(f"Moved: {filename}")
    else:
        print(f"File not found: {filename}")
