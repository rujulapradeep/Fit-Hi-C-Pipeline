import os
import subprocess
import gzip
import shutil
from glob import glob
import argparse

#function to create contact files using juicertools
def create_contact_files(hic_file, chromosome1, chromosome2, juicer_jar_path, output_dir, log_file):
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{chromosome1}_{chromosome2}_counts.txt')
    
    cmd = f"java -jar {juicer_jar_path} dump observed VC {hic_file} {chromosome1} {chromosome2} BP 500000 -p {output_file}"
    with open(log_file, 'a') as log:
        try:
            subprocess.run(cmd, shell=True, check=True, stdout=log, stderr=log)
        except subprocess.CalledProcessError as e:
            print(f"Error running Juicer Tools for {chromosome1}-{chromosome2}: {e}")

#function to modify contact file to correct format for fithic
def modify_contact_file(hic_file, chromosome1, chromosome2, output_dir):
    input_file = os.path.join(output_dir, f'{chromosome1}_{chromosome2}_counts.txt')
    modified_output_dir = os.path.join(output_dir, 'modified')
    os.makedirs(modified_output_dir, exist_ok=True)
    output_file = os.path.join(modified_output_dir, f'{chromosome1}_{chromosome2}_counts.txt')
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                fields = line.strip().split()
                outfile.write(f'{chromosome1} {fields[0]} {chromosome2} {fields[1]} {fields[2]}\n')
    except Exception as e:
        print(f"Error modifying contact file for {chromosome1}-{chromosome2}: {e}")

#function to create fragments file for fithic
def create_fithic_fragments(chromsizes_file, output_file, fithic_path, resolution, log_file):
    cmd = f"python {fithic_path}/fithic/utils/createFitHiCFragments-fixedsize.py --chrLens {chromsizes_file} --outFile {output_file} --resolution {resolution}"
    with open(log_file, 'a') as log:
        try:
            subprocess.run(cmd, shell=True, check=True, stdout=log, stderr=log)
        except subprocess.CalledProcessError as e:
            print(f"Error creating Fit-Hi-C fragments: {e}")

#function to run fithic
def run_fithic(contact_file, fragments_file, output_dir, resolution, fithic_path, log_file):
    contact_type = 'interOnly'
    parts = os.path.basename(contact_file).split('_')
    if len(parts) == 3 and parts[0] == parts[1]:
        contact_type = 'intraOnly'
    
    zipped_file = f"{contact_file}.gz"
    try:
        with open(contact_file, 'rb') as f_in:
            with gzip.open(zipped_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        cmd = [
            'python', f"{fithic_path}/fithic/fithic.py",
            '--interactions', zipped_file,
            '--fragments', fragments_file,
            '--outdir', output_dir,
            '--resolution', str(resolution),
            '--contactType', contact_type
        ]
        with open(log_file, 'a') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        os.remove(zipped_file)
    except subprocess.CalledProcessError as e:
        print(f"Error running Fit-Hi-C for {contact_file}: {e}")
        raise
    except Exception as e:
        print(f"Error running Fit-Hi-C for {contact_file}: {e}")
        raise

#function to combine all inter and intra significant interaction into one file 
def combine_significant_interactions(input_dirs, output_file):
    combined_significant_interactions = []

    for contact_dir in input_dirs:
        subdirectories = [d for d in os.listdir(contact_dir) if os.path.isdir(os.path.join(contact_dir, d))]
        for subdir in subdirectories:
            subdir_path = os.path.join(contact_dir, subdir)
            files_in_dir = os.listdir(subdir_path)
            for file_name in files_in_dir:
                if file_name.endswith('.txt.gz'):
                    significant_file = os.path.join(subdir_path, file_name)
                    try:
                        with gzip.open(significant_file, 'rt') as f:
                            for line in f:
                                combined_significant_interactions.append(line)
                    except Exception as e:
                        print(f"Error reading {significant_file}: {e}")

    with open(output_file, 'w') as out_file:
        out_file.writelines(combined_significant_interactions)
        
def main(hic_samples_file, output_dir, chromsizes_file, juicer_jar_path, fithic_path, resolution):
    log_file = os.path.join(output_dir, 'log.txt')
    os.makedirs(output_dir, exist_ok=True)

    with open(hic_samples_file, 'r') as f:
        hic_samples = [line.strip().split() for line in f.readlines()]
    
    chromosomes = [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY'] #change depending on genome
    
    for hic_file, sample_name in hic_samples:
        print(f"Creating contact count files for {sample_name}...")
        sample_output_dir = os.path.join(output_dir, sample_name)
        contact_counts_dir = os.path.join(sample_output_dir, 'contactcounts')
        chr_sig_interactions_dir = os.path.join(sample_output_dir, 'chr_sig_interactions')
        
        os.makedirs(contact_counts_dir, exist_ok=True)
        os.makedirs(chr_sig_interactions_dir, exist_ok=True)
        
        for i in range(len(chromosomes)):
            for j in range(i, len(chromosomes)):
                create_contact_files(hic_file, chromosomes[i], chromosomes[j], juicer_jar_path, contact_counts_dir, log_file)
                modify_contact_file(hic_file, chromosomes[i], chromosomes[j], contact_counts_dir)
    
    fragments_file = 'fragments.txt'
    print("Creating Fit-Hi-C fragments file...")
    create_fithic_fragments(chromsizes_file, fragments_file, fithic_path, resolution, log_file)
    
    for hic_file, sample_name in hic_samples:
        print(f"Running Fit-Hi-C for {sample_name}...")
        contact_counts_dir = os.path.join(output_dir, sample_name, 'contactcounts', 'modified')
        modified_contact_files = glob(os.path.join(contact_counts_dir, '*.txt'))
        for contact_file in modified_contact_files:
            with open(contact_file, 'r') as f:
                lines = f.readlines()
            filtered_lines = [line for line in lines if 'Infinity' not in line.split()[4]]
            with open(contact_file, 'w') as f:
                f.writelines(filtered_lines)
            
            contact_file_name = os.path.basename(contact_file)
            contact_output_dir = os.path.join(output_dir, sample_name, 'chr_sig_interactions', contact_file_name)
            os.makedirs(contact_output_dir, exist_ok=True)
            
            try:
                run_fithic(contact_file, fragments_file, contact_output_dir, resolution, fithic_path, log_file)
            except Exception:
                print(f"Error running Fit-Hi-C for {sample_name} for {contact_file_name}. Stopping execution.")
                return
    
    significant_interactions_dir = os.path.join(output_dir, 'Significant_Interactions')
    os.makedirs(significant_interactions_dir, exist_ok=True)

    for hic_file, sample_name in hic_samples:
        print(f"Combining significant interactions for {sample_name}...")
        chr_sig_interactions_dir = os.path.join(output_dir, sample_name, 'chr_sig_interactions')

        significant_output_file = os.path.join(significant_interactions_dir, sample_name, 'significant_interactions.txt')
        os.makedirs(os.path.dirname(significant_output_file), exist_ok=True)

        combine_significant_interactions([chr_sig_interactions_dir], significant_output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the Fit-Hi-C pipeline for Hi-C data.")
    parser.add_argument("hic_samples_file", help="Path to the text file containing Hi-C sample paths and names")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument("chromsizes_file", help="Path to the chromosome sizes file")
    parser.add_argument("juicer_jar_path", help="Path to Juicer Tools .jar file")
    parser.add_argument("fithic_path", help="Path to the Fit-Hi-C directory")
    parser.add_argument("resolution", type=int, help="Resolution for the Fit-Hi-C analysis")

    args = parser.parse_args()
    
    main(args.hic_samples_file, args.output_dir, args.chromsizes_file, args.juicer_jar_path, args.fithic_path, args.resolution)
