import os
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import re
import numpy as np
import multiprocessing
import shutil
import PyPDF2
import subprocess
import tempfile


project_name = snakemake.config['project_name']

PROTON = 1.007276466622
ELECTRON = 0.00054858
H = 1.007825035
C = 12.0
O = 15.99491463
N = 14.003074

N_TERMINUS = H
C_TERMINUS = O + H
CO = C + O
CHO = C + H + O
NH2 = N + (H * 2)
H2O = (H * 2) + O
NH3 = N + (H * 3)

NEUTRAL_LOSSES = [NH3,H2O,]
LOSS_NAMES = {'': '', 'NH3': '*', 'H2O': '°'}
LOSS_WEIGHTS = {'': 0, 'NH3': NH3, 'H2O': H2O}
RESIDUE_WEIGHTS = {'A': 71.037114,'R': 156.101111,'N': 114.042927,'D': 115.026943, 'C': 103.009185,'E': 129.042593,'Q': 128.058578,'G': 57.021464,'H': 137.058912,'I': 113.084064,'L': 113.084064,'K': 128.094963,'M': 131.040485,'F': 147.068414, 'P': 97.052764,'S': 87.032028,'T': 101.047679,'W': 186.079313,'Y': 163.06332,'V': 99.068414}


KNOWN_PTM_WEIGHTS = {
    'Deamidated (N)': 0.984016,
    'Deamidated (NQ)': 0.984016,
    'Deamidation (NQ)': 0.984016,
    'Deamidation (N)': 0.984016,
    'Deamidation (Q)': 0.984016,
    'Oxidation (M)': 15.994915,
    'Acetyl (N-term)': 42.010565,
    'Acetylation (N-term)': 42.010565,
    'Acetyl (Protein N-term)': 42.010565,
    'Phospho (Y)': 79.966331,
    'Phospho (ST)': 79.966331,
    'Phospho (STY)': 79.966331,
    'Phosphorylation (STY)': 79.966331,
    'Carbamidomethyl (C)': 57.021464,
    'Carbamidomethylation': 57.021464,
    'unknown': 0.0
}


ION_OFFSET = {
    'a': N_TERMINUS - CHO,
    'b': N_TERMINUS - H,
    'c': N_TERMINUS + NH2,
    'x': C_TERMINUS + CO - H,
    'y': C_TERMINUS + H,
    'z': C_TERMINUS - NH2,
}



# Function to compute potential molecular weights
def compute_potential_mws(sequence, modifications, reverse, ptm_id_weights):
    sequence_length = len(sequence)
    n_fragments = sequence_length - 1
    mzs = np.empty(n_fragments)

    if modifications and isinstance(modifications, str) and modifications != 'nan' and modifications != 'unknown':
        ptms_list = modifications.split(".")
        mods_list = [int(mod) for mod in ptms_list[1]]

        if reverse:
            ptm_start = int(ptms_list[2])
            ptm_end = int(ptms_list[0])
            mods_list = mods_list[::-1]
        else:
            ptm_start = int(ptms_list[0])
            ptm_end = int(ptms_list[2])
    else:
        mods_list = None
        ptm_start = 0
        ptm_end = 0

    if reverse:
        sequence = sequence[::-1]

    if ptm_start:
        tracking_mw = ptm_id_weights[ptm_start]
    else:
        tracking_mw = 0

    for idx in range(n_fragments):
        tracking_mw += RESIDUE_WEIGHTS[sequence[idx]]
        if mods_list is not None and mods_list[idx]:
            tracking_mw += ptm_id_weights[mods_list[idx]]
        mzs[idx] = tracking_mw

    tracking_mw += RESIDUE_WEIGHTS[sequence[n_fragments]]
    if mods_list is not None and mods_list[n_fragments]:
        tracking_mw += ptm_id_weights[mods_list[n_fragments]]
    if ptm_end:
        tracking_mw += ptm_id_weights[ptm_end]

    return mzs, tracking_mw

# Function to get ion masses
def get_ion_masses(sequence, ptm_id_weights, modifications=None):
    sub_seq_mass, total_residue_mass = compute_potential_mws(
        sequence=sequence,
        modifications=modifications,
        reverse=False,
        ptm_id_weights=ptm_id_weights
    )

    rev_sub_seq_mass, _ = compute_potential_mws(
        sequence=sequence,
        modifications=modifications,
        reverse=True,
        ptm_id_weights=ptm_id_weights
    )

    prosit_ions = {}
    prosit_ions['b'] = ION_OFFSET['b'] + sub_seq_mass
    prosit_ions['y'] = ION_OFFSET['y'] + rev_sub_seq_mass

    return prosit_ions, (total_residue_mass + N_TERMINUS + C_TERMINUS)

# Function to match m/z values
def match_mz(base_mass, frag_z, observed_mzs, loss=0.0):
    fragment_mz = (base_mass + (frag_z * PROTON) - loss) / frag_z
    matched_mz_ind = np.argmin(np.abs(observed_mzs - fragment_mz))
    return observed_mzs[matched_mz_ind] - fragment_mz, matched_mz_ind


def extract_spectrum_for_mgf_file(mgf_file_path, mgf_dict):
    spectra = {}

    scan_pattern = re.compile(r'scan=(\d+)')
    charge_pattern = re.compile(r'(\d+)\+')
    data_pattern = re.compile(r"(\d+\.\d+)\s+(\d+)")

    with open(mgf_file_path, 'r') as file:
        current_scan_number = None
        mz_values = []
        intensities = []
        charge = None
        is_reading_spectrum = False

        for line in file:
            if 'TITLE' in line or 'NativeID' in line:
                match = scan_pattern.search(line)
                if match:
                    if current_scan_number is not None:
                        # Store the previous scan's data
                        spectra[current_scan_number] = {
                            'mz_values': mz_values,
                            'intensities': intensities,
                            'charge': charge
                        }

                    # Start new scan data
                    current_scan_number = int(match.group(1))
                    if current_scan_number not in mgf_dict[os.path.basename(mgf_file_path)]:
                        continue
                    mz_values = []
                    intensities = []
                    charge = None
                    is_reading_spectrum = True

            if is_reading_spectrum:
                if 'CHARGE' in line:
                    charge_match = charge_pattern.search(line)
                    if charge_match:
                        charge = int(charge_match.group(1))

                data_match = data_pattern.match(line)
                if data_match:
                    mz_values.append(float(data_match.group(1)))
                    intensities.append(int(data_match.group(2)))

                if line.strip() == 'END IONS':
                    # Store the last scan's data
                    spectra[current_scan_number] = {
                        'mz_values': mz_values,
                        'intensities': intensities,
                        'charge': charge
                    }
                    current_scan_number = None
                    is_reading_spectrum = False

        # Store the final spectrum if the file doesn't end with 'END IONS'
        if current_scan_number is not None:
            spectra[current_scan_number] = {
                'mz_values': mz_values,
                'intensities': intensities,
                'charge': charge
            }

    return spectra


def plot_matched_ions(peptide, scan_number, mgf_file_path, mz_values, intensities, charge, ax2):
    ax1 = ax2.inset_axes([0, 0.75, 1, 0.25])  # Top 25% for peptide fragment pattern

    modifications = None
    ion_data, _ = get_ion_masses(peptide, KNOWN_PTM_WEIGHTS, modifications)

    proton_charge = PROTON

    matched_b_ions = [
        (index + 1, ion_mass, charge, loss_name,
         (ion_mass + (charge * proton_charge) - loss_weight) / charge,
         np.interp((ion_mass + (charge * proton_charge) - loss_weight) / charge, mz_values, intensities))
        for ion_type, ions in ion_data.items()
        for index, ion_mass in enumerate(ions)
        for charge in range(1, 3)
        for loss_name, loss_weight in LOSS_WEIGHTS.items()
        if abs(match_mz(ion_mass, charge, mz_values, loss=loss_weight)[0]) < 0.02 and ion_type == 'b'
    ]

    matched_y_ions = [
        (index + 1, ion_mass, charge, loss_name,
         (ion_mass + (charge * proton_charge) - loss_weight) / charge,
         np.interp((ion_mass + (charge * proton_charge) - loss_weight) / charge, mz_values, intensities))
        for ion_type, ions in ion_data.items()
        for index, ion_mass in enumerate(ions)
        for charge in range(1, 3)
        for loss_name, loss_weight in LOSS_WEIGHTS.items()
        if abs(match_mz(ion_mass, charge, mz_values, loss=loss_weight)[0]) < 0.02 and ion_type == 'y'
    ]

    mgf_file_name = os.path.basename(mgf_file_path)
    title = f'Peptide: {peptide}, Scan Number: {scan_number}, File: {mgf_file_name}, Charge: {charge}'
    ax1.set_title(title, fontweight='bold', fontsize=12)

    for i, aa in enumerate(peptide):
        ax1.text(i + 1, 1, aa, ha='center', va='center', color='black', fontsize=10)

        if any(b_ion[0] == i + 1 for b_ion in matched_b_ions):
            ax1.plot([i + 1.5, i + 1.5], [0.8, 1], color='blue', linewidth=2)
            ax1.plot([i + 1.3, i + 1.5], [0.8, 0.8], color='blue', linewidth=2)
            ax1.text(i + 1.5, 0.75, f'b{i+1}', ha='center', va='top', color='blue', fontsize=8)

        if any(y_ion[0] == len(peptide) - i for y_ion in matched_y_ions):
            position_y = i + 0.5
            ax1.plot([position_y, position_y], [1, 1.2], color='red', linewidth=2)
            ax1.plot([position_y, position_y + 0.2], [1.2, 1.2], color='red', linewidth=2)
            ax1.text(position_y + 0.2, 1.25, f'y{len(peptide) - i}', ha='center', va='bottom', color='red', fontsize=8)

    ax1.set_xlim(0, len(peptide) + 1)
    ax1.set_ylim(0.7, 1.3)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.tick_params(axis='both', which='both', length=0)
    ax1.spines['top'].set_color('none')
    ax1.spines['right'].set_color('none')
    ax1.spines['left'].set_color('none')
    ax1.spines['bottom'].set_color('none')

    ax2.stem(mz_values, intensities, linefmt="k-", markerfmt=" ", basefmt="k-")

    for b_ion in matched_b_ions:
        ax2.stem([b_ion[4]], [b_ion[5]], linefmt="b-", markerfmt=" ", basefmt="b-")
        ax2.annotate(f'b{b_ion[0]}{LOSS_NAMES[b_ion[3]]}{"" if b_ion[2] == 1 else "+" if b_ion[2] == 2 else f"{b_ion[2]}^"}',
                     (b_ion[4], b_ion[5]),
                     textcoords="offset points",
                     xytext=(0, 2 + b_ion[0] * 1),
                     ha='center', color='blue', fontsize=8)

    for y_ion in matched_y_ions:
        ax2.stem([y_ion[4]], [y_ion[5]], linefmt="r-", markerfmt=" ", basefmt="r-")
        ax2.annotate(f'y{y_ion[0]}{LOSS_NAMES[y_ion[3]]}{"" if y_ion[2] == 1 else "+" if y_ion[2] == 2 else f"{y_ion[2]}^"}',
                     (y_ion[4], y_ion[5]),
                     textcoords="offset points",
                     xytext=(0, 0 + y_ion[0] * 1),
                     ha='center', color='red', fontsize=8)

    ax2.set_xlabel('m/z', fontsize=12)
    ax2.set_ylabel('Intensity', fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=10)


def process_spectra(mgf_file_name, mgf_folder, proteasome_db):
    mgf_file_path = os.path.join(mgf_folder, mgf_file_name)
    mgf_scan_numbers = set()
    with open(mgf_file_path, 'r') as mgf_file:
        for line in mgf_file:
            match = re.search(r'scan=(\d+)', line)
            if match:
                scan_number = int(match.group(1))
                mgf_scan_numbers.add(scan_number)

    matching_entries = proteasome_db[proteasome_db['MSfile'] == mgf_file_name]
    matching_entries = matching_entries[matching_entries['scanNum'].isin(mgf_scan_numbers)]

    if not matching_entries.empty:
        print(f"Processing {mgf_file_name}...")
        grouped = matching_entries.groupby('pepSeq')['scanNum'].apply(list)
        return [(peptide, scan_numbers, mgf_file_path) for peptide, scan_numbers in grouped.items()]
    else:
        print(f"No matching entry found for scans in {mgf_file_name}. Skipping...")
        return []


def convert_spectra_data_to_dict(spectra_data):
    spectra_dict = {}
    mgf_dict ={}
    for peptide, scan_numbers, mgf_file_path in spectra_data:
        if peptide in spectra_dict:
            spectra_dict[peptide].append((scan_numbers, mgf_file_path))
        else:
            spectra_dict[peptide] = [(scan_numbers, mgf_file_path)]
        mgf_filename=os.path.basename(mgf_file_path)
        if mgf_filename in mgf_dict:
            mgf_dict[mgf_filename].update(scan_numbers)
        else:
            mgf_dict[mgf_filename]=set(scan_numbers)

    return spectra_dict, mgf_dict

def plot_single_peptide(peptide, spectra_dict, mgf_data, pdf_filename, temp_folder):
    pdf_path = os.path.join(temp_folder, "temp_" + peptide + "_" + pdf_filename)
    with PdfPages(pdf_path) as pdf:
        for scan_numbers, mgf_file_path in spectra_dict[peptide]:
            num_plots = len(scan_numbers)
            num_pages = (num_plots + 1) // 2
            for page in range(num_pages):
                fig, axs = plt.subplots(2, 1, figsize=(16, 10))
                fig.subplots_adjust(hspace=0.5)

                for i in range(2):
                    plot_index = page * 2 + i
                    if plot_index >= num_plots:
                        axs[i].axis('off')
                        continue

                    scan_number = scan_numbers[plot_index]
                    mz_values = mgf_data[os.path.basename(mgf_file_path)][scan_number]['mz_values']
                    intensities = mgf_data[os.path.basename(mgf_file_path)][scan_number]['intensities']
                    charge = mgf_data[os.path.basename(mgf_file_path)][scan_number]['charge']

                    if mz_values and intensities:
                        plot_matched_ions(peptide, scan_number, mgf_file_path, mz_values, intensities, charge, axs[i])
                    else:
                        axs[i].axis('off')

                pdf.savefig(fig)
                plt.close(fig)
    print(f"Spectra plots saved as '{pdf_path}'")


def combine_pdfs_into_single_pdf(input_folder, output_file):
    # Initialize PDF merger
    pdf_merger = PyPDF2.PdfMerger()

    # Get list of PDF files starting with 'temp_' in input_folder
    pdf_files = [file for file in os.listdir(input_folder) if file.startswith('temp_') and file.lower().endswith('.pdf')]

    # Sort PDF files to ensure they are combined in order
    pdf_files.sort()

    # Add each PDF file to the merger
    for pdf_file in pdf_files:
        pdf_merger.append(os.path.join(input_folder, pdf_file))

    # Write merged PDF to output_file
    with open(output_file, 'wb') as merged_pdf_file:
        pdf_merger.write(merged_pdf_file)

    print(f"Combined {len(pdf_files)} PDFs into {output_file}")




def process_spectra_and_plot(proteasome_db, mgf_folder, output_file_path, num_processes=15):
    mgf_files = [f for f in os.listdir(mgf_folder) if f.endswith(".mgf")]
    pdf_filename = "spectra_plots.pdf"

    spectra_data = []

    for mgf_file_name in mgf_files:
        new_data=process_spectra(mgf_file_name, mgf_folder, proteasome_db)
        spectra_data.extend(new_data)

    spectra_dict, mgf_dict=convert_spectra_data_to_dict(spectra_data)

    mgf_data = {}
    for mgf_file_name in mgf_files:
        mgf_file_path = os.path.join(mgf_folder, mgf_file_name)
        if mgf_file_name in mgf_dict.keys():
            mgf_data[os.path.basename(mgf_file_name)]=extract_spectrum_for_mgf_file(mgf_file_path, mgf_dict)
        
    # if not os.path.exists(temp_folder):
    #     os.makedirs(temp_folder)
    # else:
    #     shutil.rmtree(temp_folder)
    #     os.makedirs(temp_folder)

    # Assuming spectra_dict.keys() contains the list of peptides
    peptides = list(spectra_dict.keys())

    # Number of processes to use (adjust according to your system capabilities)
    # num_processes = multiprocessing.cpu_count()  # Use number of available CPU cores
    # num_processes=15
    with tempfile.TemporaryDirectory() as temp_folder:
        # Initialize a pool of processes
        with multiprocessing.Pool(processes=num_processes) as pool:
            # Map the plotting function to each peptide using the pool
            pool.starmap(plot_single_peptide, [(peptide, spectra_dict, mgf_data, pdf_filename, temp_folder) for peptide in peptides])
        combine_pdfs_into_single_pdf(temp_folder, output_file_path)

    # print(f"Spectra plots saved as '{pdf_filename}'")

def process_proteasome_db(proteasome_db_path, sample_list_db_path):
  proteasome_df = pd.read_csv(proteasome_db_path)
  proteasome_df['runID'] = proteasome_df['runID'].apply(lambda x: x[x.find('F'):])

  sample_df = pd.read_csv(sample_list_db_path)


  filename_to_MSfile = dict(zip(sample_df['filename'].str.split('.').str[0], sample_df['MSfile']))

  # Extract filename without extension
  proteasome_df['filename_base'] = proteasome_df['runID'].str.split('.').str[0]

  # Add a new 'MSfile' column by matching filename (without extension) with the corresponding MSfile
  proteasome_df['MSfile'] = proteasome_df['filename_base'].map(filename_to_MSfile).fillna('Not Found')

  # Drop the intermediate column
  proteasome_df.drop(columns=['filename_base'], inplace=True)

  # Return the modified DataFrame
  return proteasome_df



proteasome_db_path=snakemake.input.ProteasomeDB
sample_list_db_path=snakemake.input.sample_list
mgf_folder=snakemake.params.mgf_folder

#output_folder=snakemake.output.output_folder
#output_filename=snakemake.output.output_filename
spectra_plot=snakemake.output.spectra_plot
num_processes = 15 #[TODO:] take this from snakemake

proteasome_db=process_proteasome_db(proteasome_db_path, sample_list_db_path)

#temp_folder = os.path.join(output_folder, 'temp')
#temp_folder = f'OUTPUT/{project_name}/temp'

output_file_path = spectra_plot

process_spectra_and_plot(proteasome_db, mgf_folder, output_file_path, num_processes)



# combine_pdfs_into_single_pdf(temp_folder, output_file_path)

