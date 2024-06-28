# invitroSPI
Spliced peptide identification from in vitro digestions of polypeptides with purified proteasomes  
<img src="invitroSPI_white.png" width="400">

Please cite the following publication if you are using invitroSPI for your research:

> Roetschke, H.P., Rodriguez-Hernandez, G., Cormican, J.A. et al. InvitroSPI and a large database of proteasome-generated spliced and non-spliced peptides. Sci Data 10, 18 (2023). https://doi.org/10.1038/s41597-022-01890-6

## overview
The invitroSPI pipeline consists of six main steps that are implemented in a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow:
1. parsing of search result files and creation of a preliminary MS database (*MSDB*) containing all peptide-spectrum matches (PSMs), mapping of peptides to their substrate origins and potential product type re-assignment
2. scoring and filtering of PSMs using Mascot's ion score and q-value as well as the delta score described in the manuscript
3. identification and, optionally, removal of synthesis errors using the control runs
4. mapping of peptides to the substrate sequence accounting for potential multi-mappers
5. database statistics and diagnostic information
6. MS2 spectra plotting with annotation of peptide fragment ions.

Additionally, code for the computation of all possible spliced and non-spliced peptides which is used for the Mascot search is provided in `SOURCE/_generateSearchDB_polypeptides.R`. We also include a script containing useful functions for downstream analyses (`invitroSPI_utils.R`). 

## execution
invitroSPI relies on [Conda](https://docs.conda.io/en/latest/) and Snakemake.
In order to install Conda, click on this [link](https://docs.conda.io/en/latest/miniconda.html) and follow the installation guidelines for your respective operating system.  
After installing Conda, you need to install Snakemake. The Snakemake installation procedure is described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Briefly, open the terminal on your computer and paste the following lines sequentially:  
`conda install -n base -c conda-forge mamba`  
`conda activate base`  
`mamba create -c conda-forge -c bioconda -n snakemake snakemake`  
Additionally, you might need to run `conda update conda`.

Download this repository as a .zip file (click on *Code* at the upper right corner of this repository --> Download ZIP), move the .zip file in the desired directory on your computer and unpack it.
**Make sure to edit `INPUT/sample_list.csv` and `INPUT/config.yaml` before starting the data processing** (see below for a more detailed description)
Open the terminal in this directory and enter: `conda activate snakemake`.

The pipeline can be executed by pasting `snakemake --use-conda --cores all -R createMSDB` into the terminal. The progress of the pipeline execution should appear in your terminal window.
In case you have installed an older version of Conda/Snakemake and encounter an error when executing the pipeline, try executing
`snakemake --use-conda --cores all -R createMSDB --conda-frontend conda`.

After your jobs finished, enter `conda deactivate` in order to terminate your Conda environment.

## input
invitroSPI identifies spliced and non-spliced peptides from Mascot search result files. Therefore, the user must provide a table `sample_list.csv` in the `INPUT/` folder containing information about:
- project name
- substrate ID
- substrate sequence
- time point
- search result file name
- replicate
- MSfile (optional)
- metainformation (optional)
- folder with .mgf files 

!!!
In case you are creating/editing the `sample_list.csv` file in Microsoft Excel, make sure to save it as actual comma-separated file. I.e., `Save as...`
 --> `Comma-separated Values (.csv)`
To check that the sample list is in the correct format, you can open it using a text editor and verify that the columns are separated by commas (and NOT semicolons).
!!!

Additionally, the user must provide **search result** files deposited in the folder `INPUT/search_results/` and list their names in the `sample_list.csv` table.

An example of the `sample_list.csv` table is given below and can be modified accordingly by the user:

| project_name | substrateID | substrateSeq | digestTime | filename | replicate | MSfile | metainformation |
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| test_data | TSN2 | TSN2.fasta | 4 | F029125.csv | 1 | | INPUT/metainformation_testData.csv |
| test_data | TSN89 |	RTKAWNRQLYPEW	| 4	| F029129.csv |	1 | | INPUT/metainformation_testData.csv |
| test_data | TSN2 | VSRQLRTKAWNRQLYPEWTEAQR |	CTRL |	F029123.csv |	1 | | INPUT/metainformation_testData.csv |
| test_data | TSN89 |	RTKAWNRQLYPEW |	CTRL |	F029127.csv |	1 | | INPUT/metainformation_testData.csv |

You can either paste the substrate sequence in the `sample_list` directly or put the name of a single-entry .fasta file containing the substrate sequence. This file should be located in `INPUT/sequences`.

Note that also the filenames of the control files must be provided. **For the control files, put `CTRL` in the `digestTime` column**. Please provide all other time points in hours (*e.g.*, 0.5 for 30 min).  
`MSfile` is an optional column that might be helpful to keep track of the .raw/.mgf files that were used to generate the respective search result files. It will not be used to create the database.  
In case you would like to include additional information in the final database (for instance species, location, instrument, substrate origin, proteasome isotype, ...), you could do so by creating a .csv table containing all this information. The .csv table must also contain a column `filename` corresponding to the `filename` in the sample list (see `INPUT/metainformation_SciData2022` for an example). Provide the name of the metainformation table in the `metainformation` column of the sample list.

The invitroSPI workflow is constructed in such a way that the user can keep appending projects and search result files to the `sample_list`. However, only the samples of the current project (which is specified in `INPUT/config.yaml`) are processed and stored separately in `OUTPUT/project_name` subdirectories.

## output
The pipeline provides a final database of identified spliced and non-spliced peptides in .csv format (`OUTPUT/project_name/ProteasomeDB.csv`) , MS2 spectra plots in .pdf format(`OUTPUT/project_name/ms2_spectra.pdf`) as well as all intermediate files in binary format(`OUTPUT/project_name/tmp/`).

Additionally, some diagnostic plots , database statistics are being produced which can also be found in the `OUTPUT/` folder. They comprise:
- number of unique peptides
- number and frequencies of peptides and PSMs at each time point
- length distribution of peptides, splice reactants and intervening sequences at each time point


## parameters that can be modified by the user
We are using the following default parameters that are specified in the `INPUT/config.yaml` file and that can be changed by the user:
- `project_name`: indicates which files listed in the `sample_list.csv` should be processed (enter valid directory names only! - no spaces or German umlauts)
- `delta_score`: 0.3
- `ion_score`: 20
- `q_value`: 0.05
- `include_scanNum`: "yes"
- `keep_synErrors`: "no"
- `mgf_folder` : "path_to_your_mgf_folder"

Processing of scan numbers is only possible if .mgf files were created with **msconvert or Mascot Distiller**. In case you provide search results in another format (not recommended), please set `include_scanNum` to "no".
In case you would like to include synthesis errors (labelled as such) in the final *ProteasomeDB*, change the `keep_synErrors` flag accordingly.
