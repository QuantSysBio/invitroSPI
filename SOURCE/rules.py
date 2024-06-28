wildcard_constraints:
	project_name = config['project_name']


rule createMSDB:
	input:
		sample_list = "INPUT/sample_list.csv"
	output:
		MSDB = "OUTPUT/{project_name}/tmp/MSDB.RData"
	params:
		project_name = config['project_name'],
		include_scanNum = config['include_scanNum']
	conda:
		"dependencies.yaml"
	script:
		"1_createMSDB.R"


rule scoringAndFiltering:
	input:
		sample_list = "INPUT/sample_list.csv",
		MSDB = "OUTPUT/{project_name}/tmp/MSDB.RData"
	output:
		allPSMs = "OUTPUT/{project_name}/tmp/allPSMs.RData",
		allPSMs_Delta = "OUTPUT/{project_name}/tmp/allPSMs_Delta.RData"
	params:
		project_name = config['project_name'],
		q_value = config['q_value'],
		ion_score = config['ion_score'],
		delta_score = config['delta_score']
	conda:
		"dependencies.yaml"
	script:
		"2_scoringAndFiltering.R"


rule removeSynErrors:
	input:
		allPSMs = "OUTPUT/{project_name}/tmp/allPSMs.RData"
	output:
		PSMs = "OUTPUT/{project_name}/tmp/extractedPSMs.RData"
	params:
		keep_synErrors = config["keep_synErrors"]
	conda:
		"dependencies.yaml"
	script:
		"3_removeSynErrors.R"


rule mapping:
	input:
		PSMs = "OUTPUT/{project_name}/tmp/extractedPSMs.RData",
		allPSMs_Delta = "OUTPUT/{project_name}/tmp/allPSMs_Delta.RData"
	output:
		ProteasomeDB = "OUTPUT/{project_name}/ProteasomeDB.csv",
		DeltaPeptides = "OUTPUT/{project_name}/tmp/DeltaPeptides.csv"
	conda:
		"dependencies.yaml"
	script:
		"4_mapping.R"


rule output_statistics:
	input:
		ProteasomeDB = "OUTPUT/{project_name}/ProteasomeDB.csv"
	output:
		DB_stats = "OUTPUT/{project_name}/DB_stats.pdf",
		number_of_products = "OUTPUT/{project_name}/number_of_products.pdf",
		length_distributions =  "OUTPUT/{project_name}/length_distributions.pdf"
	conda:
		"dependencies.yaml"
	script:
		"5_output_statistics.R"


rule plot_spectra:
	input:
		ProteasomeDB = "OUTPUT/{project_name}/ProteasomeDB.csv", 
		sample_list = "INPUT/sample_list.csv"
		
	output:
		spectra_plot="OUTPUT/{project_name}/ms2_spectra.pdf"
	params:
		mgf_folder = config['mgf_folder']
	conda:
		"dependencies.yaml"
	script:
		"plot_spectra.py"