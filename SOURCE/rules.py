
rule createMSDB:
	input:
		sample_list = "INPUT/sample_list.csv"
	output:
		MSDB = "OUTPUT/tmp/MSDB.RData"
	conda:
		"dependencies.yaml"
	script:
		"1_createMSDB.R"


rule scoringAndFiltering:
	input:
		sample_list = "INPUT/sample_list.csv",
		MSDB = "OUTPUT/tmp/MSDB.RData"
	output:
		allPSMs = "OUTPUT/tmp/allPSMs.RData"
	params:
		q_value = config['q_value'],
		ion_score = config['ion_score'],
		delta_score = config['delta_score']
	conda:
		"dependencies.yaml"
	script:
		"2_scoringAndFiltering.R"


rule removeSynErrors:
	input:
		allPSMs = "OUTPUT/tmp/allPSMs.RData"
	output:
		PSMs = "OUTPUT/tmp/extractedPSMs.RData"
	params:
		keep_synErrors = config["keep_synErrors"]
	conda:
		"dependencies.yaml"
	script:
		"3_removeSynErrors.R"


rule mapping:
	input:
		PSMs = "OUTPUT/tmp/extractedPSMs.RData"
	output:
		ProteasomeDB = "OUTPUT/ProteasomeDB.csv"
	conda:
		"dependencies.yaml"
	script:
		"4_mapping.R"

