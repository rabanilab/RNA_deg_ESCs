
run_dg_estimate:
	./matlab_run.sh $(TPM_FILE) $(OUTPUT_DIR);

run_WT_example:
	make run_dg_estimate TPM_FILE=example_wt_tpm.txt OUTPUT_DIR=example_wt;

run_KO_example:
	make run_dg_estimate TPM_FILE=example_ko_tpm.txt OUTPUT_DIR=example_ko;
