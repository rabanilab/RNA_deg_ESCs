update_permissions:
	chmod +x *.sh; \
	chmod +x run_degradation_estimate;

run_dg_estimate:
	./matlab_run.sh $(TPM_FILE) $(OUTPUT_DIR);

run_WT_example:
	make run_dg_estimate TPM_FILE=example_wt_tpm.txt OUTPUT_DIR=example_wt;

run_KO_example:
	make run_dg_estimate TPM_FILE=example_ko_tpm.txt OUTPUT_DIR=example_ko;

###########################################################################
#                              kmer analysis                              #
###########################################################################

-include cont_kmer_analysis/makefile

download_kmer_analysis_package:
	git clone https://github.com/rabanilab/cont_kmer_analysis;

unzip_example:
	cd cont_kmer_analysis; \
	make unzip_kmer_files ORGANISM=mmusculus_3UTR_ensembl99; \
	cd ..;

convert_file_format:
	Rscript deg_out_to_kmer.R $(path);

convert_WT_file_example:
	make convert_file_format path=example_wt/degradation_rates.rsq.txt;

convert_KO_file_example:
	make convert_file_format path=example_ko/degradation_rates.rsq.txt;

kmer_WT_example:
	cd cont_kmer_analysis; \
	make run_ks_test path_to_kmer_matrices=mmusculus_3UTR_ensembl99/3utr/kmer_matrices_comp PARAMETERS_FILE=../example_wt/degradation_rates.rsq.kmer.tsv output_path=../example_wt term=all; \
	make ks_tests_plots ks_test_output_folder=../example_wt term=all output_path=../example_wt; \
	cd ..;

kmer_KO_example:
	cd cont_kmer_analysis; \
	make run_ks_test path_to_kmer_matrices=mmusculus_3UTR_ensembl99/3utr/kmer_matrices_comp PARAMETERS_FILE=../example_ko/degradation_rates.rsq.kmer.tsv output_path=../example_ko term=all; \
	make ks_tests_plots ks_test_output_folder=../example_ko term=all output_path=../example_ko; \
	cd ..;

run_complete_analysis:
	echo "Running degradation estimate..."; \
	make run_dg_estimate TPM_FILE=$(input_tpm_file).txt OUTPUT_DIR=$(out_dir); \
	echo "downloading kmer package..."; \
	make download_kmer_analysis_package; \
	cd cont_kmer_analysis; \
	make unzip_kmer_files ORGANISM=$(ensembl_organism); \
	cd ..; \
	echo "Running kmer analysis..."; \
	make convert_file_format path=$(out_dir)/degradation_rates.rsq.txt; \
	cd cont_kmer_analysis; \
	make run_ks_test path_to_kmer_matrices=$(ensembl_organism)/3utr/kmer_matrices_comp PARAMETERS_FILE=../$(out_dir)/degradation_rates.rsq.kmer.tsv output_path=../$(out_dir) term=all; \
	make ks_tests_plots ks_test_output_folder=../$(out_dir) term=all output_path=../$(out_dir); \
	cd ..; \
	echo "Done!";
