
% run per sample
run_degradation_estimate('example_wt_tpm.txt','example_WT');
run_degradation_estimate('example_ko_tpm.txt','example_KO');

% compare WT and KO
load('example_data.mat','KO','KOtm','WT','WTtm','gids');
degradation_estimate_pairwise(gids,WTtm,WT,KOtm,KO,'example_pairwise');
