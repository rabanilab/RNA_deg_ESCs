function run_degradation_estimate(tpmfile,output)

X = importdata(tpmfile);
gids = X.textdata(2:end,1);
tm = X.data(1,2:end);
TPM = X.data(2:end,2:end);

degradation_estimate(gids,tm,TPM,output);
