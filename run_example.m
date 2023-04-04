
load('example_data.mat','KO','KOtm','WT','WTtm','gids');

degradation_estimate(gids,WTtm,WT,'example_WT');
degradation_estimate(gids,KOtm,KO,'example_KO');
degradation_estimate_pairwise(gids,WTtm,WT,KOtm,KO,'example_pairwise');
