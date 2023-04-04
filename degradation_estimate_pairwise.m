function degradation_estimate_pairwise(gids,Xtm,X,Ytm,Y,result_dir)
% Calculate degradation rates in a transcriptional shut-off experiment
% Using a large set of stable genes for temporal normalization
%
% X = control experiment (with transcriptional shut-off)
% Y = treatment experiment (with transcriptional shut-off, after a
%     treatment or another perturbation)

linear_fit_type = 1; % generalized linear model
mkdir(result_dir);

% ----------------------------------------------------------------------
% minimal expression & log transform
% ----------------------------------------------------------------------
minEXP = 1e-2;

logX = log2(X);
logX(X < minEXP) = log2(minEXP);

logY = log2(Y);
logY(Y < minEXP) = log2(minEXP);

fprintf('input: %d genes\n', size(X,1));

% ----------------------------------------------------------------------
% select stable genes
% ----------------------------------------------------------------------
minR = 0.5;
maxD = -0.05;

[Dx,Xx,Rx] = expression_dg(Xtm,gids,logX,[],linear_fit_type,-1.5,1.5);
[Dy,Xy,Ry] = expression_dg(Ytm,gids,logY,[],linear_fit_type,-1.5,1.5);

fprintf('slopes (min, max):');
[{'Treatment' 'Control'}' num2cell([min([Dy Dx]);max([Dy Dx])])']
fprintf('log X(0) (min, max):');
[{'Treatment' 'Control'}' num2cell([min([Xy Xx]);max([Xy Xx])])']
fprintf('r squared (min, max):');
[{'Treatment' 'Control'}' num2cell([min([Ry Rx]);max([Ry Rx])])']

X0 = mean(logX(:,Xtm==0),2);
g = (Dx<=maxD).*(Rx>=minR).*(X0>1)==1;
fprintf('STABLE GENES: %d genes, Control experiment\n', sum(g));
if (sum(g) == 0)
    fprintf('NO STABLE GENES -- EXIT\n');
    return;
end

write_text_file([result_dir '/stable_genes.txt'],...
    [gids(g) num2cell([Dy(g,:) Dx(g,:)]) num2cell([logY(g,:) logX(g,:)])]);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
ylim = [-3 13];

subplot(2,2,1);
plot_gene_time(Xtm,logX(g,:),ylim,[0 0 0]);
title('stable genes: control');
subplot(2,2,2);
plot_gene_time(Ytm,logY(g,:),ylim,[0 0 0]);
title('stable genes: treatment');
subplot(2,2,3);
Ex = expression_norm(Xtm,logX,g);
plot_gene_time(Xtm,Ex(g,:),ylim,[0 0 0]);
title('stable genes: control (normalized)');
subplot(2,2,4);
Ey = expression_norm(Ytm,logY,g);
plot_gene_time(Ytm,Ey(g,:),ylim,[0 0 0]);
title('stable genes: treatment (normalized)');
saveas(h,[result_dir '/stable_genes.jpg'],'jpg');

close all;

% ----------------------------------------------------------------------
% normalize to stable genes, linear fit
% ----------------------------------------------------------------------
minD = 1e-2;
minR = 0.5;
FOLD = 1.5;

logXn = expression_norm(Xtm,logX,g);
[Dx,Xx,Rx] = expression_dg(Xtm,gids,logXn,[],linear_fit_type);
logDx = log2(Dx);
logDx(Dx<minD) = log2(minD);

logYn = expression_norm(Ytm,logY,g);
[Dy,Xy,Ry] = expression_dg(Ytm,gids,logYn,[],linear_fit_type);
logDy = log2(Dy);
logDy(Dy<minD) = log2(minD);

logD = [logDx logDy];
HL = log(2)./2.^logD;
Rsq = [Rx Ry];

T = {'id' 'half-life (C)' 'half-life (T)' 'r-sq (C)' 'r-sq (T)'};
write_text_file([result_dir '/degradation_rates.txt'],...
    [T;[gids num2cell([HL Rsq])]]);
o = sum(Rsq>=minR,2) == size(Rsq,2);
fprintf('Regression: %d genes with R_square > %.1f\n', sum(o), minR);
write_text_file([result_dir '/degradation_rates.rsq.txt'],...
    [T;[gids(o) num2cell([HL(o,:) Rsq(o,:)])]]);

q1 = (logDy-logDx >= log2(FOLD));
k = o.*(q1>0) == 1;
fprintf('Fold change: %d genes with degradation fold change > %.1f\n', sum(k), log2(FOLD));
write_text_file([result_dir '/degradation_rates.fasterdg.txt'],...
    [T;[gids(k) num2cell([HL(k,:) Rsq(k,:)])]]);
q2 = (logDy-logDx <= -1*log2(FOLD));
k = o.*(q2>0) == 1;
fprintf('Fold change: %d genes with degradation fold change > %.1f\n', sum(k), log2(FOLD));
write_text_file([result_dir '/degradation_rates.slowerdg.txt'],...
    [T;[gids(k) num2cell([HL(k,:) Rsq(k,:)])]]);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
dlim = [-6 -1];
xlim = [-2 12];

subplot(1,2,1);
x = dlim(1):0.2:dlim(2);
hold on;
y1 = hist(logDx,x);
plot(x,y1./sum(y1),'-','marker','.','markersize',20,'linewidth',1.5);
y2 = hist(logDy,x);
plot(x,y2./sum(y2),'-','marker','.','markersize',20,'linewidth',1.5);
hold off;
xlabel('degradation rate');
ylabel('frequency');
set(gca,'ylim',[0 0.25]);
legend({'Control' 'Treatment'},'box','off');
subplot(1,2,2);
x = xlim(1):0.5:xlim(2);
hold on;
y1 = hist(Xx,x);
plot(x,y1./sum(y1),'-','marker','.','markersize',20,'linewidth',1.5);
y2 = hist(Xy,x);
plot(x,y2./sum(y2),'-','marker','.','markersize',20,'linewidth',1.5);
hold off;
xlabel('X0');
ylabel('frequency');
set(gca,'ylim',[0 0.25]);
legend({'Control' 'Treatment'},'box','off');
saveas(h,[result_dir '/degradation_rates.jpg'],'jpg');

clf;
subplot(2,2,1);
plot_corr(dlim,logDx(o),logDy(o),'Control (rate)','Treatment (rate)');
subplot(2,2,2);
plot_corr(xlim,Xx(o),Xy(o),'Control (x0)','Treatment (x0)');
subplot(2,2,3);
plot_corr(dlim,logDx,logDy,'Control (rate)','Treatment (rate)');
subplot(2,2,4);
plot_corr(xlim,Xx,Xy,'Control (x0)','Treatment (x0)');
saveas(h,[result_dir '/degradation_rates.xy.jpg'],'jpg');

close all;


