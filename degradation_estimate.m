function degradation_estimate(gids,Xtm,X,result_dir)
% Calculate degradation rates in a transcriptional shut-off experiment
% Using a large set of stable genes for temporal normalization
%
% X = experiment (with transcriptional shut-off)

linear_fit_type = 1; % generalized linear model
mkdir(result_dir);

% ----------------------------------------------------------------------
% minimal expression & log transform
% ----------------------------------------------------------------------
minEXP = 1e-2;

logX = log2(X);
logX(X < minEXP) = log2(minEXP);

fprintf('input: %d genes\n', size(X,1));

% ----------------------------------------------------------------------
% select stable genes
% ----------------------------------------------------------------------
minR = 0.5;
maxD = -0.05;

[Dx,Xx,Rx] = expression_dg(Xtm,gids,logX,[],linear_fit_type,-1.5,1.5);

[[{'slopes (min, max):'} num2cell([min(Dx);max(Dx)])'];
 [{'log X(0) (min, max):'} num2cell([min(Xx);max(Xx)])'];
 [{'r squared (min, max):'} num2cell([min(Rx);max(Rx)])']]

X0 = mean(logX(:,Xtm==0),2);
g = (Dx<=maxD).*(Rx>=minR).*(X0>1)==1;
fprintf('STABLE GENES: %d genes, Control experiment\n', sum(g));
if (sum(g) == 0)
    fprintf('NO STABLE GENES -- EXIT\n');
    return;
end

write_text_file([result_dir '/stable_genes.txt'],...
    [gids(g) num2cell(Dx(g,:)) num2cell(logX(g,:))]);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
ylim = [-3 13];

subplot(1,2,1);
plot_gene_time(Xtm,logX(g,:),ylim,[0 0 0]);
title('stable genes');
subplot(1,2,2);
Ex = expression_norm(Xtm,logX,g);
plot_gene_time(Xtm,Ex(g,:),ylim,[0 0 0]);
title('stable genes (normalized)');
saveas(h,[result_dir '/stable_genes.jpg'],'jpg');

close all;

% ----------------------------------------------------------------------
% normalize to stable genes, linear fit
% ----------------------------------------------------------------------
minD = 1e-2;
minR = 0.5;

logXn = expression_norm(Xtm,logX,g);
[Dx,Xx,Rx] = expression_dg(Xtm,gids,logXn,[],linear_fit_type);
logDx = log2(Dx);
logDx(Dx<minD) = log2(minD);

logD = logDx;
HL = log(2)./2.^logD;
Rsq = Rx;

T = {'id' 'half-life' 'r-sq'};
write_text_file([result_dir '/degradation_rates.txt'],...
    [T;[gids num2cell([HL Rsq])]]);
o = sum(Rsq>=minR,2) == size(Rsq,2);
fprintf('Regression: %d genes with R_square > %.1f\n', sum(o), minR);

write_text_file([result_dir '/degradation_rates.rsq.txt'],...
    [T;[gids(o) num2cell([HL(o,:) Rsq(o,:)])]]);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
dlim = [-6 -1];
xlim = [-2 12];

subplot(1,2,1);
x = dlim(1):0.2:dlim(2);
y = hist(logDx,x);
plot(x,y./sum(y),'-','marker','.','markersize',20,'linewidth',1.5);
hold on;
y = hist(logDx(o),x);
plot(x,y./sum(y),'-','marker','.','markersize',20,'linewidth',1.5);
hold off;
xlabel('degradation rate');
ylabel('frequency');
set(gca,'ylim',[0 0.25]);
legend({'All' sprintf('R-sq > %.1f',minR)},'box','off');
subplot(1,2,2);
x = xlim(1):0.5:xlim(2);
y = hist(Xx,x);
plot(x,y./sum(y),'-','marker','.','markersize',20,'linewidth',1.5);
hold on;
y = hist(Xx(o),x);
plot(x,y./sum(y),'-','marker','.','markersize',20,'linewidth',1.5);
hold off;
xlabel('X0');
ylabel('frequency');
set(gca,'ylim',[0 0.25]);
legend({'All' sprintf('R-sq > %.1f',minR)},'box','off');
saveas(h,[result_dir '/degradation_rates.jpg'],'jpg');
close all;

