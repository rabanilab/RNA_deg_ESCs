function [DG,X0,RS] = expression_dg(T,G,logEXP,cnt_i,model,LB,UB,output_pref,nplot)

if (nargin < 4)
    cnt_i = [];
end
if (nargin < 5)
    model = 0; % default, polynomial fit
end
if (nargin < 6)
    LB = 1e-2;
end
if (nargin < 7)
    UB = 4;
end
if (nargin < 8)
    output_pref = [];
end
if (nargin < 9)
    nplot = 10;
end

% normalize
n = size(logEXP,1);
if (isempty(cnt_i))
    nEXP = logEXP;
else
    nEXP = expression_norm(T,logEXP,cnt_i);
end

% calculate degradation
DG = zeros(n,1);
X0 = zeros(n,1);
RS = zeros(n,1);
for i = 1:n
    if (model == 1)
        [DG(i),X0(i),RS(i)] = fit_glm(T,nEXP(i,:),LB,UB);
    else
        [DG(i),X0(i),RS(i)] = fit_poly(T,nEXP(i,:),LB,UB);
    end
end

% plot examples
if (~isempty(output_pref))
    h = figure;
    fk = sortrows([(1:size(RS,1))' RS],-2);
    for j = fk(1:nplot,1)'
        clf;
        hold on;
        plot(T,nEXP(j,:),'or','markersize',10);
        z = polyval([-1*DG(j) X0(j)],T);
        plot(T,z,'-k','linewidth',1.5);
        hold off;
        set(gca,'ylim',[0 15]);
        xlabel('time');
        ylabel('log TPM');
        hl = log(2)./DG(j);
        title(sprintf('%s (half-life = %.1f hr, rsq = %.2f)',G{j},hl,RS(j)));
        saveas(h,[output_pref '.' G{j} '.jpg'],'jpg');
    end
    close all;
end

function [dg,x0,rs] = fit_poly(x,y,LB,UB)

P = polyfit(x,y,1);
if (-1*P(1)<LB)
    dg = LB;
    x0 = nanmean(y,2);
elseif (-1*P(1)>UB)
    dg = UB;
    x0 = P(2);
else
    dg = -1*P(1);
    x0 = P(2);
end
ym = polyval([-1*dg x0],x);
rs = r_squared(y,ym);

function [dg,x0,rs] = fit_glm(x,y,LB,UB)

[P,s] = glmfit(x,exp(y),'poisson','link','log');
if (-1*P(2)<LB)
    dg = LB;
    x0 = nanmean(y,2);
elseif (-1*P(2)>UB)
    dg = UB;
    x0 = P(1);
else
    dg = -1*P(2);
    x0 = P(1);
end
ym = glmval(P,x,'log');
rs = r_squared(y,log(ym'));

function r2 = r_squared(y, ym)
% y = observation
% ym = model predictions

[m,n] = size(y);
mu = repmat(mean(y,2),1,n);

SSerr = sum((y - ym).^2, 2);
SStot = sum((y - mu).^2, 2);
r2 = ones(m,1) - SSerr./SStot;

r2((SStot==0)+(SSerr==0)==2) = 1;
r2((SStot==0)+(SSerr~=0)==2) = 0;
r2(abs(r2)>1) = r2(abs(r2)>1)./abs(r2(abs(r2)>1));
