function plot_corr(dlim,D1,D2,n1,n2)

j = (~isnan(D1)).*(~isnan(D2)) == 1;
if (sum(j)>=5)
    hold on;
    x1 = D1(j);
    x1(x1<dlim(1)) = dlim(1);
    x1(x1>dlim(2)) = dlim(2);
    x2 = D2(j);
    x2(x2<dlim(1)) = dlim(1);
    x2(x2>dlim(2)) = dlim(2);
    dscatter(x1,x2,'MARKER','o','MSIZE',35);
    line(dlim,dlim,'LineStyle','-','color','k','linewidth',1);
    hold off;
    axis square;
    set(gca,'ylim',dlim,'xlim',dlim);
    xlabel(n1);
    ylabel(n2);
end
r = corr(x1,x2);
title(sprintf('r = %.2f (n = %d)',r,sum(j)));
