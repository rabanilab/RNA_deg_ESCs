function plot_gene_time(tm,Z,ylim,nclr)

hold on;
x = [];
y = [];
for i = 1:size(Z,1)
    x = [x tm];
    y = [y Z(i,:)];
end
scatter(x,y,'filled','MarkerFaceColor',nclr);
alpha(0.2);
for i = 1:size(Z,1)
    [u,~,s] = unique(tm);
    m = accumarray(s,Z(i,:));
    n = accumarray(s,1);
    plot(u,m./n,'-','linewidth',1.2,'Color',[nclr 0.2]);
end
ut = unique(tm);
n = max(size(ut));
M = zeros(1,n);
for i = 1:n
    x = Z(:,tm == ut(i));
    M(i) = mean(x(:));
end
plot(ut,M,'-r','linewidth',2);
hold off;
set(gca,'ylim',ylim);
xlabel('time');
ylabel('log EXP');
