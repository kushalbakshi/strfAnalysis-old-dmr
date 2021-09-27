

MTF = abs(MTF(:));
[c_max, r_max] = max(MTF(:));
[col, row] = ind2sub([158 401], r_max);


ax2 = axes('Position', [.1 .1 .07 .65]);
plot(MTF(:,row), specmod)
set(gca, 'xdir', 'reverse', 'xlim', [0 6000], 'ylim', [0 5], 'Visible', 'off')
fig4 = figure;
plot(tempmod, MTF(col, :))

ax2 = findobj(fig2, 'type', 'axes');
ax3 = findobj(fig3, 'type', 'axes');
ax4 = findobj(fig4, 'type', 'axes');

ch3 = get(ax3, 'children');
ch4 = get(ax4, 'children');

copyobj(ch3, ax2);
copyobj(ch4, ax2);