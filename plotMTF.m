function fig2 = plotMTF(MTF, tempmod, specmod)

fig2 = figure('Color', 'w', 'Position', [30 50 1200 800]);
ax1 = axes('Position', [.15 .1 .75 .65]);
s2 = surface(ax1, tempmod,specmod,abs(MTF));
set(s2, 'FaceColor','interp','EdgeColor','interp');
xlim([-100 100])
ylim([0 10])
ylim([0 5])
cbar = colorbar;
colormap jet

MTF = abs(MTF);
[c_max, r_max] = max(MTF(:));
[col, row] = ind2sub([158 401], r_max);


ax2 = axes('Position', [.05 .1 .07 .65]);
plot(MTF(:,row), specmod)
title(ax2, 'sMTF')
cbar_lim = get(cbar, 'Limits');
cutoff = cbar_lim(2)*0.1;
set(gca, 'xdir', 'reverse', 'xlim', [cutoff cbar_lim(2)],...
    'ylim', [0 5], 'Visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')


ax3 = axes('Position', [.15 .78 .65 .07]);
plot(tempmod, MTF(col, :))
title(ax3, 'tMTF')
set(gca, 'xlim', [-100 100], 'ylim', [cutoff cbar_lim(2)], 'Visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
%sgtitle(['Site ',num2str(site), ' Chn ',num2str(channel),' MUA MTF'], 'Fontweight', 'bold')