figure; hold on;
AA('0.0_com', '^g');
AA('0.5_com', 'ob');
AA('1.0_com', '<m');
AA('1.5_com', 'sr');
AA('2.0_com', '>c');
AA('2.5_com', 'vk');
AA('3.0_com', 'xg');
xlim([0, 3.5]);
ylim([0.15, 0.55]);
xlabel('Initial KE (J)')
ylabel('Final Depth (m)')
legend(
    's = 0.0',
    's = 0.5',
    's = 1.0',
    's = 1.5',
    's = 2.0',
    's = 2.5',
    's = 3.0', 'Location', 'southeast'
)
saveas(gcf, 'new_coarse_depths.png')
hold off;
