function visResiduals_woshow(We, outlWe, t, Residuals, ktraj, folderName, fileName)

co = [0 0.4470 0.7410;
       0.8500 0.3250 0.0980];
FontSizeA = 18;
FontSizeB = 16;
LineWidth = 2;
MarkerSize = 12;

if ~iscell(folderName)
    folderName = {folderName, folderName};
end

if ~isempty(We)
    h1 = figure('Visible', 'off');
    if ~isempty(t); xlab = '$t$ (s)'; else; xlab = '$t$ (au)'; end
    if isempty(t); t = 1:numel(We); end
    plot(t, We(:), 'Color', co(1,:), 'LineWidth', LineWidth, 'LineStyle', '-')
    hold on
    if ~isempty(outlWe)
        plot(t(outlWe), We(outlWe), '*', 'Color', co(2,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
    end
    if t(end) > t(1); xlim([t(1), t(end)]); end
    xlabel(xlab, 'Interpreter', 'latex', 'FontSize', FontSizeA)
    ylabel('Weight', 'Interpreter', 'latex', 'FontSize', FontSizeA)
    grid on
    set(gca, 'FontSize', FontSizeA)
    legendEntries = {'Weights'};
    if any(outlWe(:) ~= 0); legendEntries{2} = 'Outliers'; end
    legend(legendEntries, 'Location', 'SouthEast', 'FontSize', FontSizeA - 4)
    set(h1, 'Color', 'w', 'Position', [100, 100, 1200, 800])

    if ~isempty(folderName{1}) && ~isempty(fileName)
        if ~exist(folderName{1}, 'dir'); mkdir(folderName{1}); end
        savePath1 = fullfile(folderName{1}, [fileName, '.png']);
        print(h1, savePath1, '-dpng', '-r300');
    end
    close(h1);
end

if ~isempty(Residuals)
    h2 = figure('Visible', 'off');
    if ~isempty(ktraj)
        kMin = min(ktraj, [], 1);
        kMax = max(ktraj, [], 1);
        for n = 1:2; kGrid{n} = kMin(n):kMax(n); end
        kGrid{1} = kGrid{1}';
        kGrid{2} = repmat(kGrid{2}, [1 size(Residuals, 3)]);
        for m = 1:2; Residuals = fftshift(Residuals, m); end
        Residuals = log(Residuals(:,:));
        imagesc('XData', kGrid{2}, 'YData', kGrid{1}, 'CData', Residuals);
    else
        imagesc('CData', Residuals);
    end
    colormap(jet)
    axis image
    title('log($r$)', 'Interpreter', 'latex', 'FontSize', FontSizeB)
    xlabel('$k_3$', 'Interpreter', 'latex', 'FontSize', FontSizeB)
    ylabel('$k_2$', 'Interpreter', 'latex', 'FontSize', FontSizeB, 'Rotation', 0)
    set(h2, 'Color', 'w', 'Position', [100, 100, 1200, 800])

    if ~isempty(folderName{2}) && ~isempty(fileName)
        if ~exist(folderName{2}, 'dir'); mkdir(folderName{2}); end
        savePath2 = fullfile(folderName{2}, [fileName, '.png']);
        print(h2, savePath2, '-dpng', '-r300');
    end
    close(h2);
end
end
