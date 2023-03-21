function RIDE_imagesc(results,comp_names)


[d1,d2] = size(results.erp);
n = length(comp_names);

for j = 1:n
    eval(['temp(:,:,j) = results.',comp_names{j},';']);
end
colorscal = [min(temp(:)),max(temp(:))];
    



for j = 1:n
    subplot(1,n,j); imagesc(1:d2,linspace(results.cfg.epoch_twd(1),results.cfg.epoch_twd(2),d1),temp(:,:,j));colormap(jet);
    caxis(colorscal);

    ylabel('Time after stimulus (ms)');
    xlabel('Channel');
    title(comp_names{j});
end

