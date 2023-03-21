function f = align_c(c)

latency = zeros(size(c,3),1);
temp5 = [];temp11 = [];
for iter = 1:5
    grand_c = mean(move3(c,-latency),3);
    for j = 1:size(c,3)
        for ch = 1:size(c,2)
            temp11(:,ch) = xcorr(c(:,ch,j),grand_c(:,ch),'coeff');
        end
        temp5(:,j) = mean(temp11,2);
    end
    figure;subplot(1,2,1);imagesc(temp5);subplot(1,2,2);imagesc(grand_c);
    latency = find_peak2(temp5);
    latency = round(latency-median(latency));
end

f = move3(c,-latency);