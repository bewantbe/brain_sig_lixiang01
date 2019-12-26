%

amp_cov = amp_s * amp_s' / size(amp_s,2);
figure(10);
imagesc(amp_cov);
colorbar
title('cov matrix');

amp_mean = mean(amp_s');
amp_stdvar = std(amp_s');

figure(14);
errorbar(amp_mean, amp_stdvar);
title('channel mean and stdvar');

amp_s_zscore = zscore(amp_s, 0, 2);
amp_cor = amp_s_zscore * amp_s_zscore' / size(amp_s_zscore,2);
figure(15);
imagesc(amp_cor);
colorbar
title('cor matrix');

