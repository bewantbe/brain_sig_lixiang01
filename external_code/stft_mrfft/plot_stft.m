% plot STFT

function [xtf, tc, fc] = plot_stft(x, sz_wnd, sz_hop, sz_fft, wnd_name, sampling_rate)

wnd = select_window(wnd_name, sz_wnd).';
wnd = wnd / (wnd' * wnd);
wnd = 2 * wnd / sum(wnd);  % normalize to 0dB = 1 V-peak amplitude (sine wave)

x = x - mean(x);

[xtf, tc, fc] = stft_wnd(x, wnd, sz_hop, sz_fft);

amp2db = @(x) 20*log10(abs(x));

% Note that in even fft length, the highest frequency component should cut by half.
db = amp2db(xtf(1:ceil((end+1)/2), :));
if mod(sz_fft, 2) == 0
  db(end, :) = db(end, :) - 20*log10(2);
end
fc = fc * sampling_rate;
tc = tc / sampling_rate;
fc_show = fc(1:ceil((end+1)/2));

imagesc(tc, fc_show, db);
set(gca,'YDir','normal');

colormap(inferno());
hcb = colorbar();
set(get(hcb, 'title'), 'String','dB')

xlabel('T');
ylabel('Freq');

% vim: set expandtab shiftwidth=2 softtabstop=2:
