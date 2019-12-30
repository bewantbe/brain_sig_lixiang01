%

sampling_rate = header.frequency_parameters.amplifier_sample_rate;

id_ch = 19;
x = amp_s(id_ch, (1:sampling_rate*60)+sampling_rate*100);
t_l = t_s((1:sampling_rate*60)+sampling_rate*100);

rg_show = 37*sampling_rate + (1:floor(0.5*sampling_rate));
rg_stat = 4000:length(x)-4000;  % remove head and tail for statistics

% STFT parameters
sz_wnd = 2048;
sz_hop = sz_wnd/2;

if 1
  figure(201);
  plot_stft(x, sz_wnd, sz_hop, sz_wnd, 'hanning', sampling_rate);
  caxis([-25 30]);
  %pic_output(sprintf('volt_intv=%d_ch=%d_stft', id_intv, id_ch));
end


figure(10010);
plot(t_l(rg_show) - t_l(1), x(rg_show));

% high pass filter - 
[z, p, g] = butter(10, 200 / (sampling_rate/2), 'high');
sos = zp2sos(z, p, g);
%x1 = sosfilt(sos, x);

Z = x;
for nn = 1:size(sos,1);
  %Z = filter(sos(nn,1:3),sos(nn,4:6), Z );
  Z = filtfilt(sos(nn,1:3), sos(nn,4:6), Z);
end
x1 = Z;


% for ARregression
%addpath('./external_code/GC_clean/GCcal');

%ar_od = 6;
%b = ARregressionpd(getcovzpd(x, ar_od), size(x,1));
%x2 = filter([1 b], 1, x);
%x1 = x2;

figure(211);
rms_thres = 4*std(x1, 1);
rms_thres = 50;
plot(t_l - t_l(1), abs(x1) > rms_thres);

figure(10011);
plot(t_l(rg_show) - t_l(1), x1(rg_show), ...
     t_l(rg_show([1 end])) - t_l(1), rms_thres * [1 1], ...
     t_l(rg_show([1 end])) - t_l(1), -rms_thres * [1 1]);
ylim(150*[-1 1])

if 1
  figure(202);
  plot_stft(x1, sz_wnd, sz_hop, sz_wnd, 'hanning', sampling_rate);
  caxis([-25 30]);
  %pic_output(sprintf('volt_intv=%d_ch=%d_stft', id_intv, id_ch));
end

%% get var from mean and var
%m1 = mean(x1(rg_stat));
%v1 = var(x1(rg_stat));

% get var from quantile (looks better)
q = quantile(x1(rg_stat), [0.2 0.8]);
q01 = norminv([0.2 0.8]);
m1 = mean(q);
v1 = (diff(q) / diff(q01)) .^ 2;

figure(210)
x_hist = linspace(-150, 150, 1000);
hist(x1(rg_stat), x_hist);
hold on
dx = x_hist(2) - x_hist(1);
plot(x_hist, length(rg_stat)*dx * 1/sqrt(2*pi*v1) * exp(-(x_hist - m1).^2 / (2*v1)), '-r')
hold off


