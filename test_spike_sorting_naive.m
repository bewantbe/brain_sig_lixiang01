% test spike_sorting_naive
% run test_read first

sampling_rate = header.frequency_parameters.amplifier_sample_rate;
hz_h = 200;
rms_thres_std_factor = 3.0;

id_ch = 19;
rg_stat = sampling_rate*100 + (1:sampling_rate*60);
x = amp_s(id_ch, rg_stat);
rms_x = std(x, 1);

[st1, st2] = spike_sorting_naive(x, rms_thres_std_factor, sampling_rate, hz_h);

rg_show = 1*sampling_rate + (1:floor(0.5*sampling_rate));

% show curve - spike train comparison
figure(100);
t_demo = (rg_show - rg_show(1))/sampling_rate;
plot(t_demo, 2*rms_x * st1(rg_show), ...\
     t_demo,-2*rms_x * st2(rg_show), ...
     t_demo, x(rg_show), 'g');

