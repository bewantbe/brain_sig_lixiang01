% statistics of spike-train frequency band
include_header

%  Delta波：0.5-4Hz;   0-4
%  Theta波：4-8Hz;     4-8
%  Alpha波：8-13Hz;    8-12
%  Beta波：13-32Hz;   12-40
%  Gamma波：>32Hz；   40-100
% https://cloud.tencent.com/developer/article/1532651
% https://www.jianshu.com/p/065b1ddeb8f3

hdpath = '/media/xyy/sam2t/home_extra/pmq/20191217_13_R_1841um_hip_shijue_191217_152712/info.rhd';

% load header
header = read_Intan_RHD2000_header(hdpath);

% load data
few1 = 1-eps;
t_range = 10+[0, 500*few1];
%amp_s = read_Intan_RHD2000_type('amplifier', header, t_range);

sampling_rate = header.frequency_parameters.amplifier_sample_rate;

% extract partial data
id_ch = 19;
rg_stat = sampling_rate*100 + (1:sampling_rate*60);
x = amp_s(id_ch, rg_stat);

hz_h = 200;
rms_thres_std_factor = 3.0;
st = spike_sorting_naive(x, rms_thres_std_factor, sampling_rate, hz_h);

bw_method = 2;
switch bw_method
  case 1
    st_alpha = sosfiltx(st, 20, [8 13] / (sampling_rate/2), 'bandpass');
    st_beta  = sosfiltx(st, 20, [13 32] / (sampling_rate/2), 'bandpass');
    st_gamma = sosfiltx(st, 20, [32 100] / (sampling_rate/2), 'bandpass');

    % STFT parameters
    sz_wnd = 16382;
    sz_hop = sz_wnd/2;

    figure(200);
    plot_stft(st_alpha, sz_wnd, sz_hop, sz_wnd, 'hanning', sampling_rate);
    caxis([-70 30]);
  case 2
    sr_local = 100;
    st100 = resample(st, sr_local, sampling_rate);  % down-sampling to 100 Hz
    st_alpha = sosfiltx(st100, 10, [8 13] / (sr_local/2), 'bandpass');
    st_beta  = sosfiltx(st100, 10, [13 32] / (sr_local/2), 'bandpass');
    st_gamma = sosfiltx(st100, 10, 32 / (sr_local/2), 'high');

    % STFT parameters
    sz_wnd = 128;
    sz_hop = sz_wnd/2;

    figure(201);
    plot_stft(st_alpha, sz_wnd, sz_hop, sz_wnd, 'hanning', sr_local);
    caxis([-70 30]);

    figure(202);
    plot_stft(st_beta , sz_wnd, sz_hop, sz_wnd, 'hanning', sr_local);
    caxis([-70 30]);

    figure(203);
    plot_stft(st_gamma, sz_wnd, sz_hop, sz_wnd, 'hanning', sr_local);
    caxis([-70 30]);
    
    figure(301);
    plot((1:length(st100))/sr_local, st_alpha);
    xlabel('time/s');
    ylabel('alpha band');
    
    figure(302);
    plot((1:length(st100))/sr_local, st_beta);
    xlabel('time/s');
    ylabel('alpha band');
    
    figure(303);
    plot((1:length(st100))/sr_local, st_gamma);
    xlabel('time/s');
    ylabel('alpha band');
    
end

