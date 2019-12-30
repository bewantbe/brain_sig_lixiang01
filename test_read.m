%
include_header

hdpath = '/media/xyy/sam2t/home_extra/pmq/20191217_13_R_1991um_hip_191217_164233/info.rhd';

few1 = 1-eps;
t_range = 10+[0, 500*few1];

header = read_Intan_RHD2000_loader(hdpath);
tic
t_s = read_Intan_RHD2000_type('time', header, t_range);
amp_s = read_Intan_RHD2000_type('amplifier', header, t_range);
ana_s = read_Intan_RHD2000_type('analogin',  header, t_range);
dig_s = read_Intan_RHD2000_type('digitalin', header, t_range);
tocs('(load data)');

%amp_cov = amp_s * amp_s' / size(amp_s,2);

%figure(10);
%imagesc(amp_cov);
%colorbar
%title('cov matrix');

%figure(9);
%plot(t_s, ana_s);
%title('analogin');
