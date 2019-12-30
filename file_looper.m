% file loop checker

dir_base = '/media/xyy/sam2t/home_extra/pmq/';  % the tailing '/' is important

file_list = dir(dir_base);

for k = 1 : length(file_list)
  if ~(filelist(k).isdir && filelist(k).name(1) ~= '.')
    continue  % skip non-data directory
  end
  fprintf('checking: %s\n', filelist(k).name);
  hdpath = [dir_base filelist(k).name '/info.rhd'];
  header = read_Intan_RHD2000_header(hdpath, false);
  dig_s = read_Intan_RHD2000_type('digitalin', header);
  % check digitalin channel
  fprintf('dig_s(1)            : %f\n', dig_s(1));
  % zero means data do not change
  fprintf('sum dig_s - dig_s(1): %f\n', sum(dig_s - dig_s(1)));
  fprintf('\n');
end

