% read data file named type_str

function val = read_Intan_RHD2000_type(type_str, header, t_range)

is_full_range = ~exist('t_range', 'var');  % default to full range

[~, type_str, ext] = fileparts(type_str);  % allow adding suffix like amplifier.raw

switch type_str
  case 'time'
    dsize = [1, Inf];
    dtype = '*int32';
    dtypesize = 4;
    sr = header.frequency_parameters.amplifier_sample_rate;
    postproc = @(v) double(v) / sr;
  case 'amplifier'
    dsize = [length(header.amplifier_channels), Inf];
    dtype = '*int16';
    dtypesize = 2;
    sr = header.frequency_parameters.amplifier_sample_rate;
    postproc = @(v) 0.195 * double(v);
  %case 'auxiliary'
  %case 'supply'
  case 'analogin'
    dsize = [length(header.board_adc_channels), Inf];
    dtype = '*uint16';
    dtypesize = 2;
    sr = header.frequency_parameters.board_adc_sample_rate;
    postproc = @(v) 0.000050354 * double(v);  % convert to volts
    % For intan Recording Controller
    % postproc = @(v, h) double(v -32768) * 0.0003125; % convert to volts (Recording Controller)
  case 'digitalin'
    dsize = [1, Inf];
    dtype = '*uint16';
    dtypesize = 2;
    sr = header.frequency_parameters.board_dig_in_sample_rate;
    postproc = @(v) v;
    % usage
    % digital_input_ch = (bitand(digital_word, 2^ch) > 0); % ch has a value of 0-15 here
  otherwise
    error('invalid type');
end

% guess a path
fpath = [header.path, filesep, type_str, '.dat'];

if ~is_full_range
  % get byte range from time range
  fileinfo = dir(fpath);
  i_max = floor(fileinfo.bytes / dtypesize / dsize(1));
  i_range = convert_clamp_range(t_range, sr, i_max);
  byte_begin = (i_range(1)-1) * dsize(1) * dtypesize;
  dsize(2) = diff(i_range) + 1;
end

% read in binary data
fid = fopen(fpath, 'r');
if ~is_full_range
  fseek(fid, byte_begin, 'bof');
end
val = fread(fid, dsize, dtype);
fclose(fid);

if ~strcmp(ext, '.raw')
  val = postproc(val);
end

end

% get time range in unit of sample.
% the range is inclusive.
% idx=1 <=> t = 0
% idx=2 <=> t = 1 / sr
function i_rg = convert_clamp_range(t_range, sr, i_max)
  i_rg = floor(t_range(1:2) * sr) + 1;
  if i_rg(1) < 1
    i_rg(1) = 1;
  end
  if i_rg(2) > i_max
    i_rg(2) = i_max;
  end
end
