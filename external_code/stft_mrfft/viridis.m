function cm = viridis();
  pathdir0 = fileparts(mfilename('fullpath'));
  cm = csvread([pathdir0 filesep 'viridis.csv']);
end
