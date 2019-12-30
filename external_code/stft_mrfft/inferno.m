function cm = inferno();
  pathdir0 = fileparts(mfilename('fullpath'));
  cm = csvread([pathdir0 filesep 'inferno.csv']);
end
