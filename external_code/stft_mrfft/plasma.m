function cm = plasma();
  pathdir0 = fileparts(mfilename('fullpath'));
  cm = csvread([pathdir0 filesep 'plasma.csv']);
end
