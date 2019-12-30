function cm = magma();
  pathdir0 = fileparts(mfilename('fullpath'));
  cm = csvread([pathdir0 filesep 'magma.csv']);
end
