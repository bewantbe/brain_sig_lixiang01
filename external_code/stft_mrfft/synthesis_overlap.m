% Inverse stft, by 
% use overlap and add method
% 1-dim is frequency, 2-dim is time;

function x = synthesis_overlap(xtf, wnd, sz_hop)

sz_dat = length(wnd);     % data(time) window size
sz_fft = sz_dat;          % FFT size
if ~exist('sz_hop', 'var')
  sz_hop = floor(sz_dat/2);  % hop size
end

xw = real(ifft(xtf));  % assume real signal

ids = matrep((1:sz_hop).', 1);

id0 = 1:sz_hop;
for k = 0:n_hop-1
  id = id0;
  x_slice = zeros(1, sz_hop);
  j = 0;
  while id(end) < 1
    x_slice = x_slice + xw(id(id>=1 & id < sz_dat), k+j+1);
    id = id - sz_hop;
    j = j + 1;
  end
  id0 = id0 + sz_hop;
  x(sz_hop*k+1 : sz_hop*k+sz_hop) = x_slice;
end


end

% vim: set expandtab shiftwidth=2 softtabstop=2:
