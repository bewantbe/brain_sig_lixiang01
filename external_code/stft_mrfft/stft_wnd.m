% stft (with phase)
% stft of 1-dim of x
% xtf: 1-dim is frequency, 2-dim is time, 3-dim is trial.

function [xtf, tc, fc] = stft_wnd(x, wnd, sz_hop, sz_fft)

% Safe guard about dimension of x
if size(x, 2) > size(x, 1)
  if size(x, 1) == 1
    x = x.';
  else
    warning('1-dim of x is time!');
  end
end

% Safe guard about dimension of wnd
if size(wnd, 1) == 1
  wnd = wnd.';
else
  if size(x,2) ~= size(wnd,2)
    if size(x,2) == size(wnd,1)
      wnd = wnd.';
    else
      error('Dimension of data taper (window function) not recognize.');
    end
  end
end

sz_dat = length(wnd);        % data(time) window size
if ~exist('sz_fft', 'var')
  sz_fft = sz_dat;           % FFT size
end
if ~exist('sz_hop', 'var')
  sz_hop = floor(sz_dat/2);  % hop size
end

% zero padding to make (length(x) - sz_dat) the multiple of sz_hop
x = [x; zeros(sz_hop - 1 - mod(size(x, 1) - sz_dat - 1, sz_hop), size(x,2))];

n_hop = (size(x,1)-sz_dat)/sz_hop + 1;
xtf = zeros(sz_fft, n_hop, size(x, 2));

% STFT (naive implementation)
for l = 0 : n_hop - 1
  u = x(1 + (0:sz_dat-1) + l * sz_hop, :);
  xtf(:, l+1, :) = fft(u .* wnd, sz_fft);
end

% time point of each frame (center)
tc = (0 : n_hop - 1) * sz_hop + sz_dat/2;

% frequency in each frame
fc = (1:sz_fft) - ceil(sz_fft/2);
mid = floor((sz_fft-1)/2);
fc = [fc(mid+1:end), fc(1:mid)];
fc = fc / sz_fft;

% vim: set expandtab shiftwidth=2 softtabstop=2:
