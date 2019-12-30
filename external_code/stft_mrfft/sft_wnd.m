%

function [aveS, fqs] = sft_wnd(X, wnd, sz_hop, sz_fft, spectrum_normalize)
[~, n_var, n_trials] = size(X);

sz_dat = length(wnd);
if ~exist('sz_fft', 'var') || isempty(sz_fft)
  sz_fft = sz_dat;           % FFT size
end
if ~exist('sz_hop', 'var') || isempty(sz_hop)
  sz_hop = floor(sz_dat/2);  % hop size
end
if ~exist('spectrum_normalize', 'var')
  spectrum_normalize = 'statistics';
end

n_hop = floor((size(X,1)-sz_dat)/sz_hop) + 1;

aveS = zeros(sz_fft,n_var,n_var);
S    = zeros(sz_fft,n_var,n_var);
% average over trials
for i_trial=1:n_trials
  for l = 0:n_hop-1
    bg = round(l*sz_hop);
    % windowed Fourier transform
    Jk = fft(bsxfun(@times, wnd, X(bg+1:bg+sz_fft,:,i_trial)));
    % get cross spectrum of one slice
    for chan1=1:n_var
      S(:, chan1, chan1) = Jk(:,chan1).*conj(Jk(:,chan1));
      for chan2=chan1+1:n_var
        S(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
        S(:, chan2, chan1) = conj(S(:, chan1, chan2));
      end
    end
    aveS = aveS + S;
  end
end

switch spectrum_normalize
  case {'', 'statistics'}
    fact = wnd'*wnd;
  case {'audio'}
    fact = sum(wnd)^2 / 4;
  otherwise
    error('invalid value');
end

aveS = aveS / (n_trials * n_hop * fact);

fqs = [0:floor(sz_fft/2), -floor(sz_fft/2)+mod(sz_fft+1,2):-1]'/sz_fft;

