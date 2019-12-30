% common header
if exist('included_header', 'var') && included_header
  return
end
included_header = true;

pkg load signal

% for plot_stft
addpath('./external_code/stft_mrfft/');

global pic_output
global pic_output_eps

exp_name = 'M_AD';

pic_output     = @(st) print('-dpng',   ['pic_tmp/' exp_name '_' st '.png']);
pic_output_eps = @(st) print('-depsc2', ['pic_tmp/' exp_name '_' st '.eps']);

tocs = @(st) fprintf('Time: %7.3f sec: %s\n', toc(), st);

tocids = @(st, tid) fprintf('Time: %7.3f sec: %s\n', toc(tid), st);

