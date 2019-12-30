% time series filter
% y = sosfiltb(x, 10, 200 / (sampling_rate/2), 'high');

function y = sosfiltx(x, od, freq, filtmode)

% high pass filter - 
[z, p, g] = butter(od, freq, filtmode);
sos = zp2sos(z, p, g);
%x1 = sosfilt(sos, x);

y = x;
for nn = 1:size(sos,1);
  %Z = filter(sos(nn,1:3),sos(nn,4:6), Z );
  y = filtfilt(sos(nn,1:3), sos(nn,4:6), y);
end

end
