% naive spike sorting

function [st1, st2] = spike_sorting_naive(x, rms_thres_std_factor, sampling_rate, hz_h)

% high pass filter
[z, p, g] = butter(10, hz_h / (sampling_rate/2), 'high');
sos = zp2sos(z, p, g);

Z = x;
for nn = 1:size(sos,1);
  %Z = filter(sos(nn,1:3),sos(nn,4:6), Z );
  Z = filtfilt(sos(nn,1:3), sos(nn,4:6), Z);
end
x1 = Z;

rms_thres = rms_thres_std_factor * std(x1(2000:end-2000), 1);

st1 = x1 > rms_thres;    % positive "spike"
st2 = x1 < -rms_thres;   % negative "spike"

if nargout <= 1
  st1 = st1 + st2;
end

end
