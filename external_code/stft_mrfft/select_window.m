% Return a window function
%
% wnd_name : Name of window funciton. See the code for details.
%            If there is parameter, pass a cell to it. See the example.
% len      : Length of window function
%            When omited, function handler is returned.
%            Omission is possible if there is explicit function form.
%
% Usage Example:
%
%  plot( [0.193*select_window({'parzen'}, 100); select_window({'dpss'; 3.6}, 100)']' )

function f_wnd = select_window(wnd_name, len, b_ret_handler_if_possible)
  if iscell(wnd_name)
    wnd_param = wnd_name(2:end);
    wnd_name = wnd_name{1};
  end
  if ~exist('b_ret_handler_if_possible','var')
    b_ret_handler_if_possible = false;
  end

  wnd_name = lower(wnd_name);
  switch wnd_name
    case {'rect', 'rectangular', 'dirichlet', 'daniell'}
      f_wnd = @(x) ones(size(x));
    case {'bartlett', 'triangle', 'bartlet', 'fejer'}
      f_wnd = @(x) 1 - 2*abs(x);
    case 'parzen'
      f_wnd_2 = @(x) (abs(x)<0.5).*(1-6*x.^2+6*abs(x).^3) + (abs(x)>=0.5).*2.*(1-abs(x)).^3;
      f_wnd = @(x) f_wnd_2(2*x);
    case 'rect**6'
      f_wnd_1 = @(x) (2<=x & x<3).*((3-x).^5/120) + (1<x & x<2).*(51+5*x.*(15+x.*(-42+x.*(30+x.*(-9+x)))))/120 + (0<=x & x<=1).*(33-5*(6+(x-3).*x.^2).*x.^2)/60;
      f_wnd = @(x) f_wnd_1(6*abs(x));
    case {'hanning', 'hann', 'hanning1'}
      f_wnd = @(x) 0.5+0.5*cos(2*pi*x);
    case {'hamming', 'hamming1'}
      f_wnd = @(x) 0.54+0.46*cos(2*pi*x);
    case 'blackman'
      a0=0.42; a1=0.50; a2=0.08;
      f_wnd = @(x) a0+a1*cos(2*pi*x)+a2*cos(2*pi*2*x);
    case 'exact_blackman'
      a0=7938/18608; a1=9240/18608; a2=1430/18608;
      f_wnd = @(x) a0+a1*cos(2*pi*x)+a2*cos(2*pi*2*x);
    case 'blackman-harris3'
      a0=0.42323; a1=0.49364; a2=0.07922;
      f_wnd = @(x) a0+a1*cos(2*pi*x)+a2*cos(2*pi*2*x);
    case 'blackman-harris3v2'
      a0=0.44959; a1=0.49364; a2=0.05677;
      f_wnd = @(x) a0+a1*cos(2*pi*x)+a2*cos(2*pi*2*x);
    case 'blackman-harris4'
      a0=0.35875; a1=0.48829; a2=0.14128; a3=0.01168;
      f_wnd = @(x) a0+a1*cos(2*pi*x)+a2*cos(2*pi*2*x)+a3*cos(2*pi*3*x);
    case 'blackman-harris4v2'
      a0=0.40217; a1=0.49703; a2=0.09392; a3=0.00183;
      f_wnd = @(x) a0+a1*cos(2*pi*x)+a2*cos(2*pi*2*x)+a3*cos(2*pi*3*x);
    case 'gauss2.5'
      a=2.5;
      f_wnd = @(x) exp(-0.5*(a*2*x).^2);
    case 'gauss3'
      a=3;
      f_wnd = @(x) exp(-0.5*(a*2*x).^2);
    case 'gauss5'
      a=5;
      f_wnd = @(x) exp(-0.5*(a*2*x).^2);
    case 'gauss'
      a = wnd_param{1};
      f_wnd = @(x) exp(-0.5*(a*2*x).^2);
    case {'dolph-chebyshev','chebyshev'}  %Fixme: seems this is broken
      if exist('wnd_param','var')
        a = wnd_param{1};
      else
        a = 2.0;
      end
      b = cosh(acosh(10^a)/len);
      %aacos = @(x) (abs(x)<=1.0).*(pi/2-atan(x./sqrt(1-x.*x)))...
                 %+ (abs(x)>1.0) .*log(x+sqrt(x.*x-1));
      wk = (-1).^(0:len-1).*cos(len*acos(b*cos(pi*(0:len-1)/len)))/cosh(len*acosh(b));
      f_wnd = real(ifft(wk));
    case {'kaiser', 'kaiser-bessel'}
      if exist('wnd_param','var')
        a = wnd_param{1};
      else
        a = 2.2;
      end
      f_wnd = @(x) besseli(0,pi*a*sqrt(1-(2*x).^2))/besseli(0,pi*a);
    case 'dpss'
      % Possible form of wnd_param include
      %                        nHBW  n_taper
      % {} or nonexist         2.2   1
      % {nHBW}                 nHBW  1
      % {nHBW, []}             nHBW  [2*nHBW]
      % {nHBW, n_taper}        nHBW  n_taper
      if exist('wnd_param','var') && ~isempty(wnd_param)
        nHBW = wnd_param{1};
      else
        nHBW = 2.2;
      end
      n_taper = 1;
      if length(wnd_param)>1
        n_taper = wnd_param{2};
        if isempty(n_taper)
          n_taper = floor(2*nHBW);
        end
      end
      f_wnd = dpss(len, nHBW, n_taper);
    otherwise
      error('no this type of window function');
  end

  if exist('len','var') && len>1 && length(f_wnd)==1
    f_wnd = f_wnd( ((1:len)-0.5)/len - 0.5 );
  end
  if length(f_wnd)>1 && b_ret_handler_if_possible
    f_wnd = f_wnd / sqrt(sum(f_wnd.^2));
  end

end  % funciton
