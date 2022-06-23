function rml = rmlinesc(data, params, p, plt, f0)
% RMLINESC Removes significant sine waves from data (continuous data)
%
% Syntax: 
%   data = rmlinesc(data, params, p, plt, f0)
%
% Inputs:  
%   data        - (data in [N,C] i.e. time x channels/trials or a single
%                  vector) - required.
%   params      - (structure containing parameters) - params has the 
%                 following fields: tapers, Fs, fpass, pad
%               tapers  : precalculated tapers from dpss or in the one of
%                         the following forms: 
%                         (1) A numeric vector [TW K] where TW is the 
%                             time-bandwidth product and K is the number of
%                             tapers to be used (less than or equal to
%                             2TW-1).
%                         (2) A numeric vector [W T p] where W is the 
%                             bandwidth, T is the duration of the data and
%                             p is an integer such that 2TW-p tapers are
%                             used. In this form there is no default i.e.
%                             to specify the bandwidth, you have to specify
%                             T and p as well. Note that the units of W and
%                             T have to be consistent: if W is in Hz, T
%                             must be in seconds and vice versa. Note that
%                             these units must also be consistent with the
%                             units of params.Fs: W can be in Hz if and
%                             only if params.Fs is in Hz. The default is to
%                             use form 1 with TW=3 and K=5 
%               Fs      : (sampling frequency) -- optional. Defaults to 1.
%               fpass   : (frequency band to be used in the calculation in
%                         the form [fmin fmax])- optional. Default all
%                         frequencies between 0 and Fs/2. Note that units
%                         of Fs and fpass have to be consistent.
%               pad		: (padding factor for the FFT) - optional (can take
%                         values -1,0,1,2...). -1 corresponds to no
%                         padding, 0 corresponds to padding to the next
%                         highest power of 2 etc. e.g. For N = 500, if PAD
%                         = -1, we do not pad; if PAD = 0, we pad the FFT
%                         to 512 points, if pad=1, we pad to 1024 points
%                         etc. (0 by defaults).
%   p           - (P-value for F-test) - optional. Defaults set to 0.05.
%                 The F-distribution has two parameters 2 and 2K-2, where K
%                 is the number of DPSS (Slepian) tapers.
%   plt         - (y/n for plot and no plot respectively) - optional.
%                 Default set to 'n'. If no output assigned, plt = 'y',
%                 ignoring the user's input of plt.
%   f0          - frequencies at which you want to remove the lines - if
%                 unspecified the program uses the f statistic to determine
%                 appropriate lines.
%
% Outputs: 
%   rml         - (data with significant lines removed)

% Modified by Richard J. Cui.
% $Revision: 0.4 $  $Date: Wed 04/18/2018 11:27:50.087 AM$
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% Input parses
% =========================================================================
data=change_row_to_column(data);
C=size(data, 2);
% user_specified_pval=0;
if nargin < 2 || isempty(params)
    params=[]; 
end
if nargin < 3 || isempty(p)
    p=0.05; 
% else
%     user_specified_pval=1; 
end
if nargin < 4 || isempty(plt)
    plt='n'; 
end
if nargin < 5
    f0=[];
else
    f0 = change_row_to_column(f0);
end
if nargout == 0
    plt = 'y';
end % if
% if isempty(f0) && user_specified_pval==1
%     p=p/N;  % note that p is corrected by N
% end

% check params
% [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
% clear pad fpass err trialave
[~, ~, ~, ~, ~, ~, params] = getparams(params);

% =========================================================================
% line removal
% =========================================================================
[datafit, Amps, freqs, Fval, sig] = fitlinesc(data, params, p, 'n', f0);
datan=data-datafit;

%params.tapers=dpsschk(tapers,N,Fs); % calculate the tapers

% =========================================================================
% plot the figure
% =========================================================================
if strcmp(plt,'y') 
   [S1, f] = mtspectrumc(data, params);
   S2 = mtspectrumc(datan, params);
   for k = 1:C
       % print line information
       fprintf('The significant lines for channel %d (Amplitude, Phase, Frequency)\n', k)
       printOneChannel(Amps{k}, freqs{k})
       % plot figure
       Fval_k = Fval(:, k);
       data_k = data(:, k);
       datan_k = datan(:, k);
       S1_k = S1(:, k);
       S2_k = S2(:, k);
       fign_k = sprintf('Channel %d', k);
       plotOneChannel(f, params.Fs, sig, Fval_k, data_k, datan_k, S1_k,...
           S2_k, fign_k)
   end %for
end
if nargout == 1
    rml = datan;
end % if

end % function

% =========================================================================
% subroutines
% =========================================================================
function printOneChannel(amp, freq)

for k = 1:length(amp)
    fprintf('%8.4f, %8.4f, %8.4f\n', 2*abs(amp(k)), angle(amp(k)), freq(k))
end % for
fprintf('\n')

end % function

function plotOneChannel(f, Fs, sig, Fval, data, datan, S1, S2, figname)

N = length(data);

figure('Name', figname)
subplot(3, 2, 1)
plot((0:N-1)/Fs, data)
set(gca, 'XLim', [0, N-1]/Fs)
ax = axis;
set(gca, 'Box', 'Off')
title('Original data')

subplot(3, 2, 2)
plot(f, 10*log10(S1))
set(gca, 'Box', 'Off')
ylabel('Spectrum dB')
title('Original spectrum')

subplot(3, 2, 3)
plot((0:N-1)/Fs,datan)
axis(ax)
set(gca, 'Box', 'Off')
title('Cleaned data');

subplot(3, 2, 4)
plot(f, Fval)
line(get(gca,'xlim'), [sig sig], 'Color','r')
set(gca, 'Box', 'Off')
ylabel('F-statistic')
title('F-test')

subplot(3, 2, 5)
plot((0:N-1)/Fs, data-datan)
axis(ax)
xlabel('time s')
set(gca, 'Box', 'Off')
title('Line signals');

subplot(3, 2, 6)
plot(f,10*log10(S1),f,10*log10(S2))
set(gca, 'Box', 'Off')
legend('Original', 'Cleaned')
legend('Boxoff')
xlabel('frequency Hz')
ylabel('Spectrum dB')
title('Original and cleaned spectra')

end % function


% [EOF]
