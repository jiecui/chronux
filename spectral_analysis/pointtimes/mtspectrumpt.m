function [S, f, R, Serr] = mtspectrumpt(data, params, fscorr, t)
% Multi-taper spectrum for point process times
%
% Syntax:
%   [S, f, R, Serr] = mtspectrumpt(data, params, fscorr, t)
% 
% Input(s): 
%   data        - structure array of spike times with fields of channels/trials,
%                 i.e. data.channel1 = {spiketime1,...,spiketimeM},...,;
%                 also accepts 1d array of spike times (required)
%   movingwin   - length of moving window and step size = [window, winstep]
%   params      - structure array of MT parameters with fields of tapers, pad, 
%                 Fs, fpass, err, trialave (optional)
%                 tapers: precalculated tapers from dpss or in the one of the following
%                         forms:
%                         (1) A numeric vector [TW K] where TW is the
%                             time-bandwidth product and K is the number of
%                             tapers to be used (less than or equal to
%                             2TW-1).
%                         (2) A numeric vector [W T p] where W is the
%                             bandwidth, T is the duration of the data and p
%                             is an integer such that 2TW-p tapers are used. In
%                             this form there is no default i.e. to specify
%                             the bandwidth, you have to specify T and p as
%                             well. Note that the units of W and T have to be
%                             consistent: if W is in Hz, T must be in seconds
%                             and vice versa. Note that these units must also
%                             be consistent with the units of params.Fs: W can
%                             be in Hz if and only if params.Fs is in Hz.
%                         The default is to use form 1 with TW=3 and K=5. 
%                         Note that T has to be equal to movingwin(1).
%	              pad:    padding factor for the FFT (can take values -1,0,1,2...) (optional).
%                         -1 corresponds to no padding, 0 corresponds to
%                         padding to the next highest power of 2 etc. e.g. 
%                         For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%                         to 512 points, if pad=1, we pad to 1024 points etc.
%			      	      Defaults to 0.
%                 Fs:     sampling frequency; default is 1 (optional).
%                 fpass:  frequency band to be used in the calculation in the form
%                         [fmin fmax]; Default all frequencies between 0 and Fs/2
%                         (optional).
%                 err:    error calculation 
%                         [1 p] - Theoretical error bars; 
%                         [2 p] - Jackknife error bars
%                         [0 p] or 0 - no error bars) - optional. Default 0.
%                         where p is err bar width and 1-p/2 confidence
%                         interval.
%                 trialave: when 1, average over trials/channels 
%                           when 0, don't average (default) (optional). 
%                 minmaxtime:  max and min time for moving windows - optional,
%                              Defalult = max/min spiketimes
%   fscorr      - finite size corrections, 0 (don't use finite size corrections) or
%                 1 (use finite size corrections) - optional
%                 (available only for spikes). Defaults 0.
%   t           - time grid over which the tapers are to be calculated:
%                 this argument is useful when calling the spectrum
%                 calculation routine from a moving window spectrogram
%                 calculation routine. If left empty, the spike times
%                 are used to define the grid.
% 
% Output(s):
%   S           - power spectrum with dimensions frequency x channels/trials if
%                 trialave = 0, or with dimension frequency if trialave = 1
%   f           - frequencies
%   R           - rate
%   Serr        - spectrum error bars, only if params.err(1) >= 1

% Revised by Richard J. Cui : Thu 06/07/2012  9:29:38.074 AM
% $Revision: 0.3 $  $Date: Wed 10/22/2014 10:41:09.288 AM $
%
% Barrow Neurological Institute
% 350 W Thomas Road
% Phoenix AZ 85013, USA
%
% Email: jie@neurocorrleate.com

% parse inputs
if nargin < 1
    error('Need data'); 
end

if nargin < 2
    params = []; 
end

[tapers, pad, Fs, fpass, err, trialave] = getparams(params);
data=change_row_to_column(data);

% if nargout > 3 && err(1)==0; error('cannot compute error bars with err(1)=0; change params and run again'); end;
if nargin < 3 || isempty(fscorr); fscorr=0;end;
if nargin < 4 || isempty(t);
   [mintime,maxtime]=minmaxsptimes(data);
   dt=1/Fs; % sampling time
   t=mintime-dt:dt:maxtime+dt; % time grid for prolates
end

N = length(t); % number of points in grid for dpss
nfft = max(2^(nextpow2(N) + pad), N); % number of points in fft of prolates
[f, findx] = getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
tapers = dpsschk(tapers, N, Fs); % check tapers

% FFT
[J, Msp, Nsp] = mtfftpt(data, tapers, nfft, t, f, findx); % mt fft for point process times

% get power stpectrum
S = squeeze(mean(conj(J).*J,2));    % average across tapers
if trialave == 1
    S = squeeze(mean(S,2));     % average across channel/trials
    Msp = mean(Msp);
end

R = Msp * Fs;
if nargout == 4
    if err(1) ~= 0
        if fscorr == 1
            Serr = specerr(S,J,err,trialave,Nsp);
        else
            Serr = specerr(S,J,err,trialave);
        end
    else
        Serr = [];
    end
end

% [EOF]