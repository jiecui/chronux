function [S, t, f, R, Serr] = mtspecgrampt(data, movingwin, params, fscorr)
% Multi-taper time-frequency spectrogram - point process times
%
% Usage:
%
% [S, t, f, R, Serr] = mtspecgrampt(data, movingwin, params, fscorr)
%
% Input(s):
%   data        - structure array with dimension of channels/trials,
%                 i.e. data(ch).spiketimes = {spk1,...,spkM},...,;
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
%                         where p is the p value, i.e. err bar width and
%                         1-p/2 confidence interval.
%                 trialave: when 1, average over trials/channels
%                           when 0, don't average (default) (optional).
%                 minmaxtime:  max and min time for moving windows - optional,
%                              Defalult = max/min spiketimes
%   fscorr      - finite size corrections, 0 (don't use finite size corrections) or
%                 1 (use finite size corrections) - optional
%                 (available only for spikes). Defaults 0.
%
% Output(s):
%   S           - spectrogram with dimensions time x frequency x channels/trials if trialave = 0;
%                 dimensions time x frequency if trialave = 1
%   t           - times
%   f           - frequencies
%   R           - firing rate
%   Serr        - error bars - only if err(1) >= 1

% Revised by Richard J. Cui : Thu 06/07/2012  9:29:38.074 AM
% $Revision: 0.3 $  $Date: Wed 10/22/2014 11:42:38.361 AM $
%
% Barrow Neurological Institute
% 350 W Thomas Road
% Phoenix AZ 85013, USA
%
% Email: jie@neurocorrleate.com

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;

if ~isempty(params) && isfield(params, 'tapers')
    if length(params.tapers)==3 && movingwin(1)~=params.tapers(2);
        error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
    end
end % if

[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
data=change_row_to_column(data);

if isstruct(data)
    Ch = length(data);
else
    Ch = size(data, 2);     % 1d array
end
if nargin < 4 || isempty(fscorr); fscorr=0; end;
% if nargout > 4 && err(1)==0; error('Cannot compute errors with err(1)=0'); end;

if isfield(params, 'minmaxtime')
    mintime = params.minmaxtime(1);
    maxtime = params.minmaxtime(2);
else
    [mintime,maxtime]=minmaxsptimes(data);
end

tn=(mintime+movingwin(1)/2:movingwin(2):maxtime-movingwin(1)/2);
Nwin=round(Fs*movingwin(1)); % number of samples in window
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
f=getfgrid(Fs,nfft,fpass); Nf=length(f);
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
nw=length(tn);

if trialave
    S = zeros(nw,Nf);
    R = zeros(nw,1);
    if nargout==4; Serr=zeros(2,nw,Nf); end;
else
    S = zeros(nw,Nf,Ch);
    R = zeros(nw,Ch);
    if nargout==4; Serr=zeros(2,nw,Nf,Ch); end;
end

for n=1:nw
    t=linspace(tn(n)-movingwin(1)/2,tn(n)+movingwin(1)/2,Nwin);
    datawin=extractdatapt(data,[t(1) t(end)]);
    if nargout==5;
        [s,f,r,serr]=mtspectrumpt(datawin,params,fscorr,t);
        if err(1) ~= 0
            Serr(1,n,:,:)=squeeze(serr(1,:,:));
            Serr(2,n,:,:)=squeeze(serr(2,:,:));
        else
            Serr = [];
        end
    else
        [s,f,r]=mtspectrumpt(datawin,params,fscorr,t);
    end
    S(n,:,:)=s;
    R(n,:)=r;
end
t=tn;
S=squeeze(S); R=squeeze(R); if nargout==5; Serr=squeeze(Serr);end

% [EOF]
