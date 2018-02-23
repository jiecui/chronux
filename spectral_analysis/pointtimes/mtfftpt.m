function [J,Msp,Nsp]=mtfftpt(data,tapers,nfft,t,f,findx)
% Multi-taper fourier transform for point process given as times
%
% Usage:
% [J,Msp,Nsp]=mtfftpt (data,tapers,nfft,t,f,findx) - all arguments required
% Input: 
%       data        (struct array of times with dimension channels/trials; 
%                   also takes in 1d array of spike times as a column vector) 
%       tapers      (precalculated tapers from dpss) 
%       nfft        (length of padded data) 
%       t           (time points at which tapers are calculated)
%       f           (frequencies of evaluation)
%       findx       (index corresponding to frequencies f) 
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
%       Msp (number of spikes per sample in each channel)
%       Nsp (number of spikes in each channel)

% Revised by Richard J. Cui : Thu 06/07/2012 10:37:09.142 AM
% $Revision: 0.1 $  $Date: Thu 06/07/2012 10:37:09.142 AM $
%
% Visual Neuroscience Lab (Dr. Martinez-Conde)
% Barrow Neurological Institute
% 350 W Thomas Road
% Phoenix AZ 85013, USA
%
% Email: jie@neurocorrleate.com

% Revised by Richard J. Cui : Thu 06/07/2012  9:29:38.074 AM
% $Revision: 0.2 $  $Date: Sun 10/19/2014  5:06:05.048 PM $
%
% Barrow Neurological Institute
% 350 W Thomas Road
% Phoenix AZ 85013, USA
%
% Email: jie@neurocorrleate.com

if nargin < 6; error('Need all input arguments'); end
if isstruct(data)
    C = length(data); 
else
    C = size(data, 2); 
end % number of channels/trials

K=size(tapers,2); % number of tapers
nfreq=length(f); % number of frequencies
if nfreq~=length(findx); error('frequency information (last two arguments) inconsistent'); end;
H=fft(tapers,nfft,1);  % fft of tapers
H=H(findx,:); % restrict fft of tapers to required frequencies
w=2*pi*f; % angular frequencies at which ft is to be evaluated
Nsp=zeros(1,C); Msp=zeros(1,C);

J = zeros(nfreq, K, C);     % frequency index x taper index x channel/trials
for ch = 1:C
  if isstruct(data);
     % eval(['dtmp=data(ch).' fnames{1} ';'])
     fnames_ch = fieldnames(data(ch));
     dtmp = data(ch).(fnames_ch{1});
     % indx = find(dtmp>=min(t) & dtmp <= max(t));
     % if ~isempty(indx)
     %     dtmp = dtmp(indx);
     % end
  else
     dtmp = data(:, ch);
     % indx = find(dtmp >= min(t)&dtmp<=max(t));
     % if ~isempty(indx); dtmp=dtmp(indx);
     % end
  end
  indx = dtmp >= min(t) & dtmp <= max(t);
  dtmp = dtmp(indx);
  
  Nsp(ch)=length(dtmp);
  Msp(ch)=Nsp(ch)/length(t);
  if Msp(ch)~=0;
      data_proj=interp1(t',tapers,dtmp);
      exponential=exp(-1i*w'*(dtmp-t(1))');
      J(:,:,ch)=exponential*data_proj-H*Msp(ch);
  else
      J(1:nfreq,1:K,ch)=0;
  end
end

% [EOF]
