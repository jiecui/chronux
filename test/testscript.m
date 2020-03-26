%% TESTSCRIPT
% This script runs a sequence of analysis steps using the test data
% contained in directory 'data'. The data consists of a single tetrode
% recording from macaque area LIP during a memory saccade experiment. The
% data are already separated into spikes and LFPs. LFPs are contained in
% variable 'dlfp'. Spikes from two neurons are in in a struct array 'dsp'.
% Event information is in the following set of variables:
% 
% * |trialtimes| - start times of trials
% * |fixon| - fixation light comes on
% * |fixacq| - fixation acquired
% * |targon| - target light on
% * |targoff| - target light off
% * |fixoff| - fixation off
% * |saccade| - saccade
%
% Note that spikes and event times are in seconds and the sampling
% frequency for the LFP in this experiment was 1kHz. The following
% parameters are involved in this script:
% 
% * |pname| - path name on your computer where the data file LIPdata is
%                 stored.
% * |direction| - target direction to be analysed (0-7)
%
% The remaining parameters control various computations and are discussed
% in chronux.m - type Help chronux.m - or chronux manual for more
% information.
% 
% * |movingwin| - moving window size for time-frequency analysis, [winsize,
%                 winstep] in units consistent with Fs 
% * |segave| - 1: average over segment; 0: no average
% * |wintrig|
% * |winseg|
% 
% *References*
% 
% * Pesaran, B., Pezaris, J. S., Sahani, M., Mitra, P. P., & Andersen, R.
% A. (2002). Temporal structure in neuronal activity during working memory
% in macaque parietal cortex. Nature Neuroscience, 5(8), 805-811.
% * Bokil, H., Andrews, P., Kulkarni, J. E., Mehta, S., & Mitra, P. P.
% (2010). Chronux: A platform for analyzing neural signals. Journal of
% Neuroscience Methods, 192(1), 146-151. doi:10.1016/j.jneumeth.2010.06.020
% * Mitra, P. P., & Pesaran, B. (1999). Analysis of dynamic brain imaging
% data. Biophysical Journal, 76(2), 691-708.
% * Mitra, P., & Bokil, H. (2008). Observed brain dynamics. New York:
% Oxford University Press.
%
% _Modified by Richard J. Cui (richard.jie.cui@gmail.com) on Fri 04/20/2018
% 6:08:05.471 PM_

%% Set parameters and options
pname = 'data'; % path name
fname = 'LIPdata.mat'; % data file name
direction = 5; % 0-7

movingwin = [0.5 0.05]; % winsize 500 ms, winstep 5 ms
segave = 1;
wintrig = [5*movingwin(1) 5*movingwin(1)];
winseg = 2*movingwin(1);

params.Fs = 1000; % sampling frequency
params.fpass = [10 100]; % band of frequencies to be kept
params.tapers = [3 5]; % taper parameters; TW = 3
params.pad = 2; % pad factor for fft
params.err = [2 0.05];
params.trialave = 1;

%% Load data
full_name = fullfile(pname,fname);
load(full_name)

%% Compute spectrum of the first few seconds of LFP channels 1-2
NT = round(params.Fs*10*movingwin(1)); % number of points in 5 seconds data
data = dlfp(1:NT,:); % samples x channels
data1 = data(:,1:2); % only for channel 1-2
%%
% Fs = 1000 Hz; T = 5s, TW = 3, W (bandwidth) = 3/5 = .6 Hz;
% Jackknife error bar with p-value = .05
[S, f, Serr] = mtspectrumc(data1, params); % Multi-taper spectrum - continuous process
figure
ph = plot(f,10*log10(Serr(1,:)), f,10*log10(Serr(2,:)), f, 10*log10(S));
ph(1).Color = 'g';
ph(2).Color = 'g';
ph(3).Color = 'b';
ph(3).LineWidth = 2;
legend([ph(3), ph(1)], {'spectrum', 'CI'})
title('Averaged power spectrum of Channel 1-2')
xlabel('Frequency Hz')
ylabel('Spectrum dB')

%% Compute derivative of the spectrum for the same data
% Mirta, 2008, pp.203-207
phi = [0 pi/2];
[dS, f] = mtdspectrumc(data1, phi, params);
figure
plot(f, dS(1,:), f, dS(2,:));
xlabel('frequency Hz')
ylabel('Derivatives')
legend('Time','Frequency')
title('Time and frequency derivatives of spectrum of Chaanel 1-2')

%% compute coherency between channels 1-2 and  3-4
data2=data(:,3:4); % only for channel 3-4
[C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(data1,data2,params); 
figure;plot(f,C,f,Cerr(1,:),f,Cerr(2,:));xlabel('frequency'); ylabel('Coherency');

% coherency matrix of data1
[C,phi,S12,f,confC,phierr,Cerr]=cohmatrixc(data1,params);

% compute spectrogram for 1-2
[S,t,f,Serr]=mtspecgramc(data1,movingwin,params);
figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar

% compute time-frequency derivative of the spectrogram for 1-2
phi=[0 pi/2];
[dS,t,f]=mtdspecgramc(data1,movingwin,phi,params);
% pause
% figure;subplot(211);imagesc(t,f,squeeze(dS(1,:,:))'); axis xy; colorbar;
% subplot(212);imagesc(t,f,squeeze(dS(2,:,:))'); axis xy; colorbar;

% compute coherogram between 1-2 and 3-4
NT=round(movingwin(1)*Fs);
[C,phi,S12,S1,S2,t,f,confC,phierr,Cerr]=cohgramc(data1,data2,movingwin,params);
figure;imagesc(t,f,C'); axis xy; colorbar;

% compute segmented spectrum of 1 
NT=10*round(winseg*Fs);
data1=dlfp(1:NT,1);
[S,f,varS,C,Serr]=mtspectrumsegc(data1,winseg,params,segave);
figure; subplot(211);plot(f,10*log(S));
imagesc(f,f,C); axis xy;colorbar;

% compute segmented coherency between 1 and 2
NT=10*round(winseg*Fs);
data1=dlfp(1:NT,1); data2=dlfp(1:NT,2);
[C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencysegc(data1,data2,winseg,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1)); subplot(313);plot(f,10*log10(S2))

% compute spectrum of channel 1 triggered to events E
E1=targon(find(targets==direction)); 
data1=dlfp(:,1); 
[S,f,Serr]=mtspectrumtrigc(data1,E1,wintrig,params);
figure;plot(f,10*log10(S)'); axis xy; colorbar;

% compute spectrogram of channel 1 triggered to events E
E1=targon(find(targets==direction)); 
data1=dlfp(:,1); 
[S,t,f,Serr]=mtspecgramtrigc(data1,E1,wintrig,movingwin,params);
figure; imagesc(t,f,10*log10(S)'); axis xy; colorbar;

%
% Analysis - point process stored as times
%

% dsp contains 2 channels of spikes 

data=extractdatapt(dsp,[20 30]); % extract spikes occurring between 20 and 30 seconds and compute their spectrum
[S,f,R,Serr]=mtspectrumpt(data,params);
figure; plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]);

%
% Compute the derivative of the spectrum
%
phi=[0 pi/2];
[dS,f]=mtdspectrumpt(data,phi,params);
figure; plot(f,dS);

%
% Compute the derivative of the time-frequency spectrum
%
data=extractdatapt(dsp,[20 30]);
data1=data(1); data2=data(2);fscorr=[];t=[];
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypt(data1,data2,params);
figure; plot(f,C);

%
% Compute event triggered average spectrum for one of the directions
%
E1=targon(find(targets==direction));
data=dsp(1); 
[S,f,R,Serr]=mtspectrumtrigpt(data,E1,wintrig,params);
figure;plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:))); line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]); 

%
% Compute the matrix of coherencies
%
data=extractdatapt(dsp,[20 30]);
[C,phi,S12,f,zerosp,confC,phierr,Cerr]=cohmatrixpt(data,params,fscorr);

%
% Event triggered spectrogram - first way way
%
data=createdatamatpt(dsp(1),E1,wintrig);
[S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params);
figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar;

%
% Derivative of the time-frequency spectrum
%
data=createdatamatpt(dsp(1),E1,wintrig);
phi=[0 pi/2];
[dS,t,f]=mtdspecgrampt(data,movingwin,phi,params);
figure; subplot(211); imagesc(t,f,squeeze(dS(1,:,:))'); axis xy; colorbar;
subplot(212); imagesc(t,f,squeeze(dS(2,:,:))'); axis xy; colorbar;

%
% Coherogram between the two spike trains
%
data1=createdatamatpt(dsp(1),E1,wintrig);
data2=createdatamatpt(dsp(2),E1,wintrig);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgrampt(data1,data2,movingwin,params,fscorr);
figure;imagesc(t,f,C');axis xy; colorbar

%
% Event Triggered spectrogram another way
%
data=dsp(1);
[S,t,f,R,Serr]=mtspecgramtrigpt(data,E1,wintrig,movingwin,params);
imagesc(t,f,10*log10(S)'); axis xy; colorbar
%
% Segmented spectrum
%
data=extractdatapt(dsp,[20 30]);
data=data(1);
[S,f,R,varS]=mtspectrumsegpt(data,winseg,params);
plot(f,10*log10(S)); line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]);
%
% Segmented coherency
%
data=extractdatapt(dsp,[20 30]);
data1=data(1);data2=data(2);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpt(data1,data2,winseg,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))
%
% Analysis - hybrid: one continous and one point process stored as times
%
offset=1;
data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpt(data1,data2,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))


data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegcpt(data1,data2,winseg,params,segave,fscorr);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))


data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgramcpt(data1,data2,movingwin,params,fscorr);
figure; subplot(311); imagesc(t,f,C');axis xy; colorbar; subplot(312);imagesc(t,f,10*log10(S1)');axis xy; colorbar; subplot(313); imagesc(t,f,10*log10(S2)');axis xy; colorbar


%  Analysis: Binned spike counts

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds

data=dN;
[S,f,R,Serr]=mtspectrumpb(data,params);
plot(f,10*log10(S),f,10*log10(Serr(1,:)), f,10*log10(Serr(2,:))); %line(get(gca,'xlim'),[10*log10(R) 10*log10(R)]);

data=dN;
phi=[0 pi/2];
[dS,f]=mtdspectrumpb(data,phi,params);
figure; plot(f,dS); 

data=dN;
data1=data(:,1); data2=data(:,2);fscorr=[];t=[];
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypb(data1,data2,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1)); subplot(313); plot(f,10*log10(S2));

E=targon(find(targets==direction)); E=E(find(E>20 & E<450)); [dN,t]=binspikes(dsp,params.Fs,[20 500]);data=dN(:,1);
[S,f,R,Serr]=mtspectrumtrigpb(data,E,wintrig,params);
figure;plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
[C,phi,S12,f,zerosp,confC,phierr,Cerr]=cohmatrixpb(data,params,fscorr);

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
[S,t,f,R,Serr]=mtspecgrampb(data,movingwin,params);
figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar

clear R Serr Cerr S C phierr dN dS data1 data2

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
phi=[0 pi/2];
[dS,t,f]=mtdspecgrampb(data,movingwin,phi,params);

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
data1=data(:,1); data2=data(:,2);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgrampb(data1,data2,movingwin,params,fscorr);

[dN,t]=binspikes(dsp,params.Fs); % extract spikes
dN=dN(:,1);
E=E(1:6);
data=dN;
[S,t,f,R,Serr]=mtspecgramtrigpb(data,E,wintrig,movingwin,params);

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
data=data(:,1);
[S,f,R,varS]=mtspectrumsegpb(data,winseg,params,segave,fscorr);
figure; plot(f,10*log10(S)); line(get(gca,'xlim'),10*log10([R R])); 

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
data1=data(:,1);data2=data(:,2);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpb(data1,data2,winseg,params);

%
% Analysis - hybrid: one continous and one point process stored as counts
%
data1=dlfp(20000:30000,:);data1=data1(:,1:2); [dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data2=dN; data2=data2(1:end,:);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpb(data1,data2,params);

data1=data1(:,1); data2=data2(:,1);
[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegcpb(data1,data2,winseg,params,segave,fscorr);


data1=dlfp(20000:30000,:); data1=data1(:,1:2);
[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data2=dN; data2=data2(1:end,:);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phierr,Cerr]=cohgramcpb(data1,data2,movingwin,params,fscorr);
