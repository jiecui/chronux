%% TESTSCRIPT
% Script for testing Chronux package

%% Introduction
% 
% This script runs a sequence of analysis steps using the test data
% contained in directory 'data'. The data consists of a single tetrode
% recording from macaque area LIP during a memory saccade experiment. The
% data are already separated into spikes and LFPs. LFPs are contained in
% variable 'dlfp'. Spikes from two neurons are in a struct array 'dsp'.
% 
% *Event information* is in the following set of variables:
% 
% * |trialtimes|    - start times of trials
% * |fixon|         - fixation light comes on
% * |fixacq|        - fixation acquired
% * |targon|        - target light on
% * |targoff|       - target light off
% * |fixoff|        - fixation light comes off
% * |saccade|       - saccade
%
% Note that spikes and event times are in seconds and the sampling
% frequency for the LFP in this experiment was 1kHz. The following
% *parameters* are involved in this script:
% 
% * |pname|         - path name on your computer where the data file 
%                     LIPdata is stored.
% * |direction|     - target direction to be analysed (0-7)
%
% The remaining *parameters* control various computations and are discussed
% in chronux.m - type Help chronux.m or refer to chronux manual for more
% information.
% 
% * |movingwin|     - moving window size for time-frequency analysis, 
%                     [winsize, winstep] in units consistent with Fs
%                     (sampling frequency)
% * |segave|        - whether to average over segments (1: average over 
%                     segment; 0: no average)
% * |wintrig|       - window size around the trigger event,
%                     [winsize_left, winsize_right] in units consistent
%                     with Fs
% * |winseg|        - window size of a segment in units consistent with Fs

%% Global parameters and options

% clear all
clearvars

% data location
pname   = 'data';           % path name
fname   = 'LIPdata.mat';    % data file name

% analysis parameters and options
Fs          = 1000;     % sampling frequency in Hz
direction   = 5;        % 0-7

movingwin   = [0.5 0.05];   % winsize 500 ms, winstep 50 ms (overlap 90%)
segave      = 1;            % average over segments
wintrig     = 5*movingwin(1)*ones(2, 1);    % 2500 ms for both left and right
winseg      = 10*movingwin(1);  % segment size = 5 second

params.Fs   = Fs;   % sampling frequency
params.fpass    = [10 100]; % band of frequencies to be kept
params.tapers   = [3 5];    % taper parameters; TW = 3
params.pad  = 2;    % pad factor for fft to decide the number of points
params.err  = [2 0.05]; % 2 = jacknife error bar, p = 0.05
params.trialave = 1;    % overage over 2nd dimension (trials or channels)

%% Data loading
full_name = fullfile(pname,fname);
load(full_name, 'dlfp', 'dsp', 'targets', 'targon')

%% Analysis of Continous Process (Local Field Potentials, LFPs)
% 
% *Frequency domain analysis*
%
% * Spectrum estimation
% 
% _Spectrum of the first 5 seconds of LFP channels 1-2_
% 
% In this analysis we set sampling rate Fs = 1000 Hz, window size T = 5s,
% TW = 3, so that the half-bandwidth W = 3/5 = 0.6 Hz, which results in the
% frequency resolution of 2 x W = 1.2 Hz.  To get error bar we use Jackknife
% method with p-value = 0.05.

% get the data and local options
NT = round(Fs*winseg); % number of points in 5 seconds data
data = dlfp(1:NT,:); % samples x channels of 1st 5 seconds
data1 = data(:,1:2); % only for channel 1-2

% do the calculation
[S_noseg, f, Serr] = mtspectrumc(data1, params);

% plot the results
figure
ph = plot(f,10*log10(Serr(1,:)), f, 10*log10(Serr(2,:)), f, 10*log10(S_noseg));
set(ph(1), 'Color', 'g', 'LineWidth', .5);
set(ph(2), 'Color', 'g', 'LineWidth', .5);
set(ph(3), 'Color', 'b', 'LineWidth', 2);
legend([ph(3), ph(1)], {'spectrum', 'C-Interval'})
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Averaged power spectrum of Channel 1&2')

%%
% _Derivative of the spectrum ([4] Mirta, 2008, pp.203-207)_

% get the data and local options
NT = round(Fs*winseg); % number of points in 5 seconds data
data = dlfp(1:NT,:); % samples x channels of 1st 5 seconds
data1 = data(:,1:2); % only for channel 1-2
phi = [0 pi/2];

% do the calculation
[dS, f] = mtdspectrumc(data1, phi, params);

% present the results
figure
plot(f, dS(1,:), f, dS(2,:));
xlabel('Frequency (Hz)')
ylabel('Derivatives')
legend('Time','Frequency')
title('Time and frequency derivatives of spectrum of Channels 1-2')

%% 
% _Spectrum of channel 1 triggered by events E_

% get the data and local options
E = targon(targets == direction); 
data1 = dlfp(:,1); 

% do the calculation
[S_noseg, f, Serr] = mtspectrumtrigc(data1, E, wintrig, params);

% plot the results
ph = plot(f,10*log10(Serr(1,:)), f, 10*log10(Serr(2,:)), f, 10*log10(S_noseg));
set(ph(1), 'Color', 'g', 'LineWidth', .5);
set(ph(2), 'Color', 'g', 'LineWidth', .5);
set(ph(3), 'Color', 'b', 'LineWidth', 2);
legend([ph(3), ph(1)], {'spectrum', 'C-Interval'})
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Event triggered power spectrum of Channel 1')

%% 
% _Segmented spectrum of channel 1&2_
%
% We segment a 50-second signal into 10 5-second segments and then average
% the spectrum.  

% get the data and local options
NT = 10*round(winseg*Fs); % number of samples in 50 second window segment
data1 = dlfp(1:NT, 1);

% do the calculation
[S, f, varS, C, Serr] = mtspectrumsegc(data1, winseg, params, segave);

% plot the results
figure
ph = plot(f, pow2db(S_noseg), f, pow2db(S));
set(ph(2), 'LineWidth', 2)
grid on
legend('non-seg', 'seg-smoothed')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Comparison of averaged power spectrum of Channel 1-2')

figure
subplot(2, 1, 1)
ph = plot(f, pow2db(S), f, pow2db(Serr(1, :)), f, pow2db(Serr(2, :)));
set(ph(1), 'Color', '#D95319', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
ylabel('Power Spectrum (dB)')
title('Segment averaged power spectrum of Channel 1-2')

subplot(2, 1, 2)
plot(f, varS)
xlabel('Frequency (Hz)')
ylabel('Variance')
title('Variance of the log spectrum')

figure
imagesc(f,f,C)
axis xy
colormap('jet'), colorbar('northoutside')
xlabel('Frequency (Hz)')
ylabel('Frequency (Hz)')
title('Covarance matrix of log spectrum')

%% 
% * Coherency estimation
% 
% _Coherency matrix between channels 1 and 2_ 

% get the data and local options
NT = round(winseg*Fs); % number of samples in 5 second window segment
data1 = dlfp(1:NT, 1:2);

% do the calculation
[C_noseg, phi, S12, f, confC, phistd, Cerr] = cohmatrixc(data1, params);

% display the results
fprintf('Confidence level between channel 1 and 2 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC(1, 2))

% plot the results
figure
subplot(2, 1, 1)
plot(f, pow2db(abs(S12(:, 1, 2))))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Cross spectrum between Channel 1 & 2')

subplot(2, 1, 2)
plot(f, angle(S12(:, 1, 2)))
xlabel('Frequency (Hz)')
ylabel('Phase')

figure
subplot(2, 1, 1)
ph = plot(f, C_noseg(:, 1, 2), f, Cerr(1, :, 1, 2), f, Cerr(2, :, 1, 2));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Coherence between Channel 1 and 2')

subplot(2, 1, 2)
ph = plot(f, phi(:, 1, 2), f, phi(:, 1, 2)-phistd(1, :, 1, 2).',...
    f, phi(:, 1, 2)+phistd(2, :, 1, 2).');
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Phase')

%%
% _Segmented coherency between channels 1 and channels 2_

% get the data and local options
NT = 10*round(winseg*Fs); % number of samples in 50 second window segment
data1 = dlfp(1:NT, 1);
data2 = dlfp(1:NT, 2);

% do the calculation
[C, phi, S12, S1, S2, f, confC, phistd, Cerr] = coherencysegc(data1, data2,...
    winseg, params);

% display the results
fprintf('Segment averaged confidence level between channel 1 and 2 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC)

% plot the results
figure
ph = plot(f, C_noseg(:, 1, 2), f, C);
ylim([.9 1])
set(ph(2), 'LineWidth', 2)
legend('no-seg', 'seg-smoothed')
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Comparison segmented average coherence between channels 1 and 2')

figure
subplot(2, 1, 1)
ph = plot(f, C, f, Cerr(1, :), f, Cerr(2, :));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Segmented average coherence between channels 1 and 2')

subplot(2, 1, 2)
ph = plot(f, phi, f, phi-phistd, f, phi+phistd);
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Phase')

figure
subplot(3, 1, 1)
plot(f, pow2db(abs(S12)))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Segment averaged cross power spectrum between channels 1 and 2')

subplot(3, 1, 2)
plot(f,10*log10(S1))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Segment averaged power spectrum of channels 1')

subplot(3, 1, 3)
plot(f,10*log10(S2))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Segment averaged power spectrum of channels 2')

%% 
% _Coherency between channels 1-2 and 3-4_
% 
% Here the two channels in each group are treated as two trials

% get the data and local options
NT      = round(Fs*winseg); % number of points in 5 seconds data
data    = dlfp(1:NT,:); % samples x channels of 1st 5 seconds
data1   = data(:, 1:2); % only for channel 1-2
data2   = data(:, 3:4); % only for channel 3-4

% do the calculation
[C, phi, S12, S1, S2, f, confC, phistd, Cerr] = coherencyc(data1, data2, params); 

% display the results
fprintf('Confidence level between channel 1&2 and 3&4 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC)

% plot the results
figure
subplot(2, 1, 1)
plot(f, pow2db(abs(S12)))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Cross spectrum between Channel 1-2 and 3-4')

subplot(2, 1, 2)
plot(f, angle(S12))
xlabel('Frequency (Hz)')
ylabel('Phase')

figure
subplot(2, 1, 1)
plot(f, pow2db(S1))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Spectrum of Channel 1-2')

subplot(2, 1, 2)
plot(f, pow2db(S2))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Spectrum of Channel 3-4')

figure
subplot(2, 1, 1)
ph = plot(f, C, f, Cerr(1,:), f, Cerr(2,:));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('frequency (Hz)')
ylabel('Coherency');
title('Coherence between Channels 1-2 and 3-4')

subplot(2, 1, 2)
ph = plot(f, phi, f, phi-phistd, f, phi+phistd);
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Phase')

%% Analysis of Point Process (Spike Trains)
% 
% *Frequency domain analysis*
% 
% * Spectrum estimation
% 
% _Spectrum of spike trains_

% get the data and local options
% dsp contains 2 channels of spikes;
% extract spikes occurring between 20 and 30 seconds and compute their spectrum
cell12 = extractdatapt(dsp, [20 30]); 

[S, f, R, Serr] = mtspectrumpt(cell12, params);

% plot the results
figure
ph = plot(f, 10*log10(S), f, 10*log10(Serr(1,:)), f, 10*log10(Serr(2,:)));
set(ph(1), 'Color', '#D95319', 'LineWidth', 2);
set(ph(2), 'Color', 'y', 'LineWidth', .5);
set(ph(3), 'Color', 'y', 'LineWidth', .5);
line(get(gca,'xlim'),[10*log10(R) 10*log10(R)], 'LineWidth', 2);
xlabel('Frequency (Hz)')
ylabel('Spike Spectrum (dB)')
title('Spike train spectrum of cell 1 within 20-30 sec')

%%
% *Compute the derivative of the spectrum of spike trains*

phi=[0 pi/2];
[dS, f] = mtdspectrumpt(cell12, phi, params);

% plot the results
figure
plot(f, dS);
xlabel('Frequency (Hz)')
ylabel('Derivatives')
legend('Time','Frequency')
title('Time and frequency derivatives of spectrum of cells 1 and 2')

%%
% *Compute segmented spectrum of spike train of cell 1*

data = extractdatapt(dsp, [20 30]);
data=data(1);
[S, f, R, varS] = mtspectrumsegpt(data, winseg, params);

% display the results
fprintf('varS = %0.2f\n', varS)

% plot the results
figure
plot(f,10*log10(S))
line(get(gca,'xlim'),[10*log10(R) 10*log10(R)])
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Segment averaged power spectrum of cell 1')

%%
% *Compute the coherency between two spike trains*

cell1 = cell12(1); 
cell2 = cell12(2);
fscorr = [];
t = [];
[C, phi, S12, S1, S2, f, zerosp, confC, phistd, Cerr] = coherencypt(cell1,...
    cell2, params, fscorr, t);

% display the results
fprintf('zerosp is %d.\n', zerosp)
fprintf('Confidence level between cell 1 and 2 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC)

% plot the results
figure
subplot(2, 1, 1)
ph = plot(f, C, f, Cerr(1, :), f, Cerr(2, :));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Coherence between spike trains of cell 1 and 2')

subplot(2, 1, 2)
ph = plot(f, unwrap(phi), f, unwrap(phi)-unwrap(phistd),...
    f, unwrap(phi)+unwrap(phistd));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Phase')

figure
subplot(3, 1, 1)
plot(f, pow2db(abs(S12)))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Cross power spectrum between cell 1 and 2')

subplot(3, 1, 2)
plot(f,10*log10(S1))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Power spectrum of cell 1')

subplot(3, 1, 3)
plot(f,10*log10(S2))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Power spectrum of cell 2')

%%
% *Compute event triggered average spectrum for spike trains*

cell1_whole = dsp(1);
[S, f, R, Serr] = mtspectrumtrigpt(cell1_whole, E1, wintrig, params);

% plot the results
figure
ph = plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));
set(ph(1), 'Color', '#D95319', 'LineWidth', 2);
set(ph(2), 'Color', 'y', 'LineWidth', .5);
set(ph(3), 'Color', 'y', 'LineWidth', .5);
line(get(gca,'xlim'),[10*log10(R) 10*log10(R)], 'LineWidth', 2)
xlabel('Frequency (Hz)')
ylabel('Spike Spectrum (dB)')
title('Event triggered spike train spectrum of cell 1')

%%
% *Compute the matrix of coherencies between cell 1 and 2*

fscorr=[];
[C, phi, S12, f, zerosp, confC, phistd, Cerr] = cohmatrixpt(cell12,...
    params, fscorr);

% display the results
fprintf('zerosp is %d.\n', zerosp)
fprintf('Confidence level between cell 1 and 2 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC(1, 2))

% plot the results
figure
subplot(2, 1, 1)
plot(f, pow2db(abs(S12(:, 1, 2))))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Cross spectrum between cells 1 and 2')

subplot(2, 1, 2)
plot(f, unwrap(angle(S12(:, 1, 2))))
xlabel('Frequency (Hz)')
ylabel('Phase')

figure
subplot(2, 1, 1)
ph = plot(f, C(:, 1, 2), f, Cerr(1, :, 1, 2), f, Cerr(2, :, 1, 2));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Coherence between cells 1 and 2')

subplot(2, 1, 2)
ph = plot(f, unwrap(phi(:, 1, 2)), f, unwrap(phi(:, 1, 2))-unwrap(phistd(1, :, 1, 2).'),...
    f, unwrap(phi(:, 1, 2))+unwrap(phistd(2, :, 1, 2).'));
set(ph(1), 'Color', '#0072BD', 'LineWidth', 2);
set(ph(2), 'Color', '#D95319', 'LineWidth', .5);
set(ph(3), 'Color', '#D95319', 'LineWidth', .5);
xlabel('Frequency (Hz)')
ylabel('Phase')

%%
% *Segmented coherency between cell 1 and 2*

data=extractdatapt(dsp,[20 30]);
data1=data(1);
data2=data(2);
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr] = ...
    coherencysegpt(data1,data2,winseg,params);

figure
subplot(311)
plot(f,C);subplot(312)
plot(f,10*log10(S1))
subplot(313)
plot(f,10*log10(S2))

%% Analysis of Hybrid Process (LFP and Spike Trains)
%
% *LFP and spike train analysis of cell 1*

offset=1;
data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(data1,data2,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))


data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencysegcpt(data1,data2,winseg,params,segave,fscorr);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))


data1=dlfp(20000:30000,1); data2=extractdatapt(dsp,[20 30],offset);data2=data2(1).times;
[C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr]=cohgramcpt(data1,data2,movingwin,params,fscorr);
figure; subplot(311); imagesc(t,f,C');axis xy; colorbar; subplot(312);imagesc(t,f,10*log10(S1)');axis xy; colorbar; subplot(313); imagesc(t,f,10*log10(S2)');axis xy; colorbar



%% 
% *Time-frequency analysis of LFP*
%
% * Compute spectrogram of Channels 1&2
[S, t, f, Serr] = mtspecgramc(data1, movingwin, params);

% plot the results
figure
subplot(3, 1, 1)
imagesc(t, f, 10*log10(S)')
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Spectrogram of channels 1&2')

subplot(3, 1, 2)
imagesc(t, f, 10*log10(squeeze(Serr(1, :, :)))')
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Lower bound of error of channels 1&2')

subplot(3, 1, 3)
imagesc(t, f, 10*log10(squeeze(Serr(2, :, :)))')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Upper bound of error of channels 1&2')

%%
% *Compute time-frequency derivative of the spectrogram for channels 1&2*

phi=[0 pi/2];
[dS, t, f] = mtdspecgramc(data1, movingwin, phi, params);

% plot the results
figure
subplot(211)
imagesc(t, f, squeeze(dS(1,:,:))')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(sprintf('Time-frequency derivative at angle %0.2f rad', phi(1)))

subplot(212)
imagesc(t, f, squeeze(dS(2,:,:))')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(sprintf('Time-frequency derivative at angle %0.2f rad', phi(2)))

%% 
% *Compute spectrogram of channel 1 triggered by events E*

[S, t, f, Serr] = mtspecgramtrigc(data1_e, E1, wintrig, movingwin, params);

% plot the results
figure
subplot(3, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(S)')
axis xy
colormap('jet'), colorbar
hold on
plot([0 0], ylim, 'w', 'LineWidth', 2)
ylabel('Frequency (Hz)')
title('Event triggered spectrogram of channels 1')

subplot(3, 1, 2)
imagesc(t-wintrig(1), f, 10*log10(squeeze(Serr(1, :, :)))')
axis xy
colormap('jet'), colorbar
hold on
plot([0 0], ylim, 'w', 'LineWidth', 2)
ylabel('Frequency (Hz)')
title('Lower bound of error of channels 1')

subplot(3, 1, 3)
imagesc(t-wintrig(1), f, 10*log10(squeeze(Serr(2, :, :)))')
axis xy
colormap('jet'), colorbar
hold on
plot([0 0], ylim, 'w', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Upper bound of error of channels 1')

%% 
% *Compute coherogram between channels 1&2 and channles 3&4*

[C, phi, S12, S1, S2, t, f, confC, phistd, Cerr] = cohgramc(data1, data2,...
    movingwin, params);

% display the results
fprintf('Confidence level of coherogram between channel 1&2 and 3&4 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC)

% plot the results
figure
subplot(2, 1, 1)
imagesc(t, f, C')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram magnitude between channels 1&2 and chnnels 3&4')

subplot(2, 1, 2)
imagesc(t, f, phi')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram phase between channels 1&2 and chnnels 3&4')

figure
subplot(2, 1, 1)
imagesc(t, f, squeeze(Cerr(1, :, :))')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram magnitude lower bound between channels 1&2 and chnnels 3&4')

subplot(2, 1, 2)
imagesc(t, f, squeeze(Cerr(2, :, :))')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram magnitude upper bound between channels 1&2 and chnnels 3&4')

figure
subplot(2, 1, 1)
imagesc(t, f, (phi-phistd)')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram phase lower bound between channels 1&2 and chnnels 3&4')

subplot(2, 1, 2)
imagesc(t, f, (phi+phistd)')
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram phase upper bound between channels 1&2 and chnnels 3&4')

figure
subplot(2, 1, 1)
imagesc(t, f, 10*log10(S1)')
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Spectrogram of channels 1&2')

subplot(2, 1, 2)
imagesc(t, f, 10*log10(S2)')
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Spectrogram of channels 3&4')

figure
subplot(2, 1, 1)
imagesc(t, f, 10*log10(abs(S12))')
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Cross spectrogram magnitude of channels 1&2 and channels 3&4')

subplot(2, 1, 2)
imagesc(t, f, angle(S12)')
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Cross spectrogram phase of channels 1&2 and channels 3&4')

%% 
% *Time-frequency analysis of spike trains*
%
% * Event triggered spectrogram of spike train - first way

cell1_trig = createdatamatpt(dsp(1), E1, wintrig);
[S, t, f, R, Serr] = mtspecgrampt(cell1_trig, movingwin, params);

% plot the results
figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(S)')
axis xy
colormap('jet'), colorbar('northoutside')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title('Event triggered spectrogram of cell 1 spike train')

subplot(2, 1, 2)
plot(t-wintrig(1), pow2db(R))
line([0 0], ylim, 'Color', 'k', 'LineWidth', 2)
axis tight
xlabel('Time around event onset (s)')
ylabel('Log spike rate (dB)')
title('Event triggered spike rate of cell 1')

figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(squeeze(Serr(1, :, :)))')
axis xy
colormap('jet'), colorbar('northoutside')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title('Event triggered spectrogram lower bound of cell 1 spike train')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, 10*log10(squeeze(Serr(2, :, :)))')
axis xy
colormap('jet'), colorbar('northoutside')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title('Event triggered spectrogram upper bound of cell 1 spike train')

%%
% *Event Triggered spectrogram of spike train - second way*

[S, t, f, R, Serr] = mtspecgramtrigpt(cell1_whole, E1, wintrig,...
    movingwin, params);

% plot the results
figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(S)')
axis xy
colormap('jet'), colorbar('northoutside')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title('Event triggered spectrogram of cell 1 spike train')

subplot(2, 1, 2)
plot(t-wintrig(1), pow2db(R))
line([0 0], ylim, 'Color', 'k', 'LineWidth', 2)
axis tight
xlabel('Time around event onset (s)')
ylabel('Log spike rate (dB)')
title('Event triggered spike rate of cell 1')

figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(squeeze(Serr(1, :, :)))')
axis xy
colormap('jet'), colorbar('northoutside')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title('Event triggered spectrogram lower bound of cell 1 spike train')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, 10*log10(squeeze(Serr(2, :, :)))')
axis xy
colormap('jet'), colorbar('northoutside')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title('Event triggered spectrogram upper bound of cell 1 spike train')

%%
% *Derivative of the time-frequency spectrum of spike train*

phi=[0 pi/2];
[dS, t, f] = mtdspecgrampt(cell1_trig, movingwin, phi, params);

% plot the results
figure
subplot(211)
imagesc(t-wintrig(1), f, squeeze(dS(1,:,:))')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title(sprintf('Dirivative spectrogram of cell 1 spike train at angle %0.1fº',...
   rad2deg(phi(1))))

subplot(212)
imagesc(t-wintrig(2), f, squeeze(dS(2,:,:))')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time around event onset (s)')
ylabel('Frquency (Hz)')
title(sprintf('Dirivative spectrogram of cell 1 spike train at angle %0.1fº',...
   rad2deg(phi(2))))

%%
% *Coherogram between the two spike trains*

cell2_trig = createdatamatpt(dsp(2),E1,wintrig);
[C, phi, S12, S1, S2, t, f, zerosp, confC, phistd, Cerr] = ...
    cohgrampt(cell1_trig, cell2_trig, movingwin, params, fscorr);

% display the results
fprintf('zerosp = %d\n', zerosp)
fprintf('Confidence level of coherogram between cell 1 and 2 at %0.2f%% leve is %0.2f.\n',...
    1-params.err(2), confC)

% plot the results
figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, C')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram magnitude between cell 1 and 2')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, phi')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram phase between cell 1 and 2')

figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, squeeze(Cerr(1, :, :))')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram magnitude lower bound between cell 1 and 2')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, squeeze(Cerr(2, :, :))')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram magnitude upper bound between cell 1 and 2')

figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, (phi-phistd)')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram phase lower bound between cell 1 and 2')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, (phi+phistd)')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherogram phase upper bound between cell 1 and 2')

figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(S1)')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Spectrogram of cell 1')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, 10*log10(S2)')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Spectrogram of cell 2')

figure
subplot(2, 1, 1)
imagesc(t-wintrig(1), f, 10*log10(abs(S12))')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Cross spectrogram magnitude of cell 1 and 2')

subplot(2, 1, 2)
imagesc(t-wintrig(1), f, angle(S12)')
line([0 0], get(gca, 'ylim'), 'Color', 'w', 'LineWidth', 2)
axis xy
colormap('jet'), colorbar
ylabel('Frequency (Hz)')
title('Cross spectrogram phase of cell 1 and 2')

%%
% Analysis: Binned spike counts

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
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencypb(data1,data2,params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1)); subplot(313); plot(f,10*log10(S2));

E=targon(find(targets==direction)); E=E(find(E>20 & E<450)); [dN,t]=binspikes(dsp,params.Fs,[20 500]);data=dN(:,1);
[S,f,R,Serr]=mtspectrumtrigpb(data,E,wintrig,params);
figure;plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));

[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data=dN;
[C,phi,S12,f,zerosp,confC,phistd,Cerr]=cohmatrixpb(data,params,fscorr);

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
[C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr]=cohgrampb(data1,data2,movingwin,params,fscorr);

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
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencysegpb(data1,data2,winseg,params);

%%
% Analysis - hybrid: one continous and one point process stored as counts
%
data1=dlfp(20000:30000,:);data1=data1(:,1:2); [dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data2=dN; data2=data2(1:end,:);
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpb(data1,data2,params);

data1=data1(:,1); data2=data2(:,1);
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencysegcpb(data1,data2,winseg,params,segave,fscorr);


data1=dlfp(20000:30000,:); data1=data1(:,1:2);
[dN,t]=binspikes(dsp,params.Fs,[20 30]); % extract spikes occurring between 20 and 30 seconds
data2=dN; data2=data2(1:end,:);
[C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr]=cohgramcpb(data1,data2,movingwin,params,fscorr);

%% References
% 
% # Pesaran, B., Pezaris, J. S., Sahani, M., Mitra, P. P., & Andersen, R.
%   A. (2002). Temporal structure in neuronal activity during working
%   memory in macaque parietal cortex. Nature Neuroscience, 5(8), 805-811.
% # Bokil, H., Andrews, P., Kulkarni, J. E., Mehta, S., & Mitra, P. P.
%   (2010). Chronux: A platform for analyzing neural signals. Journal of
%   Neuroscience Methods, 192(1), 146-151.
%   doi:10.1016/j.jneumeth.2010.06.020
% # Mitra, P. P., & Pesaran, B. (1999). Analysis of dynamic brain imaging
%   data. Biophysical Journal, 76(2), 691-708.
% # Mitra, P., & Bokil, H. (2008). Observed brain dynamics. New York:
%   Oxford University Press.
% 
% _Modified by Richard J. Cui (richard.jie.cui@gmail.com) on Thu 03/26/2020
% 10:49:56.888 AM_