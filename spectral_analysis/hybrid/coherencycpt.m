function [C, phi, S12, S1, S2, f, zerosp, confC, phistd, Cerr] = coherencycpt(data1, data2, params, fscorr, t)
    % Multi-taper coherency, cross-spectrum and individual spectra - continuous data and point process as times
    %
    % Syntax:
    %   [C, phi, S12, S1, S2, f, zerosp, confC, phistd, Cerr] = coherencycpt(data1, data2, params, fscorr, t)
    %
    % Input(s):
    %   data1       - (continuous data in time x trials form) -- required
    %   data2       - (structure array of spike times with dimension trials;
    %                  also accepts 1d array of spike times) -- required
    %   params      - structure with fields: tapers, pad, Fs, fpass, err, trialave
    %                 -- optional
    %                 tapers : precalculated tapers from dpss or in the one
    %                          of the following forms:
    %                          (1) A numeric vector [TW K] where TW is the
    %                              time-bandwidth product and K is the number of
    %                              tapers to be used (less than or equal to
    %                              2TW-1).
    %                          (2) A numeric vector [W T p] where W is the
    %                              bandwidth, T is the duration of the data and p
    %                              is an integer such that 2TW-p tapers are used. In
    %                              this form there is no default i.e. to specify
    %                              the bandwidth, you have to specify T and p as
    %                              well. Note that the units of W and T have to be
    %                              consistent: if W is in Hz, T must be in seconds
    %                              and vice versa. Note that these units must also
    %                              be consistent with the units of params.Fs: W can
    %                              be in Hz if and only if params.Fs is in Hz.
    %                              The default is to use form 1 with TW=3 and K=5
    %                 Fs     : (sampling frequency) - optional. Default 1.
    %                 fpass  : (frequency band to be used in the calculation
    %                          in the form [fmin fmax])- optional. Default
    %                          all frequencies between 0 and Fs/2. Note that
    %                          units of Fs and fpass have to be consistent.
    %	                pad	  :	(padding factor for the FFT) - optional (can take
    %                          values -1,0,1,2...). -1 corresponds to no
    %                          padding, 0 corresponds to padding to the next
    %                          highest power of 2 etc. e.g. For N = 500, if
    %                          FPAD = -1, we do not pad; if PAD = 0, we pad
    %                          the FFT to 512 points, if pad=1, we pad to 1024
    %                          points etc (Defaults to 0).
    %                 err    : (error calculation [1 p] - Theoretical error bars;
    %                          [2 p] - Jackknife error bars; [0 p] or 0 - no
    %                          error bars) - optional (Default 0).
    %                 trialave: (average over trials when 1, don't average
    %                          when 0) - optional. Default 0
    %   fscorr      - (finite size corrections), 0 (don't use finite size
    %                 corrections) or 1 (use finite size corrections)
    %                 - optional (available only for spikes) (Defaults 0).
    %   t           - (time grid over which the tapers are to be calculated:
    %                  this argument is useful when calling the spectrum
    %                  calculation routine from a moving window spectrogram
    %                  calculation routine). If left empty, the spike times
    %                  are used to define the grid.
    %
    % Output(s):
    %   C           - (magnitude of coherency - frequencies x trials if trialave=0;
    %                 dimension frequencies if trialave=1)
    %   phi         - (phase of coherency - frequencies x trials if trialave=0;
    %                 dimension frequencies if trialave=1)
    %   S12         - (cross spectrum -  frequencies x trials if trialave=0;
    %                 dimension frequencies if trialave=1)
    %   S1          - (spectrum 1 - frequencies x trials if trialave=0;
    %                 dimension frequencies if trialave=1)
    %   S2          - (spectrum 2 - frequencies x trials if trialave=0;
    %                 dimension frequencies if trialave=1)
    %   f           - (frequencies)
    %   zerosp      - (1 for trials where no spikes were found, 0 otherwise)
    %   confC       - (confidence level for C at 1-p %) - only for err(1)>=1
    %   phistd      - theoretical/jackknife (depending on err(1)=1/err(1)=2)
    %                 standard deviation for phi.  Note that phi + 2 phistd
    %                 and phi - 2 phistd will give 95% confidence bands for
    %                 phi - only for err(1)>=1
    %   Cerr        - (Jackknife error bars for C - use only for Jackknife
    %                 -  err(1)=2)

    % Modified by Richard J. Cui.
    % $Revision: 0.1 $  $Date: Thu 08/04/2022 12:17:26.769 AM$
    %
    % Rocky Creek Rd NE
    % Rochester, MN 55906 USA
    %
    % Email: richard.cui@utoronto.ca

    % parse inputs
    % ------------
    if nargin < 2
        error('Need data1 and data2')
    end

    if nargin < 3
        params = [];
    end

    [tapers, pad, Fs, fpass, err, trialave, params] = getparams(params);
    clear params

    if nargin < 4 || isempty(fscorr)
        fscorr = 0;
    end

    if nargin < 5 || isempty(t)
        [N, C] = size(data1);
        dt = 1 / Fs;
        t = 0:dt:(N - 1) * dt; % time grid for prolates
    end

    if nargout > 7 && err(1) == 0
        %   Errors computed only if err(1) is non-zero. Need to change params and run again.
        error('When errors are desired, err(1) has to be non-zero.')
    end

    if nargout > 9 && err(1) ~= 2
        error('Cerr computed only for Jackknife. Correct inputs and run again')
    end

    % compute coherencey
    % ------------------

    % * frequency grid
    [N, Ch] = check_consistency(data1, data2, 1);
    zerosp = zeros(1, Ch); % intialize the zerosp variable
    N = length(t); % number of points in grid for dpss
    nfft = max(2^(nextpow2(N) + pad), N); % number of points in fft of prolates
    [f, findx] = getfgrid(Fs, nfft, fpass);

    % * tapers
    tapers = dpsschk(tapers, N, Fs); % check tapers

    % * calculation
    J1 = mtfftc(data1, tapers, nfft, Fs); % fourier transform of continuous data
    J1 = J1(findx, :, :); % restrict to required frequencies
    [J2, Msp2, Nsp2] = mtfftpt(data2, tapers, nfft, t, f, findx); % fourier transform of discrete data
    zerosp(Nsp2 == 0) = 1; % set zerosp to 1 for trials where no spikes were found
    S12 = squeeze(mean(conj(J1) .* J2, 2)); % cross spectrum
    S1 = squeeze(mean(conj(J1) .* J1, 2)); % spectrum data 1
    S2 = squeeze(mean(conj(J2) .* J2, 2)); % spectrum data 2
    if trialave; S12 = squeeze(mean(S12, 2)); S1 = squeeze(mean(S1, 2)); S2 = squeeze(mean(S2, 2)); end;
    C12 = S12 ./ sqrt(S1 .* S2);
    C = abs(C12);
    phi = angle(C12);

    if nargout == 10;

        if fscorr == 1;
            [confC, phistd, Cerr] = coherr(C, J1, J2, err, trialave, [], Nsp2);
        else
            [confC, phistd, Cerr] = coherr(C, J1, J2, err, trialave);
        end;

    elseif nargout == 9;

        if fscorr == 1;
            [confC, phistd] = coherr(C, J1, J2, err, trialave, [], Nsp2);
        else
            [confC, phistd] = coherr(C, J1, J2, err, trialave);
        end;

    end;

    clear Msp2

end % function

% [EOF]
