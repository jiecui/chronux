function [Fval, A, f, sig, sd] = ftestc(data, params, p, plt)
    % FTESTC computes the F-statistic for sinusoidal wave in locally-white noise (continuous data).
    %
    % Syntax:
    %   [Fval, A, f, sig, sd] = ftestc(data, params, p, plt)
    %
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
    %
    % Outputs:
    %   Fval        - (F-statistic in frequency x channels/trials format)
    %   A           - (Line amplitude for X in frequency x channels/trials form)
    %	  f           - (frequencies of evaluation)
    %   sig         - (F distribution (1-p)% confidence level)
    %   sd          - (standard deviation of the amplitude A)

    % Modified by Richard J. Cui.
    % $Revision: 0.2 $  $Date: Fri 04/20/2018  3:33:08.023 PM$
    %
    % 3236 E Chandler Blvd Unit 2036
    % Phoenix, AZ 85048, USA
    %
    % Email: richard.cui@utoronto.ca

    % =========================================================================
    % Parse inputs
    % =========================================================================
    if nargin < 1
        error('Need data');
    end

    if nargin < 2 || isempty(params)
        params = [];
    end

    if nargin < 3 || isempty(p)
        % p=0.05/N;
        p = 0.05;
    end

    if nargin < 4 || isempty(plt)
        plt = 'n';
    end

    if nargout == 0
        plt = 'y';
    end % if

    % =========================================================================
    % Line F-test
    % =========================================================================
    [tapers, pad, Fs, fpass, ~, ~, params] = getparams(params);
    data = change_row_to_column(data);
    [N, C] = size(data);

    tapers = dpsschk(tapers, N, Fs); % calculate the tapers
    [N, K] = size(tapers);
    nfft = max(2^(nextpow2(N) + pad), N); % number of points in fft
    [f, findx] = getfgrid(Fs, nfft, fpass); % frequency grid to be returned
    % errorchk = 0; % set error checking to default (no errors calculated)
    % if nargout <= 3 % if called with 4 output arguments, activate error checking
    %     errorchk = 0;
    % else
    %     errorchk = 1;
    % end
    Kodd = 1:2:K;
    Keven = 2:2:K;
    J = mtfftc(data, tapers, nfft, Fs); % tapered fft of data - f x K x C
    Jp = J(findx, Kodd, :); % drop the even ffts and restrict fft to specified frequency grid - f x K x C
    tapers = tapers(:, :, ones(1, C)); % add channel indices to the tapers - t x K x C
    H0 = squeeze(sum(tapers(:, Kodd, :), 1)); % calculate sum of tapers for even prolates - K x C

    if C == 1
        H0 = H0';
    end

    Nf = length(findx); % number of frequencies
    H0 = H0(:, :, ones(1, Nf)); % add frequency indices to H0 - K x C x f
    H0 = permute(H0, [3 1 2]); % permute H0 to get dimensions to match those of Jp - f x K x C
    H0sq = sum(H0 .* H0, 2); % sum of squares of H0^2 across taper indices - f x C
    JpH0 = sum(Jp .* squeeze(H0), 2); % sum of the product of Jp and H0 across taper indices - f x C
    A = squeeze(JpH0 ./ H0sq); % amplitudes for all frequencies and channels
    Kp = size(Jp, 2); % number of even prolates
    Ap = A(:, :, ones(1, Kp)); % add the taper index to C
    Ap = permute(Ap, [1 3 2]); % permute indices to match those of H0
    Jhat = Ap .* H0; % fitted value for the fft

    num = (K - 1) .* (abs(A).^2) .* squeeze(H0sq); %numerator for F-statistic
    den = squeeze(sum(abs(Jp - Jhat).^2, 2) + sum(abs(J(findx, Keven, :)).^2, 2)); % denominator for F-statistic
    Fval = num ./ den; % F-statisitic
    sig = finv(1 - p, 2, 2 * K - 2); % F-distribution (1-p) percential
    var = den ./ (K * squeeze(H0sq)); % variance of amplitude
    sd = sqrt(var); % standard deviation of amplitude

    % plot
    if nargout == 0 || strcmp(plt, 'y')
        [S, f] = mtspectrumc(data, params);

        for k = 1:C
            figure('Name', sprintf('Channel %d', k))
            subplot(211)
            plot(f, 10 * log10(S(:, k)));
            set(gca, 'Box', 'Off')
            xlabel('frequency Hz');
            ylabel('Spectrum dB');

            subplot(212)
            plot(f, Fval(:, k));
            line(get(gca, 'xlim'), [sig sig], 'Color', 'r');
            set(gca, 'Box', 'Off')
            xlabel('frequency Hz');
            ylabel('F ratio');
        end % for

    end

    A = A * Fs;

end % function

% [EOF]
