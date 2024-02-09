function [datafit, Amps, freqs, Fval, sig] = fitlinesc(data, params, p, plt, f0)
    % FITLINESC fits significant sine waves to data (continuous data).
    %
    % Syntax:
    %   [datafit, Amps, freqs, Fval, sig] = fitlinesc(data, params, p, plt, f0)
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
    %   datafit     - (linear superposition of fitted sinusoidal waves) the
    %                 signal model is 2*|Aamps|*cos(2*pi*freqs*t+angle(Amps)).
    %   Amps        - (amplitudes at significant frequencies)
    %   freqs       - (significant frequencies)
    %   Fval        - (Fstatistic at all frequencies)
    %   sig         - (significance level for F distribution p value of p)

    % Modified by Richard J. Cui.
    % $Revision: 0.3 $  $Date: Sat 06/04/2022 12:38:44.215 AM$
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =========================================================================
    % Parse inputs
    % =========================================================================
    data = change_row_to_column(data);
    [N, C] = size(data);

    if nargin < 2 || isempty(params)
        params = [];
    end

    % [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
    % clear pad fpass err trialave;
    [tapers, ~, Fs, ~, ~, ~, params] = getparams(params);
    % if nargin < 3 || isempty(p);p=0.05/N;end;
    if nargin < 3 || isempty(p)
        p = 0.05;
    end

    if nargin < 4 || isempty(plt)
        plt = 'n';
    end

    if nargin < 5
        f0 = [];
    else
        f0 = change_row_to_column(f0);
    end

    if nargout == 0
        plt = 'y';
    end % if

    % =========================================================================
    % Line fit
    % =========================================================================
    params.tapers = dpsschk(tapers, N, Fs); % calculate the tapers
    [Fval, A, f, sig] = ftestc(data, params, p, 'n');
    dfit1ch = @(Amp, freq) 2 * real(exp(1i * 2 * pi * (0:N - 1)' * freq / Fs) * Amp);

    freqs = cell(1, C);
    Amps = cell(1, C);
    datafit = zeros(size(data));

    if isempty(f0) % if line frequency not specified
        fmax = findpeaks(Fval, sig); % chronux function

        for ch = 1:C
            fsig = f(fmax(ch).loc);
            freqs{ch} = fsig;
            Amps{ch} = A(fmax(ch).loc, ch);
            % datafit(:, ch) = exp(1i*2*pi*(0:N-1)'*fsig/Fs)*Amps{ch}...
            %     +exp(-1i*2*pi*(0:N-1)'*fsig/Fs)*conj(Amps{ch});
            datafit(:, ch) = fitplt1ch(dfit1ch, Amps{ch}, freqs{ch}, plt);
        end

    else % line frequency specified
        % indx = zeros( length(f0) );
        indx = zeros(1, length(f0));

        for n = 1:length(f0)
            [~, indx(n)] = min(abs(f - f0(n))); % find the closest frequency
        end

        fsig = f(indx);

        for ch = 1:C
            freqs{ch} = fsig;
            Amps{ch} = A(indx, ch);
            % datafit(:,ch)=exp(1i*2*pi*(0:N-1)'*fsig/Fs)*A(indx,ch)...
            %     +exp(-1i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(indx,ch));
            datafit(:, ch) = fitplt1ch(dfit1ch, Amps{ch}, freqs{ch}, plt);
        end

    end

end % function

% =========================================================================
% subroutines
% =========================================================================
function datafit = fitplt1ch(dfit1ch, A, f, plt)

    datafit = dfit1ch(A, f);
    % plot
    if plt == 'y'
        Nf = length(fsig);
        fprintf('The significant lines for channel %d (Amplitude, Phase, Frequency)\n', ch)

        for k = 1:Nf
            fprintf('%8.4f, %8.4f, %8.4f\n', 2 * abs(A(k)), angle(A(k)), f(k))
        end % for

    end % if

end % function

% [EOF]
