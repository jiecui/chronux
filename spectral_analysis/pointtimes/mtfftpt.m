function [J, Msp, Nsp] = mtfftpt(data, tapers, nfft, t, f, findx)
    % Multi-taper fourier transform - point process given as times
    %
    % Syntax:
    %   [J, Msp, Nsp] = mtfftpt(data, tapers, nfft, t, f, findx)
    %
    % Input(s):
    %   data        - (struct array of times with dimension channels/trials;
    %                 also takes in 1d array of spike times as a column vector)
    %   tapers      - (precalculated tapers from dpss)
    %   nfft        - (length of padded data)
    %   t           - (time points at which tapers are calculated)
    %   f           - (frequencies of evaluation)
    %   findx       - (index corresponding to frequencies f)
    %
    % Output:
    %   J           - (fft in form frequency index x taper index x channels/trials)
    %   Msp         - (number of spikes per sample in each channel)
    %   Nsp         - (number of spikes in each channel)

    % Modified by Richard J. Cui.
    % $Revision: 0.2 $  $Date: Mon 02/12/2024 09:39:13.829 AM $
    %
    % Rocky Creek Rd NE
    % Rochester, MN 55906 USA
    %
    % Email: richard.cui@utoronto.ca

    % parse inputs
    % ------------
    if nargin < 6
        error('Need all input arguments');
    end

    % point process
    % -------------
    % number of channels
    if isstruct(data)
        C = length(data);
    else
        C = 1;
    end

    % tapers
    % ------
    K = size(tapers, 2); % number of tapers
    nfreq = length(f); % number of frequencies

    if nfreq ~= length(findx)
        error('frequency information (last two arguments) inconsistent')
    end

    H = fft(tapers, nfft, 1); % fft of tapers
    H = H(findx, :); % restrict fft of tapers to required frequencies
    w = 2 * pi * f; % angular frequencies at which ft is to be evaluated

    % Fourier transform of point process
    % ----------------------------------
    Nsp = zeros(1, C);
    Msp = zeros(1, C);
    J = zeros(nfreq, K, C);

    for ch = 1:C

        if isstruct(data)
            fnames = fieldnames(data);
            % eval(['dtmp=data(ch).' fnames{1} ';'])
            dtmp = data(ch).(fnames{1});
            indx = find(dtmp >= min(t) & dtmp <= max(t));

            if ~isempty(indx)
                dtmp = dtmp(indx);
            end

        else
            dtmp = data;
            indx = find(dtmp >= min(t) & dtmp <= max(t));

            if ~isempty(indx)
                dtmp = dtmp(indx);
            end

        end

        Nsp(ch) = length(dtmp);
        Msp(ch) = Nsp(ch) / length(t);

        if Msp(ch) ~= 0
            data_proj = interp1(t', tapers, dtmp);
            exponential = exp(-1i * w' * (dtmp - t(1))');
            J(:, :, ch) = exponential * data_proj - H * Msp(ch);
        else
            J(1:nfreq, 1:K, ch) = 0;
        end

    end

    % EOF
