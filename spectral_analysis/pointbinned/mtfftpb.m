function [J, Msp, Nsp] = mtfftpb(data, tapers, nfft)
    % Multi-taper fourier transform - binned point process data
    %
    % Syntax:
    %   [J, Msp, Nsp] = mtfftpb(data, tapers, nfft)
    %
    % Input(s):
    %    data       - (in form samples x channels/trials or single vector)
    %    tapers     - (precalculated tapers from dpss)
    %    nfft  	    - (length of padded data)
    %
    % Output:
    %    J 	        - (fft in form frequency index x taper index x channels/trials)
    %    Msp 	    - (number of spikes per sample in each channel)
    %    Nsp 	    - (number of spikes in each channel)

    % Modified by Richard J. Cui.
    % $Revision: 0.1 $  $Date: Sat 07/09/2022 11:24:34.919 PM$
    %
    % Rocky Creek Rd NE
    % Rochester, MN 55906 USA
    %
    % Email: richard.cui@utoronto.ca

    % parse inputs
    % ------------
    if nargin < 3
        error('Need all input arguments')
    end

    % binned point process
    % --------------------
    data = change_row_to_column(data); % changes data stored as a row vector to a column vector
    [N, C] = size(data); % size of data
    Nsp = sum(data, 1); % number of spikes in each channel
    Msp = Nsp' ./ N; % mean rate for each channel

    % tapers
    % ------
    K = size(tapers, 2); % size of tapers
    tapers = tapers(:, :, ones(1, C)); % add channel indices to tapers
    H = fft(tapers, nfft, 1); % fourier transform of the tapers

    % uniform Fourier transfrom of binned point process
    % -------------------------------------------------
    data = data(:, :, ones(1, K)); % add taper indices to the data
    data = permute(data, [1 3 2]); % permute data to be of the same dimensions as H
    data_proj = data .* tapers; % multiply data by the tapers
    J = fft(data_proj, nfft, 1); % fft of projected data
    % * crucial - subtract the mean rate
    meansp = Msp(:, ones(1, K), ones(1, size(H, 1))); % add taper and frequency indices to meansp
    meansp = permute(meansp, [3, 2, 1]); % permute to get meansp with the same dimensions as H
    J = J - H .* meansp; % subtract the dc

    % [EOF]
