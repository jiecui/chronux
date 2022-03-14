function data=change_row_to_column(data)
% Helper routine to transform 1d arrays into column vectors that are needed
% by other routines in Chronux
%
% Usage: data=change_row_to_column(data)
%
% Inputs:
% data -- required. If data is a matrix, it is assumed that it is of the
% form samples x channels/trials and it is returned without change. If it
% is a vector, it is transformed to a column vector. If it is a struct
% array of dimension 1, it is again returned as a column vector. If it is a
% struct array with multiple dimensions, it is returned without change
% Note that the routine only looks at the first field of a struct array.
%
% Ouputs:
% data (in the form samples x channels/trials)
%

% Revised by Richard J. Cui : Thu 06/07/2012  9:29:38.074 AM
% $Revision: 0.2 $  $Date: Sun 10/19/2014  4:49:09.149 PM $
%
% Barrow Neurological Institute
% 350 W Thomas Road
% Phoenix AZ 85013, USA
%
% Email: jie@neurocorrleate.com

% dtmp=[];
if isstruct(data)
    C = length(data);    % number of channels/trials
    for c = 1:C
        fnames_c = fieldnames(data(c));
        dtmp = data(c).(fnames_c{1});   % assume data in the first field
        data(c).(fnames_c{1}) = dtmp(:);
    end %for
else
    [N, C] = size(data);
    if N == 1 || C == 1;
        data = data(:);
    end
end

% [EOF]
