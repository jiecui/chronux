function yhat=lfsmooth(varargin)
% LFSMOOTH Provides a simple interface to locfit function
% 
% Syntax:
%   yhat = lfsmooth(x)
%   yhat = lfsmooth(x,y)
%   yhat = lfsmooth(x, 'parameter', value,...)
% 
% Input(s):
%   The same as those in locfit(); All locfit options, except evaluation
%   structures, are valid.
%
% Output(s):
%   A vector of smoothed values, at each data point.
%
% Example
%   To smooth a time series of observations,
%
% t = (1:100)';
% y = 2*sin(t/10) + normrnd(0,1,100,1);
% plot(t,y,'.')
% hold on
% plot(t,lfsmooth(t,y,'nn',0.5))
% hold off
%
% See also locfit, locfit_all.

% Last revised by Richard J. Cui. on Thu 07/05/2012 10:49:11.329 AM
%
% Visual Neuroscience Lab (Dr. Martinez-Conde)
% Barrow Neurological Institute
% 350 W Thomas Road
% Phoenix AZ 85013, USA
%
% Email: jie@neurocorrleate.com

% 
% Minimal input validation
if nargin < 1
    error( 'At least one input argument required' );
end

% fit = locfit(varargin{:},'module','simple');  % 'module' simple doesn't
% work now
x = varargin{1};
fit = locfit(x, varargin{2:end});
% yhat = fit.fit_points.fitted_values;
yhat = predict(fit, x);

% return

end % function

% [EOF]
