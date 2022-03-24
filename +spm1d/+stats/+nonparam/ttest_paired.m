function [SnPM] = ttest_paired(yA, yB, varargin)
%__________________________________________________________________________
% Copyright (C) 2022 Todd Pataky


SnPM = spm1d.stats.nonparam.ttest(yA - yB, 0, varargin{:});
