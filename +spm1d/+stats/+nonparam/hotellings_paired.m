function [SnPM] = hotellings_paired(yA, yB, varargin)
%__________________________________________________________________________
% Copyright (C) 2016 Todd Pataky
% $Id: ttest.m 1 2016-01-04 16:07 todd $

SnPM  = spm1d.stats.nonparam.hotellings(yA - yB);

