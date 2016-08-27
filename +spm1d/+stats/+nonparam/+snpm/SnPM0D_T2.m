%__________________________________________________________________________
% Copyright (C) 2016 Todd Pataky
% $Id: SPM0D.m 1 2016-01-04 16:07 todd $


classdef SnPM0D_T2 < spm1d.stats.nonparam.snpm.ASnPM0D
    methods
        function [self] = SnPM0D_T2(z, perm, varargin)
            self@spm1d.stats.nonparam.snpm.ASnPM0D(z, perm, varargin{:})
            self.STAT = 'T2';
        end
    end
end