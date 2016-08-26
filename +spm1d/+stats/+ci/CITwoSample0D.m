%__________________________________________________________________________
% Copyright (C) 2016 Todd Pataky
% $Id: SPM0D.m 1 2016-01-04 16:07 todd $


classdef CITwoSample0D < spm1d.stats.ci.CITwoSample & spm1d.stats.ci.CI0D
    methods
        function [self]  = CITwoSample0D(spmi, mA, mB, hstar, mu)
            self@spm1d.stats.ci.CITwoSample(spmi, mA, mB, hstar, mu)
            self@spm1d.stats.ci.CI0D(spmi)
        end
        
        function plot(self)
            self.plot_multimean()
        end
        
    end
end

