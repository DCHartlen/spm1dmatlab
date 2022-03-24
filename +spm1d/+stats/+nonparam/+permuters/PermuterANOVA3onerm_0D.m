%__________________________________________________________________________
% Copyright (C) 2022 Todd Pataky



classdef PermuterANOVA3onerm_0D < spm1d.stats.nonparam.permuters.APermuterANOVA0DmultiF
    methods
        function [self] = PermuterANOVA3onerm_0D(y, A, B, C, SUBJ)
            self@spm1d.stats.nonparam.permuters.APermuterANOVA0DmultiF(y, A, B, C, SUBJ)
            self.calc           = spm1d.stats.nonparam.calculators.CalculatorANOVA3onerm(self.A, self.B, self.C, self.SUBJ);
            self.nEffects       = 7;
        end
    end
end



