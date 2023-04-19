%__________________________________________________________________________
% Copyright (C) 2022 Todd Pataky



classdef APermuterTwoSample < spm1d.stats.nonparam.permuters.APermuter
    properties
        JA            %number of responses (Group A)
        JB            %number of responses (Group B)
        J             %total number of responses
        labelsZeros   %empty labels
    end
    
    
    methods

        function [self] = APermuterTwoSample(yA, yB)
            self.Y           = self.stack(yA, yB);
            self.Q           = size(yA, 2);
            self.JA          = size(yA, 1);
            self.JB          = size(yB, 1);
            self.J           = self.JA + self.JB;
            self.labels0     = [zeros(self.JA, 1); ones(self.JB, 1)];
            self.labelsZeros = zeros(self.J, 1);
            %calculate total number of permutations
            if isinf( factorial(self.J) )
                nPerm        = inf;
            else
                nPerm        = factorial(self.J) / ( factorial(self.JA) * factorial(self.JB)  );
            end
            self.nPermTotal  = nPerm;
        end

        function [self] = build_pdf(self, iterations)
            if iterations==-1
                ONES     = nchoosek( 1:self.J, self.JA );
                n        = self.nPermTotal;
                Z        = zeros(n, self.Q);
                for i = 1:n
                    Z(i,:) = self.get_test_stat_ones( ONES(i,:)' );
                end
            else
                % % Original method of building PDF
                % n        = iterations;
                % Z        = zeros(n, self.Q);
                % for i = 1:n
                %     ONES = randperm(self.J, self.JA);
                %     Z(i,:) = self.get_test_stat_ones( ONES' );
                % end

                % DCHartlen Modification (04-2023):
                % Use of randperm presented issues at low sample counts, as
                % there is a high risk of repeated label permutations. This
                % is because method permuted indices, not labels, with
                % repeatition. Test with small numbers of permutations
                % (<1000) found upwards of 40% repeated labels, leading to
                % issues with undefined distributions and very evident
                % non-determinism between runs.
                %
                % Revised methodology creates a list of non-repeating
                % label permutations, if and only if, the maximum number of
                % label permutations is less than 10K. If more than 10K,
                % fall bacak on random method. This method uses more
                % memory, but a 10k long array is generally managable

                n = iterations;
                Z = zeros(n, self.Q);
                % Use random method when number of permutions is v. high.
                if n > 10000
                    % Random perm method
                    for i = 1:n
                        ONES = randperm(self.J, self.JA);
                        Z(i,:) = self.get_test_stat_ones( ONES' );
                    end
                % If number of permutes is low, create a list and cycle
                % through
                else
                    % Generate random combination of "ones" 
                    permutes = nchoosek(1:self.J,self.JA); 
                    % Cycle through each permutation
                    for iPerm = 1:n
                        Z(iPerm,:) = ...
                            self.get_test_stat_ones(permutes(iPerm,:)');
                    end
                end
            end
            
            self.Z         = max(Z, [], 2);
            if self.dim == 1
                self.ZZ    = Z;
            end
        end
        
        
        function [z] = get_test_stat(self, labels)
            if self.Q==1
                [yA,yB] = deal( self.Y(labels==0,:), self.Y(labels==1,:) );
            else
                [yA,yB] = deal( self.Y(labels==0,:,:), self.Y(labels==1,:,:) );
            end
            z       = self.calc.get_test_stat(yA, yB);
        end

        function [z] = get_test_stat_ones(self, ONES)
            labels       = self.labelsZeros;
            labels(ONES) = 1;
            z            = self.get_test_stat( labels );
        end

    end
    
    methods (Access = protected)
        function [y] = stack(~, yA, yB)
            y   = [yA; yB];
        end
    end
    
    
    
end



