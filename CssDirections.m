classdef CssDirections < int32
    %CSSDIRECTIONS 
    % class with the line directions definition
    enumeration
        forward (0)
        reverse (1)
    end
    methods
        function out = str(self)
            out = {};
            for i = 1:length(self)
                if self(i) == 0
                    out = [out; {'Forward'}];
                else
                    out = [out; {'Reverse'}];
                end
            end
        end
        
        function out = Flip(self)
            out = [];
            for i = 1:length(self)
                if self(i) == 0
                    out = [out; CssDirections.reverse];
                else
                    out = [out; CssDirections.forward];
                end
            end
            if ~isempty(out)
                out = reshape(out(:), size(self));
            end
        end
    end
    
    methods(Static)
        function out = ToDirection(value)
            if value == 0
                out = CssDirections.forward;
            else
                out = CssDirections.reverse;
            end
        end
    end
end

