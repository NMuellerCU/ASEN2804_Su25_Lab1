 function [value, isTerminal, direction] = stoppingPoint(t, S)
value = S(6); % This is the variable we want to check (z)
isTerminal = 1; % tells that we want to stop the integration when
direction = 1; % the sign goes negative
end


