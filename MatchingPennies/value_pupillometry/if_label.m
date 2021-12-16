function [if_xlabel,if_ylabel] = if_label(currPred, nPred, nback)

% determine whether to include axes label in the subplot

% input:
%   currPred: current predictor
%   nPred:    number of predictor
% output:
%   if_xlabel: whether to include x label
%   if_ylabel: whether to include y label

if nPred == 8  % for dQ and chosen Q
    if currPred-1 == 1 || currPred - 1 == 5
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred-1 >= 5
        if_xlabel = 1;
    else
        if_xlaebl = 0;
    end
    
elseif nPred == 9 % for control plot
    if mod(currPred-2,3) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred - 1 >= 7
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
    
elseif nPred == 7  % for RPE/chosenQ regression
    if mod(currPred-2, 2) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred-1 >= 4
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
    
elseif nPred == 2 & nback == 0% for pos/neg RPE
    if_ylabel = 1;
    if_xlabel = 1;

elseif nPred == 2 & nback == 2 % for choice and reward
    if mod(currPred-2,3) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred - 1 >= 7
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
elseif nPred == 14 % for future choice and reward
    if mod(currPred-2,4) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred - 1 >= 11
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
    
elseif nPred == 5  % for reward
    if mod(currPred-2,3) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred - 1 >= 3
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
elseif nPred == 10  % % for dQ and chosenQ plot
    if mod(currPred-2,3) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred - 1 >= 7
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
elseif nPred == 6  % % for dQ and chosenQ plot
    if mod(currPred-2,3) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred - 1 >= 4
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
elseif nPred == 12
    if mod(currPred-2,3) == 0
        if_ylabel = 1;
    else
        if_ylabel = 0;
    end
    if currPred-1 > 9
        if_xlabel = 1;
    else
        if_xlabel = 0;
    end
else
    if_ylabel = 1;
    if_xlabel = 0;
    
    
end
end