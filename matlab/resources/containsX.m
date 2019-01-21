function bool = containsX(str,exp)
%CONTAINS compares a string and expression.
%   bool = CONTAINS(str,exp) compares expressions to strings similar to the
%   contains function implemented in Matlab after 2016b. The function is
%   only necessary for versions of Matlab before 2016b.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  str, char or cell with char format strings.
%           exp, char or cell with char format strings.
%   Output: bool, boolean variable for each string.

if ischar(str) && ischar(exp)
    % Simple compare for characters
    bool = any(regexpi(str,exp));
elseif ischar(exp)
    % Compare cell with char to char
    bool = cellfun(@(x) ~isempty(x),regexpi(str,exp));
else
    % Compare cell with char to cell with char
    list = cellfun(@(x) regexpi(exp,x),str,'Un',0);
    if ischar(exp)
        bool = cellfun(@(x) ~isempty(x), list);
    else
        bool = false(size(list));
        for i = 1:length(list)
            bool(i) = any(cellfun(@(x) ~isempty(x), list{i}));
        end
    end
end
end