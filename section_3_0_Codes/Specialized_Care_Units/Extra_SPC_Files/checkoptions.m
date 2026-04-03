
%CHECKOPIONS check option structure

% $Revision: 1.1 $  $Date: 2005/02/11 16:41:56 $

function [yn,bad] = checkoptions(options, goodopt)
    % Minimal version: returns true only if all fields in options are valid
    fields = fieldnames(options);
    bad = setdiff(fields, goodopt);
    yn = isempty(bad);
end