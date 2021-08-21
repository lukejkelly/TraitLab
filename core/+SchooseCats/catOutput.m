function [catStruct] = catOutput(catVector)
    % Converting to a struct so that we can reuse Bcats functions
    catStruct = struct('i', catVector(1), 'j', catVector(2), 'k', catVector(3));
end
