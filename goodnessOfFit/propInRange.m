function [prop] = propInRange(X, cons)
% Proportion of entries in X in range [cons(1), cons(2)]
  prop = mean( (cons(1) <= X) & (X <= cons(2)) );
end
