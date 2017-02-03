function P = probability(X,Y,x)

% Computes p-value of x under distribution of X : P(X>=x). 
% (X,Y) defines the probability distribution function.
% Initial code by Stolovitsky et al. Optimized by Anne-Claire Haury, 2012.
%
% See also GRNInferenceEvaluation

dx = X(2)-X(1);
P = sum(Y(X>=x)*dx);
