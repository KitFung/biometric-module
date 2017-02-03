function [G gold_complete] = sparseGS(gold_positives)

% G = sparseGS(gold_positives) returns a sparse adjacency-like matrix given
% a two-column matrix gold_positives containing verified edges. 
%
% [G gold_complete] = sparseGS(...) additionally returns a completed gold
% standard, that is where unexistent edges are specified. 
%
% Based on code by G. Stolovitsky and R. Prill for DREAM challenge. 
% 
% See also : GRNInferenceEvaluation
%
% Anne-Claire Haury, 2012



regulators = unique(gold_positives(:,1));
targets = unique(gold_positives(:,1:2));

%% lookup matrix for positive edges
A = edgelist2sparse(gold_positives(:,1:2));

%% build gold standard matrix for positives (1) AND negatives (-1)
G=A;
mask=ones(size(A));
mask(logical(A))=0;
mask(setdiff(targets,regulators),:)=0;
mask(logical(diag(ones(size(A,1),1),0)))=0;
G(logical(mask))=-1;
%% complete gold standard edge list (positives and negatives)
% edge_count = sum(sum(G~=0));
[I J] = find(G>0);	%% positives
[K L] = find(G<0);	%% negatives
gold_complete = [ I J ones(length(I),1) ; K L zeros(length(K),1) ];