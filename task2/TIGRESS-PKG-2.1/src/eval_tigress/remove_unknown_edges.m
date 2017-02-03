function prediction_cleaned = remove_unknown_edges(prediction,G)

% prediction_cleaned = remove_unknown_edges(prediction,G) removes predicted
% edges that are not caused by a known transcription factor. 
% 
% Based on a function by G. Stolovitsky and R. Prill for DREAM challenge. 
% 
% See also : sparseGS, GRNInferenceEvaluation 
% 
% Anne-Claire Haury, 2012

ind=sub2ind(size(G),prediction(:,1),prediction(:,2));
isGS=find(abs(G(:)));
ind2=ismember(ind,isGS);
prediction_cleaned=prediction(ind2,:);
