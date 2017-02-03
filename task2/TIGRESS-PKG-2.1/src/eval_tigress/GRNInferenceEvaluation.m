function [TPR FPR PREC REC L AUROC AUPR P_AUROC P_AUPR] = GRNInferenceEvaluation(gold_positives, prediction_raw, pdf_aupr, pdf_auroc)


% DREAM5-like accuracy evaluation. Returns accuracy statistics given a
% prediction matrix prediction_raw, as obtained by, e.g. predict_network.m.
% This function and dependent function are based on DREAM5 evaluation
% functions written by G. Stolovitsky and R. Prill.
%
% See also predict_network, sparseGS, remove_unknown_edges, probability, edgelist2sparse 
% 
% Anne-Claire Haury, 2012


%% preprocess gold standard to get graphs G and H
idx = gold_positives(:,3);
gold_positives = gold_positives(logical(idx),:);		%% make sure they are really positives
G = sparseGS(gold_positives);
H = G>0;	%% just the positives

%% total, positive, negative
P = full(sum(sum(G>0)));
N = full(sum(sum(G<0)));
T = P + N;

%% preprocess prediction so only contains edges in the gold standard
prediction = remove_unknown_edges(prediction_raw,G);
L = size(prediction,1);

%% create the "discovery" vector, the order with which positves are discovered
I = prediction(:,1);
J = prediction(:,2);
discovery=sub2ind(size(H),I,J);
isGS=find(H(:));
discovery=ismember(discovery,isGS);


%% number of positives discovered in the prediction list of length L
TPL = sum(discovery);

%% p = the number of remaining positives / number of remaining list entries.
%% p is the uniform probablility of picking a positive by guessing.
if L < T
	p = (P - TPL) / (T - L);
else
	p = 0;	%% they were all already found
end

%% p is also the fractional number of positives discovered per guess,
random_positive_discovery = repmat(p, T-L, 1);
random_negative_discovery = repmat(1-p, T-L, 1);

%% "complete" the discovery vectors
positive_discovery = [  discovery ; random_positive_discovery ];
negative_discovery = [ ~discovery ; random_negative_discovery ];

%% true positives (false positives) at depth k
TPk = cumsum(positive_discovery);
FPk = cumsum(negative_discovery);

K = (1:T)';
%% sanity check 
if ( (P ~= round(TPk(end))) || (N ~= round(FPk(end))) )
	disp('ERROR. There is a problem with the completion of the prediction list.')
end

%% finishing touch
TPk(end) = round(TPk(end));
FPk(end) = round(FPk(end));
%% metrics
TPR = TPk / P;
FPR = FPk / N;
REC = TPR;  %% same thing
PREC = TPk ./ K;

%% faster built-in integration function
AUROC = trapz(FPR,TPR);
AUPR = trapz(REC,PREC) / (1-1/P);	%% normalized by max possible value.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calcaulate p-values for these metrics
if nargin==4
    X = pdf_aupr.X;
    Y = pdf_aupr.Y;
    P_AUPR = probability(X, Y, AUPR);

    X = pdf_auroc.X;
    Y = pdf_auroc.Y;
    P_AUROC = probability(X, Y, AUROC);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1); clf
%pn = 2;
%pm = 2;
%pj = 0;
%
%pj = pj + 1; subplot(pn,pm,pj)
%plot(FPR(1:L),TPR(1:L),'r')
%hold on
%plot(FPR(L+1:end),TPR(L+1:end),'b')
%hold off
%xlabel('FPR')
%ylabel('TPR')
%hold on
%plot([0,1],[0,1],'k--')  %% diagonal
%hold off
%axis([0,1,0,1])
%
%pj = pj + 1; subplot(pn,pm,pj)
%plot(REC(1:L),PREC(1:L),'r')
%hold on
%plot(REC(L+1:end),PREC(L+1:end),'b')
%hold off
%xlabel('Recall')
%ylabel('Precision')
%axis([0,1,0,1])
