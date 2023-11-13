% original file: compute_descriptors_with_landmarks.m
% adapted by L Leblanc
function [fct] = compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
if nargin < 4, t = 200; end;
Basis = S.evecs(:,1:numEigs);
Ev = S.evals(1:numEigs);
A = S.A;

fct = [];

if nargin <6
    computationMethod = 'WKS';
else
    computationMethod = method;
end

% keep all the descriptors from the landmarks;
if nargin > 2
    if ~isempty(landmarks)
        if strcmp(computationMethod,'WKS')
            fct = waveKernelMap(Basis, Ev, A, t, landmarks);
        elseif strcmp(computationMethod,'HKS')
            fct = heatKernelMap(Basis, Ev, A, t, landmarks);
        else
            error('Unknown method');
        end
    end
end

% Subsample descriptors (for faster computation). More descriptors is
% usually better, but can be slower.
if nargin < 5
    num_skip = 40;
end
fct = fct(:,1:num_skip:end);

% fprintf('done computing descriptors (%d with %d landmarks)\n',size(fct,2),length(landmarks)); toc;
end