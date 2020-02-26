

%%%% this script is to check the density. I am kind of just playing with
%%%% it. 

%%%% I found a code in stack overflow to make a heatmap base on the
%%%% Euclidean disntances. 


%%%% this is the original code:

a = rand(1000,3);       % Create random matrix, use your data here
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.2;              % Tolerance for (squared) distance to count as "nearby" 
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
    dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance 
    n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just 
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;


%%%% this is one I play with the f20 data

load('All_More_BrainReg2.mat');
load('f20_cleaned_idxs.mat');

a = ROI_temp2.f20(idx_rsq_test_f20short_cleaned,:);       % f20 filterd ROIs
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 1000;              % Tolerance for (squared) distance to count as "nearby" 
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
    dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance 
    n(ii) = nnz(dists < tol); % Count number of points within tolerance
    %n(ii) = nnz(dists < tol2); % get the distances
    mean_dists(ii)=mean(dists);
end

%tol2=mean(mean_dists)/2;

% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], '.');
grid on;

