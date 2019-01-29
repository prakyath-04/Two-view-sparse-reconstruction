function [norm_F,inliers_a,inliers_b] = estimateFundamentalMatrixRANSAC(pts1,pts2)
mxInliers = 0;

norm_F = zeros(3,3);
x_a = [pts1 ones(size(pts1,1),1)];
x_b = [pts2 ones(size(pts2,1),1)];

for i = 1:10000
    ind = randi(size(pts1,1), [8,1]); %% randomly selecting 8 features from all the features extracted 
    [n_pts1 ,t1]= normalize2DPoints(pts1(ind,:));
    [n_pts2 ,t2]= normalize2DPoints(pts2(ind,:)); %% normalizing 2d points 
    %% finding fundamental matrix for the randomly selected feature points
    A = zeros(8,9);
    for j = 1:8
        A(j,:) = [n_pts1(j,1)*n_pts2(j,1),n_pts1(j,1)*n_pts2(j,2),n_pts1(j,1),...
            n_pts1(j,2)*n_pts2(j,1),n_pts1(j,2)*n_pts2(j,2),n_pts1(j,2)...
            ,n_pts2(j,1),n_pts2(j,2),1];
    end
    [~,~,v] = svd(A);
    f= v(:,end);
    f = reshape(f,[3 3]);
    
    [u,d,v] = svd(f);
    new_d = diag([d(1,1) d(2,2) , 0]);
    F = u * new_d * v' ;
    
    F = t2' * F * t1;
    %% finding the error in the estimated F matrix and checking it with the 
    %% previouly estimated F matrix 
    err = sum((x_b .* (F * x_a')'),2);
    curr_Inliers = size( find(abs(err) <= 0.05) , 1);
    if (curr_Inliers > mxInliers)
        norm_F = F;
        mxInliers = curr_Inliers;
    end
end
%% sorting according to the error values  
err = sum((x_b .* (norm_F * x_a')'),2);
[~,I]  = sort(abs(err),'ascend');
inliers_a = pts1(I(1:340),:);
inliers_b = pts2(I(1:340),:);
% disp(t1);
% disp(t2);
% disp(n_pts1);
end