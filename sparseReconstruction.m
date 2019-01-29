clear all;close all;
%% initializations 
a = imread('../images/img1.png');
b = imread('../images/img2.png');
a_i = rgb2gray(a);
b_i = rgb2gray(b);
pts1 = detectSURFFeatures(a_i,'MetricThreshold',10);
pts2 = detectSURFFeatures(b_i,'MetricThreshold',10);
[ft1,vpts1] = extractFeatures(a_i,pts1);
[ft2,vpts2] = extractFeatures(b_i,pts2);
idxpairs = matchFeatures(ft1,ft2) ;
matchedpts1 = vpts1(idxpairs(:,1));
matchedpts2 = vpts2(idxpairs(:,2));
%% montage plot 
a = [matchedpts1.Location];
b = [matchedpts2.Location];
ax = axes;
showMatchedFeatures(a_i,b_i,matchedpts1,matchedpts2,'montage','Parent',ax);
title(ax, 'feature point matches');
legend(ax, 'Matched points 1','Matched points 2');
% showMatchedFeatures(a,b,matchedpts1,matchedpts2);
% legend('matched points 1','matched points 2');
%% estimating fundamental and essential matrix 
[norm_F,a1,b1] = estimateFundamentalMatrixRANSAC(a,b);
% [norm_F,ab] = estimateFundamentalMatrix(matchedpts1.Location,matchedpts2.Location,'NumTrials',300);
% disp(size(ab));
% for i= 1:size(ab)
%     if(ab(i)==1)
%         a1(i,1) = a(i,1);a1(i,2) = a(i,2);
%         b1(i,1) = b(i,1);b1(i,2) = b(i,2);
%     end
% end
[u,d,v] = svd(norm_F);
new_d = diag([d(1,1) d(2,2) , 0]);
F = u * new_d * v' ;
% F = F ./ F(3,3);
fprintf("The Fundamental Matrix is:\n");
disp((F./F(3,3)));

K =[558.7087 0.0000 310.3210 ;...
    0.0000 558.2827 240.2395 ;...
    0.0000 0.0000 1.0000 ] ;

E = K'*F*K;
[u,d,v] = svd(E);
new_d = diag([(d(1,1)+d(2,2))/2, (d(1,1)+d(2,2))/2 , 0]);
E = u * new_d * v' ;
%% finding rotation and translation matrix 
a_h = ones(size(a1,1),3);
a_h(:,1:2) = a1;
a_h = a_h';
b_h = ones(size(a1,1),3);
b_h(:,1:2) = b1;
b_h = b_h';
[R, t] = decomposeEssentialMatrix(E, a1, b1, K);
fprintf("The Rotation Matrix is:\n");
disp(R);
fprintf("The Translation Matrix is:\n");
disp(t);

T2 = K * [R t];
T1 = K * [eye(3,3) [0 0 0]'];
figure;
%% plots
% plotCameraFrustum(T2,'r', .1);hold on;
% plotCameraFrustum(T1,'b', .1);hold off;

pts_3d = algebraicTriangulation(a_h,b_h, T1, T2);
pts_3d = pts_3d ./ pts_3d(4,:);
scatter3(pts_3d(1,:),pts_3d(2,:),pts_3d(3,:),'filled');
% pts_3d = algebraicTriangulation(mhpts1, transpose(mhpts2), ProjMat_1, ProjMat_2);
% surf(pts_3d)