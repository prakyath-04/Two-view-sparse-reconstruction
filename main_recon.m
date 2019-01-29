%% reading the projection matrices and image points in different frames 
clear all;close all;clc;
load('cube_imgs.mat');load('projMatrices.mat');
proj=projMatrices;
numViews = size(proj,1);
img_pts = cell(1,numViews);
for i=1:numViews
    img_pts{1,i} = squeeze(image_pts(i,:,:))';
end
%% Reconstruction using least square approximation (SVD)
siz = size(img_pts{1},1);
pts_3d = cell(1,siz);
for j=1:siz
    A = zeros(16,4);
    for i = 1: numViews
        A(2*i-1,:) =img_pts{1,i}(j,1)*proj{i}(3,:) - proj{i}(1,:);
        A(2*i,:)=img_pts{1,i}(j,2)*proj{i}(3,:) - proj{i}(2,:);
    end
    [~,~,V] = svd(A);
    pts_3d{1,j} = V(:,end);
    pts_3d{1,j} = pts_3d{j}/pts_3d{j}(4,1); % normalizing the matrix 
    %To make the fourth co-ordinate  = 1
end
%% plotting
for i = 1:siz
    scatter3(pts_3d{1,i}(1),pts_3d{1,i}(2),pts_3d{1,i}(3),'o','y','filled'),title('reconstructed cube'),xlabel('x-axis'),ylabel('y-axis'),zlabel('z-axis');
    hold on;
end
hold off;