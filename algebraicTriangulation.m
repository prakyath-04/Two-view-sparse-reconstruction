function [pts_3d] = algebraicTriangulation(pts_a_2d, pts_b_2D, ProjMat_a, ProjMat_b)
    Pa1=ProjMat_a(1,:);
    Pa2=ProjMat_a(2,:);
    Pa3=ProjMat_a(3,:);
    Pb1=ProjMat_b(1,:);
    Pb2=ProjMat_b(2,:); 
    Pb3=ProjMat_b(3,:);
    for i=1:size(pts_a_2d,2)
        x=pts_a_2d(1,i); y=pts_a_2d(2,i); xa=pts_b_2D(1,i); ya=pts_b_2D(2,i);
        A= [(x*Pa3-Pa1);  (y*Pa3-Pa2); (xa*Pb3-Pb1);  (ya*Pb3-Pb2)];
        [~,~,V] = svd(A);
        pts_3d(:,i) = V(:,end);
    end
end

% function [pts3D] = algebraicTriangulation(pts2D_1, pts2D_2, ProjMat_1, ProjMat_2)
% 
% ProjMat_1=[ProjMat_1(1,1), ProjMat_1(1,2), ProjMat_1(1,4); ProjMat_1(2,1), ProjMat_1(2,2), ProjMat_1(2,4); ProjMat_1(3,1), ProjMat_1(3,2), ProjMat_1(3,4)];
% ProjMat_2=[ProjMat_2(1,1), ProjMat_2(1,2), ProjMat_2(1,4); ProjMat_2(2,1), ProjMat_2(2,2), ProjMat_2(2,4); ProjMat_2(3,1), ProjMat_2(3,2), ProjMat_2(3,4)];
% disp(size(ProjMat_1));
% disp(size(pts2D_1));
% pts3D=inv(ProjMat_1)*pts2D_1;
% temp=inv(ProjMat_2)*pts2D_2;
% pts3D=[pts3D temp];
% %pts3D=transpose(pts3D);
% end