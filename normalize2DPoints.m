function [normPts2D, T] = normalize2DPoints(pts)
% clear all;
% pts =[1 2;3 4;2 4;3 5;2 5];
mn = mean(pts);
% disp(mn);
siz = size(pts);
% for i=1:siz(1)
%     c_pts(i,1) = pts(i,1) - mn(1);
%     c_pts(i,2) = pts(i,2) - mn(2);
% end
% vr = var(c_pts);
% vr = sqrt(vr(1)*vr(1) + vr(2) * vr(2));
vr = mean(vecnorm(pts));
% d = 0;
% for i=1:siz(1)
%     d = d + sqrt(pts(i,1)*pts(i,1) + pts(i,1)*pts(i,1));
% end
% vr =d ;
scale = sqrt(2)/sqrt(vr); 
T = [scale 0 -scale*mn(1);0 scale -scale*mn(2);0 0 1];
normPts2D = zeros(siz(1),3);
for i = 1:siz(1)
    a = T*[pts(i,1) pts(i,2) 1]';
    normPts2D(i,:) = [a(1) a(2) a(3)];
end
end