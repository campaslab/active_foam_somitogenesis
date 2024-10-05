%% Compute centroid of given polygon
function [centroid,area]=smt_faceCentroid(fvr)

% compute determinant of each pair of vertices
fvrc=[fvr;fvr(1,:)];
det_val=zeros(size(fvr,1),1);
for vrc=1:size(det_val,1)
    det_val(vrc)=det(fvrc(vrc:vrc+1,:));
end

% compute area and perimeter. 
area=sum(det_val)/2;
centroid=sum((fvrc(1:end-1,:)+fvrc(2:end,:)).*det_val/6/area);

end