%% Compute face area and center

function [ar,cen]=smt_faceArea(fEd,eMd)

% Compute polygon area with straight edges.
fCrd=smt_faceVrtx(fEd,eMd);
cen=mean(fCrd);

% compute area
fCrd=[fCrd;fCrd(1,:)];
ar=0;
for i=1:size(fCrd,1)-1
    ar = ar + det(fCrd(i:i+1,:))/2;
end

end