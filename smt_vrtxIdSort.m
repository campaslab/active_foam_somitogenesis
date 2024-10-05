%% Sort vertex ID in counter-clockwise direction. 
% This function is used before periodic boundary condition is applied.

function vSrt=vrIdSrtFn(vLst,vCrd)    

% assign face vertex coordinates.
vCrd=vCrd(vLst,:); 

% calculate centroid position and angle about x-axis.
cCrd=[mean(vCrd(:,1)),mean(vCrd(:,2))];
angV=atan2(vCrd(:,2)-cCrd(2),vCrd(:,1)-cCrd(1));

% sort vertices in the order of angle value.    
[~,idc]=sort(angV);
vSrt=vLst(idc);    

end