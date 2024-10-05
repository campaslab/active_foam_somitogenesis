%% Compute number of neighbors.

function ne=smt_faceNei(fEd,eFa,fId)

% Find near neighbor faces
neFa=unique(eFa(abs(fEd),:));
neFa=neFa(neFa~=fId & neFa~=0);

% Count number of neighbors if they are not not extracellular space
% cell.    
ne=size(neFa,1);

end