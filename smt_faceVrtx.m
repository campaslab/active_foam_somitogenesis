%% Obtain face-vertex coordinate, including intermediate vertices

function fCrd=smt_faceVrtx(fEd,eMd)

% for each edge, use intermediate vertex coordinates to get face-vertex
% coordinates
fCrd=cell(size(fEd,2),1);
for ii=1:size(fEd,2)
    if fEd(ii)>0
        fCrd{ii}=eMd{fEd(ii)}(1:end-1,:);
    else
        fCrd{ii}=flipud(eMd{-fEd(ii)}(2:end,:));
    end
end    
fCrd=cell2mat(fCrd);

end