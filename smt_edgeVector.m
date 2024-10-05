%% Compute normal and tangential vector 
function [tv,nv,eln]=smt_edgeVector(eVs,eVf,fCn)

% compute tangential vector
tv=eVf-eVs;
md=(eVf+eVs)/2;

% normalize tangential vector
eln=norm(tv);
tv=tv/eln;

% compute normal vector, from face 1 to face 2.
nv=[-tv(2),tv(1)];
fv=md-fCn;
if dot(nv,fv)<0
    nv=-nv;
end

end