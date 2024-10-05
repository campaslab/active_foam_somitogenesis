%% Compute an individual edge length.
% If an edge is straight, calculate the distance between two end points.

function ln=smt_edgeLen(eCrd)

% compute edge length        
ln=zeros(size(eCrd,1)-1,1);
for ii=1:size(eCrd,1)-1
    ln(ii)=norm(eCrd(ii+1,:)-eCrd(ii,:));
end

end