%% Find an edge that is located at the center of the boundary

function edId=smt_boundaryCenterEdgeFind(edge,face,gm_p)

edId=zeros(gm_p.nSmt-1,1);
edFaCs=edge{1}(:,3:4);
edFaCs(edFaCs(:,1)~=0,1)=face{1}(edFaCs(edFaCs(:,1)~=0,1));
edFaCs(edFaCs(:,2)~=0,2)=face{1}(edFaCs(edFaCs(:,2)~=0,2));
edFaCs=sort(edFaCs,2);

for smc=1:gm_p.nSmt-1
    bndEdId=find(ismember(edFaCs,[smc,smc+1],'rows'));
    bndEdLoc=mean(edge{1}(bndEdId,end-1:end),2);
    [~,idx]=sort(bndEdLoc);
    idx=idx(floor(numel(bndEdId)/2));
    edId(smc)=bndEdId(idx);
end

end


