%% Compute face-edge id from face-vrtx id
function fed=smt_faceEdgeId(edge,face,gm_p,rg,fId)

fvr=face{2}{fId};
eVr=edge{1}(:,rg.ei(1):rg.ef(1));

% general caess.
if size(fvr,2)>2
    fvrc=[fvr,fvr(1)];    
    fed=zeros(1,size(fvr,2));

    for jj=1:size(fvrc,2)-1
        % find edge id based on edge-vrtx id
        if fvrc(jj)<fvrc(jj+1)
            eid=find(ismember(eVr,fvrc(jj:jj+1),'rows')==1);
            edSgn=1;
        else
            eid=find(ismember(eVr,[fvrc(jj+1),fvrc(jj)],'rows')==1);
            edSgn=-1;
        end
        
        % eid will be unique if an edge is a part of cell with 1 neighbors.
        % If eid is not unique, edge-face id should be tested to identify
        % right edge id
        if size(eid,1)==1
            fed(jj)=edSgn*eid;
        else
            efid=edge{1}(eid(1),rg.ei(2):rg.ef(2));
            efid=efid(efid==fId);
            if isempty(efid)==1
                fed(jj)=edSgn*eid(2);
            else
                fed(jj)=edSgn*eid(1);
            end               
        end
    end
    
% cells with one neighbors only
elseif size(fvr,2)==2
    
    % find edge id. There should be two edges only
    if fvr(1)<fvr(2)
        eId=find(ismember(eVr,fvr,'rows')==1 & ...
            (edge{1}(:,3)==fId | edge{1}(:,4)==fId));
        fedInd=[eId(1),-eId(2)];
    else
        eId=find(ismember(eVr,[fvr(2),fvr(1)],'rows')==1 & ...
            (edge{1}(:,3)==fId | edge{1}(:,4)==fId));
        fedInd=[-eId(1),eId(2)];
    end
    farInd=smt_faceArea(fedInd,edge{2});
    if farInd>0
        fed=fedInd;
    else
        fed=[-fedInd(2),-fedInd(1)];
    end   

% cells with no neighbors.
else
    eId=find(ismember(eVr,[fvr(1),fvr(1)],'rows')==1);
    fedInd=eId;
    farInd=smt_faceArea(fedInd,edge{2});
    if farInd>0
        fed=fedInd;
    else
        fed=-fedInd;
    end
end

end