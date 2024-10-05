%% Reassign geometic ID

function [vrtx_n,edge_n,face_n]=smt_newIdAssgn(vrtx,edge,face,rg,vVl,eVl)
    
[vrtx_n,edge_n,face_n]=deal(vrtx,edge,face);
        
% face id
for ii=1:size(face{2},1)
    fvr=face{2}{ii};
    fed=face{3}{ii};
    for jj=1:size(fvr,2)
        fvr(jj)=find(vVl==fvr(jj));
        fed(jj)=sign(fed(jj))*find(eVl==abs(fed(jj)));
    end
    face_n{2}{ii}=fvr;
    face_n{3}{ii}=fed;
end

% vertex id
for ii=1:size(vrtx,1)
    ved=vrtx_n(ii,rg.vi(1):rg.vf(1));
    for jj=1:3
        if ved(jj)>0
            vrtx_n(ii,rg.vi(1)-1+jj)=find(eVl==ved(jj));
        end
    end
end

% edge id
for ii=1:size(edge{1},1)
    for jj=1:2
        edge_n{1}(ii,rg.ei(1)-1+jj)=find(vVl==edge_n{1}(ii,rg.ei(1)-1+jj));
    end
end 

end