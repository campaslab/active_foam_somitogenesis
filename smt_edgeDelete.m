%% Merge two cells into one cell.

function [vrtx_n,edge_n,face_n,gm_p_n,eVl]=...
    smt_edgeDelete(vrtx,edge,face,rg,gm_p,deId)

%% Define output variables.
[vrtx_n,edge_n,face_n,gm_p_n]=deal(vrtx,edge,face,gm_p);
[vMx,eMx]=deal(size(vrtx_n,1),size(edge_n{1},1));

%% Find local geometric elements.
% Center vertex ID
cvId=edge_n{1}(deId,rg.ei(1):rg.ef(1));

% Local face ID
if edge_n{1}(deId,rg.ei(2))~=0 || edge_n{1}(deId,rg.ef(2))~=0
    error('This is not an edge between extracellular space.');
end

fId=[setdiff(vrtx_n(cvId(1),rg.vi(2):rg.vf(2)),0),...
    setdiff(vrtx_n(cvId(2),rg.vi(2):rg.vf(2)),0)];

% Local edge ID
eId=zeros(1,4);
for ii=[1,3]
    veid=vrtx_n(cvId(floor(ii/2)+1),rg.vi(1):rg.vf(1));
    veid=veid(veid~=0);
    eId(ii:ii+1)=veid(veid~=deId);
end

% Local vertex ID
vId=zeros(1,4);
for ii=[1,3]
    if eId(ii)~=eId(ii+1)
        vId(ii:ii+1)=...
            [setdiff(edge_n{1}(eId(ii),rg.ei(1):rg.ef(1)),cvId(floor(ii/2)+1)),...
            setdiff(edge_n{1}(eId(ii+1),rg.ei(1):rg.ef(1)),cvId(floor(ii/2)+1))];
    else
        vId(ii:ii+1)=cvId(floor(ii/2)+1);
    end
end

%% Reassign vertex-edge relation
if eId(1)~=eId(2)
    vEd=vrtx_n(vId(2),rg.vi(1):rg.vf(1));
    vEd=vEd(vEd~=eId(2));
    vrtx_n(vId(2),rg.vi(1):rg.vf(1))=sort([eId(1),vEd]);
else
    vrtx_n(vId(2),rg.vi(1):rg.vf(1))=[0,0,eId(1)];
end

if eId(3)~=eId(4)
    vEd=vrtx_n(vId(4),rg.vi(1):rg.vf(1));
    vEd=vEd(vEd~=eId(4));
    vrtx_n(vId(4),rg.vi(1):rg.vf(1))=sort([eId(3),vEd]);
else
    vrtx_n(vId(4),rg.vi(1):rg.vf(1))=[0,0,eId(3)];
end

%% Reassign edge-vertex relation
edge_n{1}(eId(1),rg.ei(1):rg.ef(1))=sort([vId(1),vId(2)]);
edge_n{1}(eId(3),rg.ei(1):rg.ef(1))=sort([vId(3),vId(4)]);

%% Recompute edge tension value.
edge_n{1}(eId(1),rg.ei(3))=...
    (edge_n{1}(eId(1),rg.ei(3))+edge_n{1}(eId(2),rg.ei(3)))/2;
edge_n{1}(eId(3),rg.ei(3))=...
    (edge_n{1}(eId(3),rg.ei(3))+edge_n{1}(eId(4),rg.ei(3)))/2;

%% Update intermediate edge coordinates.
if eId(1)~=eId(2)
    emd1=edge{2}{eId(1)};
    if edge{1}(eId(1),1)==cvId(1)
        emd1=flipud(emd1);
    end
    emd2=edge{2}{eId(2)};
    if edge{1}(eId(2),2)==cvId(1)
        emd2=flipud(emd2);
    end
    vr=[emd1;emd2(2:end,:)];
    if edge_n{1}(eId(1),1)==vId(2)
        vr=flipud(vr);
    end
else
    vr=edge{2}{eId(1)};
end

[edge_n{2}{eId(1)},edge_n{1}(eId(1),rg.ei(4))]=...
    smt_edgeMidPtAvg(vr,gm_p.edpc);
edge_n{1}(eId(1),rg.ei(5):rg.ef(5))=smt_edgeRng(edge_n,eId(1));

if eId(3)~=eId(4)
    emd1=edge{2}{eId(3)};
    if edge{1}(eId(3),1)==cvId(2)
        emd1=flipud(emd1);
    end
    emd2=edge{2}{eId(4)};
    if edge{1}(eId(4),2)==cvId(2)
        emd2=flipud(emd2);
    end
    vr=[emd1;emd2(2:end,:)];
    if edge_n{1}(eId(3),1)==vId(4)
        vr=flipud(vr);
    end
else
    vr=edge{2}{eId(3)};
end

[edge_n{2}{eId(3)},edge_n{1}(eId(3),rg.ei(4))]=...
    smt_edgeMidPtAvg(vr,gm_p.edpc);
edge_n{1}(eId(3),rg.ei(5):rg.ef(5))=smt_edgeRng(edge_n,eId(3));

%% Update edge length.     
edge_n{4}{eId(1)}=smt_edgeLen(edge_n{2}{eId(1)}); 
edge_n{4}{eId(3)}=smt_edgeLen(edge_n{2}{eId(3)}); 

%% Reassign face-vertex.
if vId(1)~=vId(2)
    fvr=face_n{2}{fId(1)};
    fvr=fvr(fvr~=cvId(1));
    face_n{2}{fId(1)}=fvr;
end

if vId(3)~=vId(4)
    fvr=face_n{2}{fId(2)};
    fvr=fvr(fvr~=cvId(2));
    face_n{2}{fId(2)}=fvr;
end

%% Reassign face-edge relation    
for ii=1:2
    face_n{3}{fId(ii)}=smt_faceEdgeId(edge_n,face_n,gm_p,rg,fId(ii));
end

%% Recalculate area value. 
for ii=1:2
    fEd=face_n{3}{fId(ii)};

    [face_n{1}(fId(ii),rg.fi(2)),face_n{1}(fId(ii),rg.fi(6):rg.ff(6))]...
        =smt_faceArea(fEd,edge_n{2});
    face_n{1}(fId(ii),rg.fi(3))=0;
    for jj=1:size(fEd,2)
        face_n{1}(fId(ii),rg.fi(3))=face_n{1}(fId(ii),rg.fi(3))+...
            sum(edge_n{4}{abs(fEd(jj))});
    end
    face_n{1}(fId(ii),rg.fi(4))=smt_faceNei(fEd,edge_n{1}...
        (:,rg.ei(2):rg.ef(2)),fId(ii));
    face_n{1}(fId(ii),rg.fi(5))=size(fEd,2);
end

%% Reassign face-vertex relation and update ID number.
vlReId=cvId;
elReId=[deId,eId(2),eId(4)];

if vId(1)==vId(2)
    vlReId=vlReId(vlReId~=vId(1));
    elReId=elReId(elReId~=eId(2));
end

if vId(3)==vId(4)
    vlReId=vlReId(vlReId~=vId(3));
    elReId=elReId(elReId~=eId(4));
end

vVl=setdiff(1:vMx,vlReId);
eVl=setdiff(1:eMx,elReId);

%% Taking only valid geometric elements.
vrtx_n=vrtx_n(vVl,:);
edge_n{1}=edge_n{1}(eVl,:);
edge_n{2}=edge_n{2}(eVl);
edge_n{3}=edge_n{3}(eVl);
edge_n{4}=edge_n{4}(eVl);

%% Reassigning new IDs for geometric elements.   
[vrtx_n,edge_n,face_n]=smt_newIdAssgn(vrtx_n,edge_n,face_n,rg,vVl,eVl);  

gm_p_n.nVr=size(vrtx_n,1);
gm_p_n.nEd=size(edge_n{1},1);

end