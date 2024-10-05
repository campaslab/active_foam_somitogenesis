%% Execute a single T1 trnasition. 

function [vrtx_n,edge_n,face_n]=smt_t1Flip(vrtx,edge,face,rg,gm_p,mc_p,t1Id)

%% Define output variables.
[vrtx_n,edge_n,face_n]=deal(vrtx,edge,face);  

%% Find vertex ID, edge ID, and face ID of local T1 configuration.
% local edge ID. 
cvId=edge_n{1}(t1Id,rg.ei(1):rg.ef(1));     
eId=[setdiff(vrtx_n(cvId(1),rg.vi(1):rg.vf(1)),t1Id),...
    setdiff(vrtx_n(cvId(2),rg.vi(1):rg.vf(1)),t1Id)];
ef1=edge_n{1}(eId(1),rg.ei(2):rg.ef(2));
ef1=ef1(ef1~=0);
if size(ef1,2)==1
    eId(1:2)=[eId(2),eId(1)];
    ef1=edge_n{1}(eId(1),rg.ei(2):rg.ef(2));
    ef1=ef1(ef1~=0);
end

ef2=edge_n{1}(eId(4),rg.ei(2):rg.ef(2));
ef2=ef2(ef2~=0);
if size(ef2,2)==1
    eId(3:4)=[eId(4),eId(3)];
    ef2=edge_n{1}(eId(4),rg.ei(2):rg.ef(2));
    ef2=ef2(ef2~=0);
end

if isempty(intersect(ef1,ef2))==0
    eId(3:4)=[eId(4),eId(3)];
end

% local vertex ID. 
vId=zeros(4,1);    
for ii=1:2
    vId(ii)=setdiff(edge_n{1}(eId(ii),rg.ei(1):rg.ef(1)),cvId(1));
end
for ii=3:4
    vId(ii)=setdiff(edge_n{1}(eId(ii),rg.ei(1):rg.ef(1)),cvId(2));
end

% local face ID.
fId=[intersect(edge_n{1}(t1Id,3:4),edge_n{1}(eId(1),3:4)),...
    intersect(edge_n{1}(t1Id,3:4),edge_n{1}(eId(2),3:4)),...
    setdiff(edge_n{1}(eId(1),3:4),edge_n{1}(t1Id,3:4)),...
    setdiff(edge_n{1}(eId(3),3:4),edge_n{1}(t1Id,3:4))];

%% Reassign vertex-edge relation.
if eId(1)~=eId(3)
    vrtx_n(cvId(1),rg.vi(1):rg.vf(1))=sort([t1Id,eId(1),eId(3)]);
else
    vrtx_n(cvId(1),rg.vi(1):rg.vf(1))=sort([0,t1Id,eId(1)]);
end

if eId(2)~=eId(4)
    vrtx_n(cvId(2),rg.vi(1):rg.vf(1))=sort([t1Id,eId(2),eId(4)]);
else
    vrtx_n(cvId(2),rg.vi(1):rg.vf(1))=sort([0,t1Id,eId(2)]);
end


%% Reassign vertex-face relation.
vrtx_n(cvId(1),rg.vi(2):rg.vf(2))=sort([fId(1),fId(3),fId(4)]);
vrtx_n(cvId(2),rg.vi(2):rg.vf(2))=sort([fId(2),fId(3),fId(4)]);                 

%% Reassign edge-vertex relation.
edge_n{1}(eId(2),rg.ei(1):rg.ef(1))=sort([cvId(2),vId(2)]);
edge_n{1}(eId(3),rg.ei(1):rg.ef(1))=sort([cvId(1),vId(3)]);

%% Reassign edge-face relation.
edge_n{1}(t1Id,rg.ei(2):rg.ef(2))=sort([fId(3),fId(4)]);

%% Adjust vertex coordinates for T1 vertices. 
vCrdTmp=vrtx_n(cvId,rg.vi(3):rg.vf(3));
t1CrdTmp=repmat(1/2*sum(vCrdTmp(:)),2)-[vCrdTmp(2,2),vCrdTmp(1,1);...
    vCrdTmp(1,2),vCrdTmp(2,1)];

emd1=edge_n{2}{eId(1)};
if edge{1}(eId(1),1)==vId(1)
    emd1(end,:)=t1CrdTmp(1,:);
    if vId(1)==cvId(2)
        emd1=emd1(2:end,:);
    end            
else
    emd1(1,:)=t1CrdTmp(1,:);
    if vId(1)==cvId(2)
        emd1=emd1(1:end-1,:);
    end
end

emd2=edge_n{2}{eId(2)};
if edge{1}(eId(2),1)==vId(2)
    emd2(end,:)=t1CrdTmp(2,:);
    if vId(2)==cvId(2)
        emd2=emd2(2:end,:);
    end 
else
    emd2(1,:)=t1CrdTmp(2,:);
    if vId(2)==cvId(2)
        emd1=emd1(1:end-1,:);
    end 
end

emdc=[emd1;emd2];
emd1=emdc(1:size(emd1,1),:);
emd2=emdc(size(emd1,1)+1:end,:);

ec=0;
for ii=1:size(emd1,1)-1
    for jj=1:size(emd2,1)-1
        [ec_cs,~]=smt_lineCross([emd1(ii:ii+1,:);emd2(jj:jj+1,:)]);
        if ec_cs==1
            ec=1;
        end
    end
end

if ec==1
    vrtx_n(cvId(1),rg.vi(3):rg.vf(3))=t1CrdTmp(2,:);
    vrtx_n(cvId(2),rg.vi(3):rg.vf(3))=t1CrdTmp(1,:);
else
    vrtx_n(cvId(1),rg.vi(3):rg.vf(3))=t1CrdTmp(1,:);
    vrtx_n(cvId(2),rg.vi(3):rg.vf(3))=t1CrdTmp(2,:);
end

%% Update intermediate edge coordinates.
eVr=vrtx_n(edge_n{1}(t1Id,rg.ei(1):rg.ef(1)),rg.vi(3):rg.vf(3));
[edge_n{2}{t1Id},edge_n{1}(t1Id,rg.ei(4))]=smt_edgeMidPtAvg(eVr,gm_p.edpc);
if eId(1)~=eId(3)
    for ii=[1,3]
        if edge{1}(eId(ii),1)==vId(ii)
            emd=edge{2}{eId(ii)}(1:end-1,:);
            if edge_n{1}(eId(ii),1)==vId(ii)
                emdc=[emd;vrtx_n(cvId(1),rg.vi(3):rg.vf(3))];
            else
                emdc=[vrtx_n(cvId(1),rg.vi(3):rg.vf(3));flipud(emd)];
            end
        else
            emd=edge{2}{eId(ii)}(2:end,:);
            if edge_n{1}(eId(ii),2)==vId(ii)
                emdc=[vrtx_n(cvId(1),rg.vi(3):rg.vf(3));emd];
            else
                emdc=[flipud(emd);vrtx_n(cvId(1),rg.vi(3):rg.vf(3))];
            end
        end            
        eMd=emdc;
        [edge_n{2}{eId(ii)},edge_n{1}(eId(ii),rg.ei(4))]...
            =smt_edgeMidPtAvg(eMd,gm_p.edpc);
    end
else
    emd=edge{2}{eId(1)}(2:end-1,:);
    emdc=[vrtx_n(cvId(1),rg.vi(3):rg.vf(3));emd;...
        vrtx_n(cvId(1),rg.vi(3):rg.vf(3))];
    eMd=emdc;
    [edge_n{2}{eId(1)},edge_n{1}(eId(1),rg.ei(4))]...
        =smt_edgeMidPtAvg(eMd,gm_p.edpc);
end

if eId(2)~=eId(4)
    for ii=[2,4]
        if edge{1}(eId(ii),1)==vId(ii)
            emd=edge{2}{eId(ii)}(1:end-1,:);
            if edge_n{1}(eId(ii),1)==vId(ii)
                emdc=[emd;vrtx_n(cvId(2),rg.vi(3):rg.vf(3))];
            else
                emdc=[vrtx_n(cvId(2),rg.vi(3):rg.vf(3));flipud(emd)];
            end
        else
            emd=edge{2}{eId(ii)}(2:end,:);
            if edge_n{1}(eId(ii),2)==vId(ii)
                emdc=[vrtx_n(cvId(2),rg.vi(3):rg.vf(3));emd];
            else
                emdc=[flipud(emd);vrtx_n(cvId(2),rg.vi(3):rg.vf(3))];
            end
        end            
        eMd=emdc;
        [edge_n{2}{eId(ii)},edge_n{1}(eId(ii),rg.ei(4))]...
            =smt_edgeMidPtAvg(eMd,gm_p.edpc);
    end
else
    emd=edge{2}{eId(2)}(2:end-1,:);
    emdc=[vrtx_n(cvId(2),rg.vi(3):rg.vf(3));emd;...
        vrtx_n(cvId(2),rg.vi(3):rg.vf(3))];
    eMd=emdc;
    [edge_n{2}{eId(2)},edge_n{1}(eId(2),rg.ei(4))]...
        =smt_edgeMidPtAvg(eMd,gm_p.edpc);
end
  
%% Calculate new edge length.
edge_n{4}{t1Id}=smt_edgeLen(edge_n{2}{t1Id});
edge_n{1}(t1Id,rg.ei(5):rg.ef(5))=smt_edgeRng(edge_n,t1Id);
for ii=1:4
    edge_n{4}{eId(ii)}=smt_edgeLen(edge_n{2}{eId(ii)});
    edge_n{1}(eId(ii),rg.ei(5):rg.ef(5))=smt_edgeRng(edge_n,eId(ii));
end

%% Reassign face-vertex relation.
% Insert vertex 2 in face 1 and vertex 1 in face 2. 
for ii=[1,2]
    if fId(ii)>0
        fvrInd=face_n{2}{fId(ii)};
        face_n{2}{fId(ii)}=fvrInd(fvrInd~=cvId(mod(ii,2)+1));
    end
end

for ii=[3,4]
    if fId(ii)>0
        fvrInd=face_n{2}{fId(ii)};
        fvrInd=smt_faceVrIdIn(cvId(mod(ii-1,2)+1),vId(ii-1),...
            cvId(mod(ii,2)+1),fvrInd);
        face_n{2}{fId(ii)}=fvrInd;
    end
end
     
%% Reassign face-edge relation.
for ii=1:4
    if fId(ii)>0
        face_n{3}{fId(ii)}=smt_faceEdgeId(edge_n,face_n,gm_p,rg,fId(ii));
    end
end        

%% Calculate new face area, number of neighbors, perimeter.
for ii=1:4
    if fId(ii)>0
        fEd=face_n{3}{fId(ii)};

        [face_n{1}(fId(ii),rg.fi(2)),face_n{1}(fId(ii),rg.fi(6):rg.ff(6))]...
            =smt_faceArea(fEd,edge_n{2});
        if face_n{1}(fId(ii),rg.fi(2))<0
            face_n{2}{fId(ii)}=fliplr(face_n{2}{fId(ii)});
            face_n{3}{fId(ii)}=smt_faceEdgeId(edge_n,face_n,gm_p,rg,fId(ii));
            fEd=face_n{3}{fId(ii)};
            [face_n{1}(fId(ii),rg.fi(2)),face_n{1}(fId(ii),rg.fi(6):rg.ff(6))]...
                =smt_faceArea(fEd,edge_n{2});
        end
        face_n{1}(fId(ii),rg.fi(3))=0;
        for jj=1:size(fEd,2)
            face_n{1}(fId(ii),rg.fi(3))=face_n{1}(fId(ii),rg.fi(3))+...
                sum(edge_n{4}{abs(fEd(jj))});
        end
        face_n{1}(fId(ii),rg.fi(4))=smt_faceNei(fEd,edge_n{1}(:,rg.ei(2):rg.ef(2)),...
            fId(ii));
        face_n{1}(fId(ii),rg.fi(5))=size(fEd,2);
    end
end 

%% Compute edge tension for T1 edge.
eTn=[smt_edgeTension(face,edge,rg,mc_p,gm_p,t1Id),...
    smt_edgeTension(face_n,edge_n,rg,mc_p,gm_p,t1Id)];

edge_n{1}(t1Id,rg.ei(3))=edge_n{1}(t1Id,rg.ei(3))-(eTn(1)-eTn(2));

end