%% Merge two cells into one cell.
function [vrtxN,edgeN,faceN,gmpN,eVl]=...
    smt_t3ReverseTransition(vrtx,edge,face,rg,gmp,deId)

%% Define output variables.
[vrtxN,edgeN,faceN,gmpN]=deal(vrtx,edge,face,gmp);
[vMx,eMx]=deal(size(vrtxN,1),size(edgeN{1},1));

%% Find local geometric elements.
% Center vertex ID
cvId=edgeN{1}(deId,rg.ei(1):rg.ef(1));

% Local face ID
% if edgeN{1}(deId,rg.ei(2))~=0 || edgeN{1}(deId,rg.ef(2))~=0
%     error('This is not an edge between extracellular space.');
% end

fId=[edgeN{1}(deId,rg.ei(2):rg.ef(2)).';
    setdiff(vrtxN(cvId(1),rg.vi(2):rg.vf(2)),edgeN{1}(deId,rg.ei(2):rg.ef(2)));...
    setdiff(vrtxN(cvId(2),rg.vi(2):rg.vf(2)),edgeN{1}(deId,rg.ei(2):rg.ef(2)))];

% Local edge ID
eId=zeros(1,4);
eId(1)=find(ismember(edge{1}(:,rg.ei(2):rg.ef(2)),...
    sort([fId(1),fId(3)]),'rows'));
eId(2)=find(ismember(edge{1}(:,rg.ei(2):rg.ef(2)),...
    sort([fId(2),fId(3)]),'rows'));
eId(3)=find(ismember(edge{1}(:,rg.ei(2):rg.ef(2)),...
    sort([fId(1),fId(4)]),'rows'));
eId(4)=find(ismember(edge{1}(:,rg.ei(2):rg.ef(2)),...
    sort([fId(2),fId(4)]),'rows'));

% for ii=[1,3]
%     veid=vrtxN(cvId(floor(ii/2)+1),rg.vi(1):rg.vf(1));
% %     veid=veid(veid~=0);
%     eId(ii:ii+1)=veid(veid~=deId);
% end

% Local vertex ID
vId=zeros(1,4);
for edc=1:4
    vId(edc)=setdiff(edgeN{1}(eId(edc),rg.ei(1):rg.ef(1)),cvId);
end
% for ii=[1,3]
%     if eId(ii)~=eId(ii+1)
%         vId(ii:ii+1)=...
%             [setdiff(edgeN{1}(eId(ii),rg.ei(1):rg.ef(1)),cvId(floor(ii/2)+1)),...
%             setdiff(edgeN{1}(eId(ii+1),rg.ei(1):rg.ef(1)),cvId(floor(ii/2)+1))];
%     else
%         vId(ii:ii+1)=cvId(floor(ii/2)+1);
%     end
% end

%% Reassign vertex-edge relation
% if eId(1)~=eId(2)
vEd=vrtxN(vId(2),rg.vi(1):rg.vf(1));
vEd=vEd(vEd~=eId(2));
vrtxN(vId(2),rg.vi(1):rg.vf(1))=sort([eId(1),vEd]);
% else
%     vrtxN(vId(2),rg.vi(1):rg.vf(1))=[0,0,eId(1)];
% end

% if eId(3)~=eId(4)
vEd=vrtxN(vId(4),rg.vi(1):rg.vf(1));
vEd=vEd(vEd~=eId(4));
vrtxN(vId(4),rg.vi(1):rg.vf(1))=sort([eId(3),vEd]);
% else
%     vrtxN(vId(4),rg.vi(1):rg.vf(1))=[0,0,eId(3)];
% end

%% Reassign vertex-face relation
fvr=faceN{2}{fId(2)};

for vrc=1:numel(fvr)
    vFa=vrtxN(fvr(vrc),rg.vi(2):rg.vf(2));
    vFa=vFa(vFa~=fId(2));
    vrtxN(fvr(vrc),rg.vi(2):rg.vf(2))=sort([fId(1),vFa]);
end


%% Reassign edge-vertex relation
edgeN{1}(eId(1),rg.ei(1):rg.ef(1))=sort([vId(1),vId(2)]);
edgeN{1}(eId(3),rg.ei(1):rg.ef(1))=sort([vId(3),vId(4)]);

%% Reassign edge-face relation
fed=abs(face{3}{fId(2)});

for edc=1:numel(fed)
    eFa=edge{1}(fed(edc),rg.ei(2):rg.ef(2));
    eFa(eFa==fId(2))=fId(1);
    edgeN{1}(fed(edc),rg.ei(2):rg.ef(2))=sort(eFa);
end

%% Recompute edge tension value.
edgeN{1}(eId(1),rg.ei(3))=...
    (edgeN{1}(eId(1),rg.ei(3))+edgeN{1}(eId(2),rg.ei(3)))/2;
edgeN{1}(eId(3),rg.ei(3))=...
    (edgeN{1}(eId(3),rg.ei(3))+edgeN{1}(eId(4),rg.ei(3)))/2;

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
    if edgeN{1}(eId(1),1)==vId(2)
        vr=flipud(vr);
    end
else
    vr=edge{2}{eId(1)};
end

% vr=ts_crdLocal(vr,gmp.bs,gmp.sstn);
[edgeN{2}{eId(1)},edgeN{1}(eId(1),rg.ei(4))]=...
    smt_edgeMidPtAvg(vr,gmp.edpc);
edgeN{1}(eId(1),rg.ei(5):rg.ef(5))=smt_edgeRng(edgeN,eId(1));

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
    if edgeN{1}(eId(3),1)==vId(4)
        vr=flipud(vr);
    end
else
    vr=edge{2}{eId(3)};
end

% vr=ts_crdLocal(vr,gmp.bs,gmp.sstn);
[edgeN{2}{eId(3)},edgeN{1}(eId(3),rg.ei(4))]=...
    smt_edgeMidPtAvg(vr,gmp.edpc);
edgeN{1}(eId(3),rg.ei(5):rg.ef(5))=smt_edgeRng(edgeN,eId(3));

%% Update edge length.     
edgeN{4}{eId(1)}=smt_edgeLen(edgeN{2}{eId(1)}); 
edgeN{4}{eId(3)}=smt_edgeLen(edgeN{2}{eId(3)}); 


%% Reassign face-vertex.

fvr1=face{2}{fId(1)};
vrPs=find(fvr1==cvId(1));
fvr1=[fvr1(vrPs:numel(fvr1)),fvr1(1:vrPs-1)];
if fvr1(2)~=cvId(2)
    fvr1=[fvr1(end),fvr1(1:end-1)];
end
fvr2=face{2}{fId(2)};
vrPs=find(fvr2==cvId(1));
fvr2=[fvr2(vrPs:numel(fvr2)),fvr2(1:vrPs-1)];
if fvr2(2)~=cvId(2)
    fvr2=[fvr2(end),fvr2(1:end-1)];
end
faceN{2}{fId(1)}=[fvr1(3:end),fvr2(3:end)];


if vId(1)~=vId(2)
    fvr=faceN{2}{fId(3)};
    fvr=fvr(fvr~=cvId(1));
    faceN{2}{fId(3)}=fvr;
end

if vId(3)~=vId(4)
    fvr=faceN{2}{fId(4)};
    fvr=fvr(fvr~=cvId(2));
    faceN{2}{fId(4)}=fvr;
end

%% Reassign face-edge relation    
for ii=[1,3,4]
    faceN{3}{fId(ii)}=smt_faceEdgeId(edgeN,faceN,gmp,rg,fId(ii));
end

%% Recalculate area value. 
for ii=[1,3,4]
    fEd=faceN{3}{fId(ii)};

    [faceN{1}(fId(ii),rg.fi(2)),faceN{1}(fId(ii),rg.fi(6):rg.ff(6))]...
        =smt_faceArea(fEd,edgeN{2});
    faceN{1}(fId(ii),rg.fi(3))=0;
    for jj=1:size(fEd,2)
        faceN{1}(fId(ii),rg.fi(3))=faceN{1}(fId(ii),rg.fi(3))+...
            sum(edgeN{4}{abs(fEd(jj))});
    end
    faceN{1}(fId(ii),rg.fi(4))=...
        smt_faceNei(fEd,edgeN{1}(:,rg.ei(2):rg.ef(2)),fId(ii));
    faceN{1}(fId(ii),rg.fi(5))=size(fEd,2);
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
fVl=setdiff(1:gmpN.nFa,fId(2));

%% Reassign face type as nSmt+1
faceN{1}(fId(1),rg.fi(1))=gmpN.nSmt+1;

%% Taking only valid geometric elements.
vrtxN=vrtxN(vVl,:);
edgeN{1}=edgeN{1}(eVl,:);
edgeN{2}=edgeN{2}(eVl);
edgeN{3}=edgeN{3}(eVl);
edgeN{4}=edgeN{4}(eVl);
faceN{1}=faceN{1}(fVl,:);
faceN{2}=faceN{2}(fVl);
faceN{3}=faceN{3}(fVl);

%% Take one cell count out
gmpN.nFa=gmpN.nFa-1;

%% Reassigning new IDs for geometric elements.   
[vrtxN,edgeN,faceN]=smt_newIdAssign(vrtxN,edgeN,faceN,rg,vVl,eVl,fVl);

end