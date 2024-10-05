%% Adjust vertex position in force direction.
% For a given set of vertex coordinates, calculate force for each vertex
% and displace vertex by force*scFac*edLnAv. 

function [vrtx_n,edge_n,face_n]=smt_iteration(vrtx,edge,face,rg,gm_p,mc_p)
    
%% Declare new variables
[vrtx_n,edge_n,face_n]=deal(vrtx,edge,face);
vrtx_n(:,rg.vi(4):rg.vf(4))=zeros(gm_p.nVr,2);

%% Compute net force for intermediate vertices
for edc=1:gm_p.nEd
    rfn=edge{1}(edc,rg.ei(4));
    vTn=zeros(rfn+2,2);
    eMd=edge{2}{edc};
    fCn=face{1}(edge{1}(edc,rg.ef(2)),rg.fi(6):rg.ff(6));
    eTf=max(0,edge{1}(edc,rg.ei(3)));

    % for each segment, compute tangential and normal force
    for jj=1:rfn+1
        [tv,nv,eln]=smt_edgeVector(eMd(jj,:),eMd(jj+1,:),fCn);
        eNf=smt_edgeNormalFrc(edge,face,rg,edc,mc_p.del,eln,gm_p);
        vTn(jj,:)=vTn(jj,:)+eTf*tv;
        vTn(jj+1,:)=vTn(jj+1,:)-eTf*tv;
        vTn(jj,:)=vTn(jj,:)+eNf*nv;
        vTn(jj+1,:)=vTn(jj+1,:)+eNf*nv;
    end
    
    % update end vertex forces
    vrtx_n(edge_n{1}(edc,rg.ei(1)),rg.vi(4):rg.vf(4))=...
        vrtx_n(edge_n{1}(edc,rg.ei(1)),rg.vi(4):rg.vf(4))+vTn(1,:);
    vrtx_n(edge_n{1}(edc,rg.ef(1)),rg.vi(4):rg.vf(4))=...
        vrtx_n(edge_n{1}(edc,rg.ef(1)),rg.vi(4):rg.vf(4))+vTn(end,:);

    edge_n{3}{edc}=vTn;        
    edge_n{2}{edc}=eMd+gm_p.dt*vTn;
end

%% Adjust end vertices and update intermediate vertex information
vrtx_n(:,rg.vi(3):rg.vf(3))=vrtx_n(:,rg.vi(3):rg.vf(3))+...
    gm_p.dt*vrtx_n(:,rg.vi(4):rg.vf(4));

% update end vertex coordinate in intermediate vertex coordinates
for edc=1:gm_p.nEd
    eMd=edge_n{2}{edc};
    eMd(1,:)=vrtx_n(edge_n{1}(edc,rg.ei(1)),rg.vi(3):rg.vf(3));
    eMd(end,:)=vrtx_n(edge_n{1}(edc,rg.ef(1)),rg.vi(3):rg.vf(3));
    edge_n{2}{edc}=eMd;        
end    

%% Update edge tension.
[eTn,eTp]=deal(zeros(size(edge{1},1),1));
for edc=1:gm_p.nEd
    [eTn(edc),eTp(edc)]=smt_edgeTension(face_n,edge_n,rg,mc_p,gm_p,edc);
end

edge_n{1}(:,rg.ei(3))=edge_n{1}(:,rg.ei(3))...
    -gm_p.dt/mc_p.taut*(edge_n{1}(:,rg.ei(3))-eTn)...
    +mc_p.mu/mc_p.taut*sqrt(gm_p.dt)/2*eTp...
    .*normrnd(0,1,gm_p.nEd,1);
%     
%% Check refinement and adjust number of intermediate vertices
for edc=1:gm_p.nEd
    emd=edge_n{2}{edc};
    eln=smt_edgeLen(emd);             
    if sum(eln)<gm_p.edpc*(edge_n{1}(edc,rg.ei(4))-1/3) || ...
            sum(eln)>gm_p.edpc*(edge_n{1}(edc,rg.ei(4))+4/3) || ...
            (max(eln)/min(eln)>2)                
        [edge_n{2}{edc},edge_n{1}(edc,rg.ei(4))]=...
            smt_edgeMidPtAvg(emd,gm_p.edpc);
        edge_n{4}{edc}=smt_edgeLen(edge_n{2}{edc});
    else
        edge_n{4}{edc}=eln;
    end
    edge_n{1}(edc,rg.ei(5):rg.ef(5))=smt_edgeRng(edge_n,edc);
end

%% calculate area, perimeter, and number of neighbors and center.
for fac=1:gm_p.nFa
    fEd=face_n{3}{fac};

    [face_n{1}(fac,rg.fi(2)),face_n{1}(fac,rg.fi(6):rg.ff(6))]...
        =smt_faceArea(fEd,edge{2});
    face_n{1}(fac,rg.fi(3))=0;
    for jj=1:size(fEd,2)
        face_n{1}(fac,rg.fi(3))=face_n{1}(fac,rg.fi(3))+...
            sum(edge{4}{abs(fEd(jj))});
    end
    face_n{1}(fac,rg.fi(4))=smt_faceNei(fEd,edge{1}(:,rg.ei(2):rg.ef(2)),fac);
    face_n{1}(fac,rg.fi(5))=size(fEd,2);
end

%     vrtx_n=crdDcFn(vrtx_n,rg,gm_p);
end