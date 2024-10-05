%% function description.
% for a given seed point. Generate Voronoi diagram with periodic boundary
% condition first and then take cells only in [0,len*nSmt]X[0,wth]

function [vrtx,edge,face,rg,gm_p,mc_p]=...
        smt_initialConfiguration(len,wth,nSmt,edpc,ften)

%% Define output variables
% vrtx=[ed_id(3),fa_id(3),coordinate(2),force(2)
% edge{1}=[vr_id(2),fa_id(2),tension(1),refinement(1),
%   range(4):x_min,x_max,y_min,y_max]
% edge{2}=intermediate vertices
% edge{3}=intermediate vertex forces
% edge{4}=edge length
% face{1}=[type(1),area(1),perimeter(1),# of nei, # of edges, center]
% face{2}=vr_id
% face{3}=ed_id

% range variable
rg=struct;
[rg.vi,rg.vf]=deal([1,4,7,9],[3,6,8,10]);
[rg.ei,rg.ef]=deal([1,3,5,6,7],[2,4,5,6,10]);
[rg.fi,rg.ff]=deal([1,2,3,4,5,6],[1,2,3,4,5,7]);

% geometric quantity structure variable
gm_p=struct;
[gm_p.wth,gm_p.len,gm_p.nSmt,gm_p.edpc]=...
    deal(wth,len,nSmt,edpc);
gm_p.bs=gm_p.len*gm_p.nSmt;
gm_p.nFa=gm_p.bs^2;

% mechanical parameter structure variable
mc_p=struct;
mc_p.ften=ften;

%% Generate an initial configuration.
% Use lloyd algorithm for relatively regular configuration
rng('shuffle');
sdPt=smt_seedLloyd(gm_p.bs,20);        

%% Voronoi diagram generation    
% copy seed points for periodic boundary condition
spCp=smt_seedCopy(sdPt,gm_p.bs);

% generate voronoi tessellation
[vr,fVrIni]=voronoin(spCp);
fVrIni=fVrIni.';

% find cells within [0,gm_p.bs]X[0,gm_p.wth]
int_cell=(spCp(:,1)>=0 & spCp(:,1)<=gm_p.bs & ...
    spCp(:,2)>=0 & spCp(:,2)<=gm_p.wth);
int_fvr=fVrIni(int_cell);

% find vertex id that belongs to internal cells
vrId=unique(cell2mat(int_fvr));

% count number of vertex occurance for internal cells
vrCnt=zeros(size(vrId,2),1);
for fac=1:size(int_fvr,2)
    for nec=1:size(int_fvr{fac},2)
        vrCnt(vrId==int_fvr{fac}(nec))=...
            vrCnt(vrId==int_fvr{fac}(nec))+1;
    end
end

% find redundant vertex at the boundary and remove it
redunVr=vrId(vrCnt==1);
while isempty(redunVr)==0
    nei_cnt=zeros(size(int_fvr,2),1);
    for fac=1:size(int_fvr,2)
        int_fvr{fac}=setdiff(int_fvr{fac},redunVr);
        nei_cnt(fac)=size(int_fvr{fac},2);
    end
    int_fvr=int_fvr(nei_cnt>2);

    vrId=unique(cell2mat(int_fvr));
    vrCnt=zeros(size(vrId,2),1);

    for fac=1:size(int_fvr,2)
        for nec=1:size(int_fvr{fac},2)
            vrCnt(vrId==int_fvr{fac}(nec))=...
                vrCnt(vrId==int_fvr{fac}(nec))+1;
        end
    end
    redunVr=vrId(vrCnt==1);
end

%% define face and vertex information
[gm_p.nFa,gm_p.nVr]=deal(size(int_fvr,2),size(vrId,2));

vrtx=zeros(gm_p.nVr,rg.vf(end));

face=cell(3,1);
face{1}=zeros(gm_p.nFa,rg.ff(end));
[face{2},face{3}]=deal(cell(gm_p.nFa,1));

% assign vrtx coordinates
vrtx(:,7:8)=vr(vrId,:);

% assign face-vrtx id 
for fac=1:gm_p.nFa
    fvrInd=zeros(1,size(int_fvr{fac},2));
    for nec=1:size(int_fvr{fac},2)
        fvrInd(nec)=find(vrId==int_fvr{fac}(nec));
    end
    face{2}{fac}=smt_vrtxIdSort(fvrInd,vrtx(:,7:8));
end

% assign vrtx-face id 
vrFa=zeros(gm_p.nVr,3);
vrCnt=ones(gm_p.nVr,1);
for fac=1:gm_p.nFa
    fvrInd=face{2}{fac};
    for nec=1:size(fvrInd,2)
        vrFa(fvrInd(nec),vrCnt(fvrInd(nec)))=fac;
        vrCnt(fvrInd(nec))=vrCnt(fvrInd(nec))+1;
    end
end

for vrc=1:gm_p.nVr
    vrtx(vrc,rg.vi(2):rg.vf(2))=sort(vrFa(vrc,:));
end

% define edge id based on vertex information
edVr=zeros(8*gm_p.nFa,2);
edc=1;

for fac=1:gm_p.nFa
    fvrInd=[face{2}{fac},face{2}{fac}(1)];
    for nec=1:size(fvrInd,2)-1
        edVr(edc,:)=sort(fvrInd(nec:nec+1));
        edc=edc+1;
    end
end

edVr=edVr(1:edc-1,:);
edVr=unique(edVr,'rows');

% define edge data cell and assign edge-vrtx id
gm_p.nEd=size(edVr,1);
edge=cell(4,1);
edge{1}=zeros(gm_p.nEd,rg.ef(end));
[edge{2},edge{3},edge{4}]=deal(cell(gm_p.nEd,1));

edge{1}(:,rg.ei(1):rg.ef(1))=edVr;

% assign vertex-edge id
vrCnt=ones(gm_p.nVr,1);
for edc=1:gm_p.nEd
    for vrc=1:2
        evId=edge{1}(edc,vrc);
        vrtx(evId,vrCnt(evId))=edc;
        vrCnt(evId)=vrCnt(evId)+1;
    end
end

%% find edge-face relation.
for edc=1:gm_p.nEd
    edge{1}(edc,rg.ei(2):rg.ef(2))=sort(intersect(vrtx(edge{1}(edc,1),...
        rg.vi(2):rg.vf(2)),vrtx(edge{1}(edc,2),rg.vi(2):rg.vf(2))));
end  

%% define intermediate point of initial configuration
for edc=1:gm_p.nEd
    evr=vrtx(edge{1}(edc,rg.ei(1):rg.ef(1)),rg.vi(3):rg.vf(3));
    [edge{2}{edc},edge{1}(edc,rg.ei(4))]=smt_edgeMidPtAvg(evr,gm_p.edpc);               
end

%% assign face-edge relation.
for fac=1:gm_p.nFa
    face{3}{fac}=smt_faceEdgeId(edge,face,gm_p,rg,fac);
end       

%% calculate edge length.
for edc=1:gm_p.nEd 
    edge{4}{edc}=smt_edgeLen(edge{2}{edc});        
end

%% calculate area, perimeter, and number of neighbors and center.
for fac=1:gm_p.nFa
    fEd=face{3}{fac};

    [face{1}(fac,rg.fi(2)),face{1}(fac,rg.fi(6):rg.ff(6))]...
        =smt_faceArea(fEd,edge{2});
    for edc=1:size(fEd,2)
        face{1}(fac,rg.fi(3))=face{1}(fac,rg.fi(3))+...
            sum(edge{4}{abs(fEd(edc))});
    end
    face{1}(fac,rg.fi(4))=smt_faceNei(fEd,edge{1}(:,rg.ei(2):rg.ef(2)),fac);
    face{1}(fac,rg.fi(5))=size(fEd,2);
end

%% assign face type.
% face type is based on x coordinates with segment length of gm_p.len
for sgc=1:gm_p.nSmt-1
    seg_id=((face{1}(:,rg.fi(6))<gm_p.len*sgc)&(face{1}(:,rg.fi(1))==0));
    face{1}(seg_id,rg.fi(1))=sgc;
end
seg_id=(face{1}(:,rg.fi(1))==0);
face{1}(seg_id,rg.fi(1))=gm_p.nSmt;

%% calculate edge tension
edge{1}((edge{1}(:,3)>0),rg.ei(3))=1;
edge{1}((edge{1}(:,3)==0),rg.ei(3))=mc_p.ften;

%% Calculate edge range
for edc=1:gm_p.nEd
    edge{1}(edc,rg.ei(5):rg.ef(5))=smt_edgeRng(edge,edc);
end

end