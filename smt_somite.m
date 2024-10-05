%% Run somite simulations
% This is open boundary simulations with confluent condition. To simulation
% somite formation, we introduce heterotypic tension between somite
% boundary

function [vrtx,edge,face,gm_p,mc_p,rg]=...
    smt_somite(len,wth,nSmt,edpc,ften,beta,taut,mu,del,hten,taus,tauh,tmLag,rpc)

%% Generate an initial configuration.
[errCnt,repCnt]=deal(0);

while errCnt==repCnt
    try
        [vrtx,edge,face,rg,gm_p,mc_p]=...
            smt_initialConfiguration(len,wth,nSmt,edpc,ften);
    catch
        errCnt=errCnt+1;
    end
    repCnt=repCnt+1;
end

clear err_cnt rep_cnt

%% Assign model parameter values to structure variables
[gm_p.eLnc,gm_p.dt]=deal(smt_edgeLenAll(edge,gm_p),0.005);
[mc_p.beta,mc_p.taut,mc_p.del,mc_p.mu,mc_p.hten,mc_p.taus,mc_p.tauh,mc_p.tmLag]=...
    deal(beta,taut,del,mu,hten,taus,tauh,tmLag);
mc_p.hten_ind=zeros(gm_p.nSmt-1,1);

%% Initial relaxation of configuration. 
shC=10;
gm_p.shEd=0.001*2*sqrt(pi)*shC;

for itc=1:10000
    [vrtx,edge,face]=smt_iteration(vrtx,edge,face,rg,gm_p,mc_p);

    if mod(itc,3)==0
        eLn=smt_edgeLenAll(edge,gm_p);
        t1Id=find(eLn<gm_p.shEd);
        if isempty(t1Id)==0
            for jj=1:size(t1Id,1)
                if gm_p.eLnc(t1Id(jj))>eLn(t1Id(jj))
                    [vrtx,edge,face]=...
                        smt_t1Flip(vrtx,edge,face,rg,gm_p,mc_p,t1Id(jj));                                                            
                end
            end
            
            zrId=find(edge{1}(:,4)==0);
            while ~isempty(zrId)
                [vrtx,edge,face,gm_p,~]=...
                    smt_edgeDelete(vrtx,edge,face,rg,gm_p,zrId(1));
                zrId=find(edge{1}(:,4)==0);
            end
        else
        
            if shC<90
                shC=shC+1;
                gm_p.shEd=0.001*2*sqrt(pi)*shC;
            end
        end        
        gm_p.eLnc=eLn;
    end
    if mod(itc,200)==0
        face=smt_isolatedFaceRsgn(edge,face,gm_p,rg);  
    end
    itc
end

%% swap face types to attain initial straightness measure. 
for tpc=1:gm_p.nSmt-1
    strtTg=0.65+0.1*2*(rand-0.5);
    strt=1;
    while strt>strtTg
        face=smt_faceTypeSwap(edge,face,tpc,tpc+1,gm_p,rg);
        strt=smt_straightness(vrtx,edge,face,gm_p,tpc,tpc+1);
    end
end

%% run simulations with heterotypic tension every 20*taut
ttMx=mc_p.taus*(gm_p.nSmt+2);
runTm=0;

gm_p.shEd=0.01*2*sqrt(pi);

[itc,exc,dtPt]=deal(1,1,4000);
ttDtPt=ttMx/gm_p.dt/10;    

[vr,ed1,ed2,ed3,ed4,fa1,fa2,fa3]=deal(cell(dtPt,1));
dtc=1;

cs_name=sprintf('somite_mu%s_ften%s_hten%s_tmLag%s',...        
        strrep(num2str(mc_p.mu),'.','d'),...
        strrep(num2str(mc_p.ften),'.','d'),...
        strrep(num2str(mc_p.hten),'.','d'),...
        strrep(num2str(mc_p.tmLag),'.','d'));
    
odir=sprintf('sm_sim/%s/',cs_name);

if ~exist(odir,'dir')
    mkdir(odir);
end

while runTm<ttMx
    
    mc_p.hten_ind=smt_hetTenUpdate(gm_p,mc_p,runTm);
    [vrtx,edge,face]=smt_iteration(vrtx,edge,face,rg,gm_p,mc_p);

    if mod(itc,3)==0
        eLn=smt_edgeLenAll(edge,gm_p);
        t1Id=find(eLn<gm_p.shEd);
        if isempty(t1Id)==0
            for jj=1:size(t1Id,1)
                if gm_p.eLnc(t1Id(jj))>eLn(t1Id(jj))
                    [vrtx,edge,face]=...
                        smt_t1Flip(vrtx,edge,face,rg,gm_p,mc_p,t1Id(jj));                   
                end
            end
            zrId=find(edge{1}(:,4)==0);
            while ~isempty(zrId)
                [vrtx,edge,face,gm_p,~]=...
                    smt_edgeDelete(vrtx,edge,face,rg,gm_p,zrId(1));
                zrId=find(edge{1}(:,4)==0);
            end
        end
        gm_p.eLnc=eLn;
    end
    
    if mod(itc,200)==0
        face=smt_isolatedFaceRsgn(edge,face,gm_p,rg);  
    end
    itc=itc+1;
    runTm=runTm+gm_p.dt;
    
    if mod(itc,10)==3
        [vr{exc},ed1{exc},ed2{exc},ed3{exc},ed4{exc},fa1{exc},...
            fa2{exc},fa3{exc}]=deal(vrtx,edge{1},edge{2},edge{3},...
            edge{4},face{1},face{2},face{3});
        exc=exc+1;
        if exc==dtPt+1
            save(sprintf('%s/%s_vr_r%d_d%d.mat',odir,cs_name,rpc,dtc),'vr');
            save(sprintf('%s/%s_ed1_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed1');
            save(sprintf('%s/%s_ed2_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed2');
            save(sprintf('%s/%s_ed3_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed3');
            save(sprintf('%s/%s_ed4_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed4');
            save(sprintf('%s/%s_fa1_r%d_d%d.mat',odir,cs_name,rpc,dtc),'fa1');
            save(sprintf('%s/%s_fa2_r%d_d%d.mat',odir,cs_name,rpc,dtc),'fa2');
            save(sprintf('%s/%s_fa3_r%d_d%d.mat',odir,cs_name,rpc,dtc),'fa3');
            dtc=dtc+1;
            exc=1;
            ttDtPt=ttDtPt-dtPt;
            [vr,ed1,ed2,ed3,ed4,fa1,fa2,fa3]=...
                deal(cell(min(dtPt,ttDtPt),1));    
        end  
        runTm
    end          
end

if exc>1    
    save(sprintf('%s/%s_vr_r%d_d%d.mat',odir,cs_name,rpc,dtc),'vr');
    save(sprintf('%s/%s_ed1_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed1');
    save(sprintf('%s/%s_ed2_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed2');
    save(sprintf('%s/%s_ed3_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed3');
    save(sprintf('%s/%s_ed4_r%d_d%d.mat',odir,cs_name,rpc,dtc),'ed4');
    save(sprintf('%s/%s_fa1_r%d_d%d.mat',odir,cs_name,rpc,dtc),'fa1');
    save(sprintf('%s/%s_fa2_r%d_d%d.mat',odir,cs_name,rpc,dtc),'fa2');
    save(sprintf('%s/%s_fa3_r%d_d%d.mat',odir,cs_name,rpc,dtc),'fa3');
end  

save(sprintf('%s/%s_gmp_r%d.mat',odir,cs_name,rpc),'gm_p');
save(sprintf('%s/%s_mcp_r%d.mat',odir,cs_name,rpc),'mc_p');            

end