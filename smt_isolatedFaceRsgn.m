%% check connectivity of somite and chnage cell type for isolated cells
function face_n=smt_isolatedFaceRsgn(edge,face,gm_p,rg)

% define output variable
face_n=face;


for smc=1:gm_p.nSmt
    
    faId=find(face{1}(:,1)==smc);
    faAggrt=cell(size(faId,1),1);
    agc=1;
    while isempty(faId)==0
        if isempty(faAggrt{agc})==1
            faAggrt{agc}=faId(1);
            faId=faId(2:end);
        else
            faAggrtNew=cell(size(faAggrt{agc},1),1);
            for fac=1:size(faAggrt{agc},1)
                faAggrtNew{fac}=unique(edge{1}(abs(face{3}{faAggrt{agc}(fac)}),3:4));                
            end
            faAggrtNew=unique(cell2mat(faAggrtNew));
            faAggrtNew=faAggrtNew(faAggrtNew~=0);
            faAggrtNew=faAggrtNew(face{1}(faAggrtNew,1)==smc);
            if isequal(faAggrt{agc},faAggrtNew)
                agc=agc+1;
            else
                faAggrt{agc}=faAggrtNew;
            end
            faId=setdiff(faId,faAggrtNew);
        end
    end
    
    if agc>1
        aggSz=zeros(agc,1);
        for ii=1:agc
            aggSz(ii)=size(faAggrt{ii},1);
        end
        [~,idx]=sort(aggSz);
        for ii=1:agc-1
            faRsgnId=faAggrt{idx(ii)};
            faNeId=cell(size(faRsgnId,1),1);
            for fac=1:size(faRsgnId,1)
                faNeId{fac}=unique(edge{1}(abs(face{3}{faRsgnId(fac)}),3:4));                
            end
            faNeId=unique(cell2mat(faNeId));
            faNeId=faNeId(faNeId~=0);
            faNeId=faNeId(face{1}(faNeId,1)~=smc);
            face_n{1}(faRsgnId,1)=face{1}(faNeId(1),1);
        end
    end
%     faAggrt{1}
end
    

% 
% % check whether there is any isolated cells and change cell type
% for fac=1:gm_p.nFa
%     fn_id=unique(edge{1}(abs(face_n{3}{fac}),rg.ei(2):rg.ef(2)));
%     fn_id=fn_id(fn_id~=fac & fn_id~=0);
%     fn_tp=face_n{1}(fn_id,1);
%     if any(fn_tp==face_n{1}(fac))==0
%         face_n{1}(fac)=fn_tp(1);
%     end
% end

end