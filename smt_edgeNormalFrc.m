%% Compute net normal force for a given edge.
function eNf=smt_edgeNormalFrc(edge,face,rg,eId,del,eln,gm_p)
    
eFa=edge{1}(eId,rg.ei(2):rg.ef(2));
eNf=0;
for ii=1:2
    if eFa(ii)~=0
        faCs=face{1}(eFa(ii),rg.fi(1));
        if faCs<=gm_p.nSmt
            eNf=eNf+del*(1/(face{1}(eFa(ii),rg.fi(2)))-1)...
                *eln/2*((-1)^(ii));
        else
            eNf=eNf+del*(1/(face{1}(eFa(ii),rg.fi(2))/2)-1)...
                *eln/2*((-1)^(ii));
        end
    end
end

end