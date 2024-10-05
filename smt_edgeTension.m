%% Compute fixed tension point

function [etn,etp]=smt_edgeTension(face,edge,rg,mc_p,gm_p,eId)

% find edge-face id and delete 0, find face types. 
efId=edge{1}(eId,rg.ei(2):rg.ef(2));
efId=efId(efId~=0);
efTp=face{1}(efId,rg.fi(1));

if size(efId,2)==2 && max(efTp)<=gm_p.nSmt
    % compute homegeneous tension of cell-cell junctions
    etn=1;
    etp=2;
    
    % compute heterogeneous tension of cell-cell junctions
    if efTp(1)~=efTp(2) 
        efTp=sort(efTp);
        etn=etn+mc_p.hten_ind(efTp(1));
%         etn=etn+mc_p.hten_ind(efTp(1))+mc_p.hten_ind(efTp(2));
    end
else
    % compute homegeneous tension of cell-extracellular space junctions
    etn=mc_p.ften;    
    etp=1;
end

end