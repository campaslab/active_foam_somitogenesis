%% update heterotypic tension at each time step

function htenInd=smt_hetTenUpdate(gm_p,mc_p,runTm)

% compute a linear function    
% hten_ind=mc_p.hten/mc_p.tauh*(runTm-((1:gm_p.nSmt)+1)*mc_p.taus)+...
%     mc_p.hten/2;

htenInd=mc_p.hten./(1+exp(-(runTm-((1:gm_p.nSmt-1)+1+mc_p.tmLag)*mc_p.taus)*5/mc_p.tauh));

% restrict the range of hten between 0 to hten
for htc=1:size(htenInd,2)
    htenInd(htc)=max(0,htenInd(htc));
    htenInd(htc)=min(mc_p.hten,htenInd(htc));
end

end