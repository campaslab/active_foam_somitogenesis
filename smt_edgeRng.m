%% Compute edge range.
% range is [x_min, x_max, y_min, y_max]
function erng=smt_edgeRng(edge,eId)

emd=edge{2}{eId};
erng=[min(emd(:,1)),max(emd(:,1)),min(emd(:,2)),max(emd(:,2))];

end