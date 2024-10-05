%% Insert vertex id between vr1 and vr2

function fvrs_n=smt_faceVrIdIn(vr1,vr2,vrin,fvrs)
    fvps=sort([find(fvrs==vr1),find(fvrs==vr2)]);
    if fvps(1)==1 && fvps(2)==size(fvrs,2)
        fvrs_n=[fvrs,vrin];
    else
        fvrs_n=[fvrs(1:fvps(1)),vrin,fvrs(fvps(2):end)];
    end
end