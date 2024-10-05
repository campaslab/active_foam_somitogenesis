%% Construct initial configuration using Lloyd algorithm
function sd=smt_seedLloyd(bs,sdPar)

% first generate random seed points inside a rectangle
nFa=bs^2;
sd=rand(nFa,2)*bs;

% displace center position to centroids by sd_parameter times
for rpc=1:sdPar
    sd_copy=smt_seedCopy(sd,bs);
    [vrtx,v_cell]=voronoin(sd_copy);
    for fac=1:nFa
        v_cell{fac}=smt_vrtxIdSort(v_cell{fac},vrtx);
        sd(fac,:)=smt_faceCentroid(vrtx(v_cell{fac},:));
    end
    sd=mod(sd,bs);
end    

end