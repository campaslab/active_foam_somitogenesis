% somiteFn(len,wth,nSmt,edpc,ften,beta,taut,mu,del,hten,taus,tauh,tmLag,rpc)

tic;
[vrtx,edge,face,gm_p,mc_p,rg]=...
    smt_somite(4,8,3,0.15,2,0,10,0.5,10,1,20,40,0.3,1);
toc;
