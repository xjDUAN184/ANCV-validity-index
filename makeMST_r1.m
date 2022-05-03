function T = makeMST_r1(dist)
    xishu = sparse(dist) ; 
    G = prim_mst(xishu) ; 
    T = full(G) ; 
end

