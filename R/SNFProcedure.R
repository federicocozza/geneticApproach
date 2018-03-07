################################################################################
SNFProcedure = function(list_dat, k, t, a,NUMC = 2:5) {
    # Normalize the features in each of the views
    dataL = lapply(list_dat, standardNormalization)
    # Calculate the distances for each view
    distL = lapply(dataL, function(x) SNFtool::dist2(x, x))
    # Construct the similarity graphs
    affinityL = lapply(distL, function(x) affinityMatrix(x, K = k, sigma = a))
    # Construct the fused network
    W = SNF(affinityL, K = k, t = t)


    nc_tip = estimateNumberOfClustersGivenGraph(W, NUMC)

    return(list(W = W, AF = affinityL, nc_tip = nc_tip))
}
