# Run first
#########################################################################################################################################
# Packages

# read packages from local directory
.libPaths("~/R_packages/")

req_packages = c(
  'cetcolor',
  'combinat',
  'data.table',
  'devtools',
  'extrafont',
  'gdata',
  'ggplot2',
  'ggtern',
  'grid',
  'gtools',
  'igraph',
  'jcolors',
  'LICORS',
  'lwgeom',
  'Matrix',
  'mclust',
  'PRIMME',
  'R.utils',
  'RColorBrewer',
  'rgeos',
  'rmapshaper',
  'rnaturalearth',
  'rstudioapi',
  'RSpectra',
  'sf',
  'testthat',
  'USAboundaries',
  'USAboundariesData'
)

new_packages = req_packages[!(req_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)){
  cat("PLEASE INSTALL THE FOLLOWING PACKAGES:\n")
  cat(new_packages, sep="\n")
  cat("\n")
}

for(package in req_packages){
  require(package, character.only = TRUE)
}

#loadfonts()

# Ghostscript
#Sys.setenv(R_GSCMD = "/usr/bin/gs")


# Clear workspace
rm(list=ls())
setwd("~/Documents/github/masters-dissertation/r/")

# Ghostscript
#Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.26/bin/gswin64c.exe")

# Generative Models
#########################################################################################################################################
sampleDSBM = function(block_sizes, sparsity_matrix){

  n = sum(block_sizes)
  k = length(block_sizes)
  cumul_sizes = c(0,cumsum(block_sizes))
  G = matrix(0, nrow=n, ncol=n)

  for(i in 1:k){
    for(j in 1:k){
      x_range = (cumul_sizes[i]+1):cumul_sizes[i+1]
      y_range = (cumul_sizes[j]+1):cumul_sizes[j+1]
      G[x_range, y_range] = rbinom(block_sizes[i]*block_sizes[j], 1, sparsity_matrix[i,j])
    }
  }

  diag(G) = 0

  return(drop0(G))
}

sampleDirectedCP = function(block_sizes, p1, p2){

  sparsity_matrix = matrix(c(p2,p1,p2,p2,
                             p2,p1,p2,p2,
                             p2,p1,p1,p1,
                             p2,p2,p2,p2), nrow=4, byrow=TRUE)

  G = sampleDSBM(block_sizes, sparsity_matrix)

  return(G)
}

sampleBSBM = function(dest_block_sizes, targ_block_sizes, p1, p2){

  sparsity_matrix = matrix(c(0,0,p1,p2,
                             0,0,p2,p1,
                             0,0,0,0,
                             0,0,0,0), nrow=4, byrow=TRUE)

  G = sampleDSBM(c(dest_block_sizes, targ_block_sizes), sparsity_matrix)

  return(G)
}


# Implementations of Other Methods
#########################################################################################################################################
runSVDEmb = function(G){

  # TODO check correct

  #n = nrow(G)
  #rows_to_keep = (1:n)[(apply(G,1,sum)>0)]
  #cols_to_keep = (1:n)[(apply(G,2,sum)>0)]
  #A = G[rows_to_keep, cols_to_keep]
  #invD1 = diag(apply(A,1,sum) ^ (-0.5))
  #invD2 = diag(apply(A,2,sum) ^ (-0.5))

  #An = invD1 %*% A %*% invD2
  #AnSVD = svds(An, 2)
  #u2 = AnSVD$u[,1]
  #v2 = AnSVD$v[,1]
  #z2 = c(invD1 %*% u2, invD2 %*% v2)

  z2 = runLapEmb(G+t(G),2,'rw')$vects[,2]

  return(z2)

}

runHermRW = function(G, l){

  # TODO if needed

  n = nrow(G)
  A = as.matrix(G - t(G)) * 1i
  Dinv = diag(apply(Mod(A),1,sum) ^ (-1))

  spect = eigen(Dinv %*% A)
  inds = 2*(1:k)

  vects = spect$vectors[,inds]
  print(vects)

  #P =
}


# Motif adjacency matrices
#########################################################################################################################################
motifAdjacency = function(G, motif_name, type=c('func','struc')){

  # Builds motif adjacency matrix for a simple graph adjacency matrix G and motif name.
  # type is either: 'func'  for finding all instances of S in G.
  #                 'struc' for finding all induced instances of S in G

  type = match.arg(type)
  IM = buildIndMats(G)

  if(type=='func'){

    motifadj = motifAdjCalcs(G, G, IM$Gd, IM$J, IM$Jn, IM$J, IM$Jd, motif_name)
  }

  else if(type=='struc'){

    motifadj = motifAdjCalcs(G, IM$Gs, IM$Gd, IM$J, IM$J0, IM$Js, IM$Jd, motif_name)
  }

  motifadj = unname(drop0(motifadj))

  return(motifadj)
}

buildIndMats = function(G){

  J  = drop0_killdiag( 1*(G > 0) )
  J0 = drop0_killdiag( 1*((G+t(G)) == 0 ) )
  Jn = drop0_killdiag( 1+0*G )
  Gs = drop0_killdiag( G*(1 - t(J)) )
  Gd = drop0_killdiag( (G+t(G)) * J * t(J) )
  Js = drop0_killdiag( 1*(Gs > 0) )
  Jd = drop0_killdiag( 1*(Gd > 0) )

  return(list(J=J, J0=J0, Jn=Jn, Gs=Gs, Gd=Gd, Js=Js, Jd=Jd))
}

motifAdjCalcs = function(G, Gs, Gd, J, J0, Js, Jd, motif_name){

  # TODO check args - G, J not needed
  # TODO make func by default

  if(motif_name == 'Ms'){
    motifadj = Gs + t(Gs)
  }

  else if(motif_name == 'Md'){
    motifadj = Gd / 2
  }

  else if(motif_name == 'M1'){
    C = t(Js)*(Js%*%Gs) + t(Js)*(Gs%*%Js) + t(Gs)*(Js%*%Js)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M2'){
    C = t(Js)*(Jd%*%Gs) + t(Js)*(Gd%*%Js) + t(Gs)*(Jd%*%Js)
    C = C + t(Js)*(Js%*%Gd) + t(Js)*(Gs%*%Jd) + t(Gs)*(Js%*%Jd)
    C = C + Jd*(Js%*%Gs) + Jd*(Gs%*%Js) + Gd*(Js%*%Js)
    motifadj = (C + t(C)) / 4
  }

  else if(motif_name == 'M3'){
    C = Js*(Jd%*%Gd) + Js*(Gd%*%Jd) + Gs*(Jd%*%Jd)
    C = C + Jd*(Jd%*%Gs) + Jd*(Gd%*%Js) + Gd*(Jd%*%Js)
    C = C + Jd*(Js%*%Gd) + Jd*(Gs%*%Jd) + Gd*(Js%*%Jd)
    motifadj = (C + t(C)) / 5
  }

  else if(motif_name == 'M4'){
    motifadj = (Jd*(Jd%*%Gd) + Jd*(Gd%*%Jd) + Gd*(Jd%*%Jd)) / 6
  }

  else if(motif_name == 'M5'){
    C = Js*(Js%*%Gs) + Js*(Gs%*%Js) + Gs*(Js%*%Js)
    C = C + Js*(Js%*%t(Gs)) + Js*(Gs%*%t(Js)) + Gs*(Js%*%t(Js))
    C = C + Js*(t(Js)%*%Gs) + Js*(t(Gs)%*%Js) + Gs*(t(Js)%*%Js)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M6'){
    C = Js*(Js%*%Gd) + Js*(Gs%*%Jd) + Gs*(Js%*%Jd)
    Cprime = Jd*(t(Js)%*%Gs) + Jd*(t(Gs)%*%Js) + Gd*(t(Js)%*%Js)
    motifadj = (C + t(C) + Cprime) / 4
  }

  else if(motif_name == 'M7'){
    C = Js*(Jd%*%Gs) + Js*(Gd%*%Js) + Gs*(Jd%*%Js)
    Cprime = Jd*(Js%*%t(Gs)) + Jd*(Gs%*%t(Js)) + Gd*(Js%*%t(Js))
    motifadj = (C + t(C) + Cprime) / 4
  }

  else if(motif_name == 'M8'){
    C = Js*(Gs%*%J0) + Gs*(Js%*%J0)
    Cprime = J0*(t(Js)%*%Gs) + J0*(t(Gs)%*%Js)
    motifadj = (C + t(C) + Cprime) / 2
  }

  else if(motif_name == 'M9'){
    C = Js*(J0%*%t(Gs)) + Gs*(J0%*%t(Js))
    C = C + J0*(Js%*%Gs) + J0*(Gs%*%Js)
    C = C + Js*(t(Gs)%*%J0) + Gs*(t(Js)%*%J0)
    motifadj = (C + t(C)) / 2
  }

  else if(motif_name == 'M10'){
    C = Js*(J0%*%Gs) + Gs*(J0%*%Js)
    Cprime = J0*(Js%*%t(Gs)) + J0*(Gs%*%t(Js))
    motifadj = (C + t(C) + Cprime) / 2
  }

  else if(motif_name == 'M11'){
    C = Jd*(Gs%*%J0) + Gd*(Js%*%J0)
    C = C + J0*(Jd%*%Gs) + J0*(Gd%*%Js)
    C = C + Js*(Gd%*%J0) + Gs*(Jd%*%J0)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M12'){
    C = Jd*(J0%*%Gs) + Gd*(J0%*%Js)
    C = C + J0*(Js%*%Gd) + J0*(Gs%*%Jd)
    C = C + Js*(J0%*%Gd) + Gs*(J0%*%Jd)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M13'){
    C = Jd*(Gd%*%J0) + Gd*(Jd%*%J0) + J0*(Jd%*%Gd)
    motifadj = (C + t(C)) / 4
  }

  else if(motif_name == 'coll'){
    C = J0*(Js%*%t(Gs)) + J0*(Gs%*%t(Js))
    motifadj = C / 2
  }

  else if(motif_name == 'expa'){
    C = J0*(t(Js)%*%Gs) + J0*(t(Gs)%*%Js)
    motifadj = C / 2
  }

  else if(motif_name == 'path'){
    C = J0*(Js%*%Gs) + J0*(Gs%*%Js)
    motifadj = (C + t(C)) / 2
  }


  return(motifadj)

}


# Clustering
#########################################################################################################################################
scorePartition = function(G, partLabels, scoreType=c('cut_size','ratio_cut','n_cut','conductance'), directed=FALSE, warnings=TRUE){

  # Returns a score for a given k-partitioned graph adjacency matrix.
  # partLabels are of the form [1,1,2,1,2,3,2] for example.

  scoreType = match.arg(scoreType)
  score = 0
  uLabels = unique(partLabels)
  n = nrow(G)

  if(length(uLabels) <= 1){
    stop('Must be at least two clusters')
  }

  if(scoreType == 'cut_size'){
    for(k in uLabels){
      score = score + sum(G[partLabels == k, partLabels != k])
    }
    if(!directed){
      score = score / 2
    }
  }

  else if(scoreType == 'ratio_cut'){
    for(k in uLabels){
      score = score + sum(G[partLabels == k, partLabels != k]) / sum(partLabels == k)
    }
    if(!directed){
      score = score / 2
    }
  }

  else if(scoreType == 'n_cut'){
    for(k in uLabels){
      volK = sum(G[partLabels == k,])
      if(volK == 0){
        if(warnings){
          warning('Cluster has volume zero')
        }
        return(Inf)
      }
      score = score + sum(G[partLabels == k, partLabels != k]) / volK
    }
    score = score/2
  }

  else if(scoreType == 'conductance'){
    for(k in uLabels){
      volK = sum(G[partLabels == k, partLabels == k])
      volKbar = sum(G[partLabels != k, partLabels != k])
      if(volK == 0 || volKbar == 0){
        if(warnings){
          warning('Cluster has volume zero')
        }
        return(Inf)
      }
      score = score + sum(G[partLabels == k, partLabels != k]) / min(volK, volKbar)
    }
  }

  return(score)

}

sweepScores = function(G, vector, scoreType=c('cut_size','ratio_cut','n_cut','conductance'), simplify=1, verbose=TRUE){

  # Returns the vector of sweep scores given an eigenvector of embedding values.

  scoreType = match.arg(scoreType)

  n = length(vector)
  sorted = sort(vector)
  scores = rep(0,n-1)
  simplify_iterator = seq(1, n-1, simplify)

  for(i in simplify_iterator){
    partLabels = rep(1,n)
    cutoff = sorted[i]
    partLabels[vector > cutoff] = 2
    scores[i] = scorePartition(G, partLabels, scoreType, warnings=FALSE)

    if(verbose){
      cat('\t\tSweep profile ',round(100*i/(n-1)), '% done         \r', sep='')
    }
  }

  if(verbose){
    cat('\n')
  }

  return(scores)
}

minSweepPartition = function(G, vector, scoreType=c('cut_size','ratio_cut','n_cut','conductance')){

  # Returns a 2-part partition (in labels vector notation) which minimises the sweep score.

  scoreType = match.arg(scoreType)

  scores = sweepScores(G, vector, scoreType)
  sorted = sort(vector)
  cutoff = sorted[which.min(scores)]
  partLabels = rep(1, length(vector))
  partLabels[vector > cutoff] = 2

  return(partLabels)

}

minSweepPartition2 = function(G, vector, scores, scoreType=c('cut_size','ratio_cut','n_cut','conductance')){

  # Returns a 2-part partition (in labels vector notation) which minimises the sweep score.
  # For use when scores are already known

  scoreType = match.arg(scoreType)

  sorted = sort(vector)
  cutoff = sorted[which.min(scores)]
  partLabels = rep(1, length(vector))
  partLabels[vector > cutoff] = 2

  return(partLabels)

}

sortMatrix = function(mat, clusts){

  # TODO remove

  # get parameters
  ans = list()
  u_clusts = unique(clusts)
  k = length(u_clusts)
  n = nrow(mat)

  # randomise inputs
  samp = sample(1:n, n, replace=FALSE)
  mat = mat[samp,samp]
  clusts = clusts[samp]

  i = 1
  for(perm in permn(u_clusts)){

    perm_func = function(x){return(perm[x])}
    p_clusts = sapply(clusts, perm_func)
    inds = order(p_clusts)
    p_matrix = mat[inds,inds]
    ans[[i]] = p_matrix
    i = i+1
  }

  return(ans)

}

clustDirectedCP = function(G, clusts){

  n = nrow(G)
  errors = c()
  perms = permn(c(1,2,3,4))

  for(perm in perms){

    block_sizes = c(sum(clusts == perm[1]),sum(clusts == perm[2]),sum(clusts == perm[3]),sum(clusts == perm[4]))
    ord = c((1:n)[clusts == perm[1]],(1:n)[clusts == perm[2]],(1:n)[clusts == perm[3]],(1:n)[clusts == perm[4]])
    truth = sampleDirectedCP(block_sizes,1,0)
    error = sum(abs(1*(G[ord,ord]>0) - truth))
    errors = c(errors, error)
  }

  perm = unlist(perms[which.min(errors)])
  newClusts = rep(0,n)

  for(i in 1:4){
    newClusts[clusts == perm[i]] = i
  }

  return(newClusts)
}

makeClusterLike = function(clusts, target_clusts){

  n = length(clusts)
  errors = c()
  u_clusts = unique(clusts)
  perms = permn(u_clusts)

  for(perm in perms){

    newClusts = rep(0,n)

    for(i in u_clusts){
      newClusts[clusts == i] = perm[i]
    }

    error = sum(target_clusts != newClusts)
    errors = c(errors, error)
  }

  perm = unlist(perms[which.min(errors)])
  newClusts = rep(0,n)

  for(i in u_clusts){
    newClusts[clusts == i] = perm[i]
  }

  return(newClusts)

}


# Helper functions
#########################################################################################################################################
dataFrameNames = function(column_names, n_rows){

  df = as.data.frame(matrix(NA, nrow=n_rows, ncol=length(column_names)))
  names(df) = column_names

  return(df)
}

drop0_killdiag = function(mat){

  ans = mat
  diag(ans) = 0
  ans = drop0(ans)

  return(ans)
}

isDirectedMotif = function(motif_name){

  if(motif_name %in% undirectedMotifs()){
    return(FALSE)
  }

  else if(motif_name %in% directedMotifs()){
    return(TRUE)
  }

  else{
    stop('Invalid motif name')
  }
}

largestComponent = function(G){

  n = nrow(G)
  Gr = graph_from_adjacency_matrix(G, weighted=TRUE)
  comps = components(Gr)
  verts_to_keep = (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)

}

buildMotif = function(motif_name){

  # Returns common motifs as graphs.

  if(motif_name == 'Ms')  {S = graph.formula(1-+2)}
  else if(motif_name == 'Md') {S = graph.formula(1++2)}
  else if(motif_name == 'M1')  {S = graph.formula(1-+2-+3-+1)}
  else if(motif_name == 'M2')  {S = graph.formula(1-+2++3-+1)}
  else if(motif_name == 'M3')  {S = graph.formula(1-+2++3++1)}
  else if(motif_name == 'M4')  {S = graph.formula(1++2++3++1)}
  else if(motif_name == 'M5')  {S = graph.formula(1-+2-+3+-1)}
  else if(motif_name == 'M6')  {S = graph.formula(1+-2-+3++1)}
  else if(motif_name == 'M7')  {S = graph.formula(1-+2+-3++1)}
  else if(motif_name == 'M8')  {S = graph.formula(1+-2-+3)}
  else if(motif_name == 'M9')  {S = graph.formula(1+-2+-3)}
  else if(motif_name == 'M10') {S = graph.formula(1-+2+-3)}
  else if(motif_name == 'M11') {S = graph.formula(1++2-+3)}
  else if(motif_name == 'M12') {S = graph.formula(1++2+-3)}
  else if(motif_name == 'M13') {S = graph.formula(1++2++3)}

  return(S)

}

adj_from_edge_list = function(edge_list){

  n_edges = nrow(edge_list)
  n_col = ncol(edge_list)
  vertices = unique(c(edge_list[,1], edge_list[,2]))
  n_vertices = length(vertices)
  adj = matrix(0, nrow=n_vertices, ncol=n_vertices)

  # create edge weights
  if(n_col == 2){
    edge_weights = rep(1, n_edges)
  }
  else if(n_col == 3){
    edge_weights = edge_list[,3]
  }
  else{
    stop('edge_list must have 2 or 3 columns')
  }

  # build adj
  for(i in 1:n_edges){
    adj[edge_list[i,1],edge_list[i,2]] = adj[edge_list[i,1],edge_list[i,2]] + edge_weights[i]
  }

  return(drop0_killdiag(adj))
}

scaleVec = function(x){

  temp = (x - min(x)) / (max(x) - min(x))
  return(temp)
}

simplify_sf = function(sf_object, tol){

  # TODO remove

  temp = methods::as(object = sf_object, Class = "Spatial")
  temp = gSimplify(temp, tol, topologyPreserve = TRUE)
  temp = st_as_sf(temp)

  return(temp)
}

sigfig = function(vec, digits){

  return(gsub("\\.$", "", formatC(signif(vec,digits=digits), digits=digits, format="fg", flag="#")))
}

# Spectral methods
#########################################################################################################################################
computeTopSpectrum = function(L, typeLap=c('comb','rw'), topk){

  # Computes eigenvectors/values of lowest 'topk' eigenvalues of a Laplacian L.
  # Returns a list of $vects and $vals

  if(typeLap == 'comb'){
    ansEigs = eigs_sym(L, topk, which = 'SM')
  }
  else if(typeLap == 'rw'){
    ansEigs = eigs(L, topk, which = 'SM')
  }

  inds = seq(topk,1,-1)
  vects = Re(ansEigs[['vectors']])[,inds]
  vals = Re(ansEigs[['values']])[inds]

  ansSpect = list()
  ansSpect[['vects']] = vects
  ansSpect[['vals']]  = vals

  return(ansSpect)

}

buildLaplacian = function(G, typeLap=c('comb','rw')){

  # Builds various types of Laplacian Matrix, given an adjacency matrix G and a type of Laplacian.

  # typeLap
  # 'comb'      combinatorial Laplacian                           L = D - G
  # 'rw'        random walk Laplacian (row-normalised)            L = D^(-1) (D - G)

  typeLap = match.arg(typeLap)

  degsG = apply(G, 1, sum)
  n = nrow(G)

  if (typeLap=='comb'){
    D = diag(degsG)
    L =  D - G
  }
  else if (typeLap == 'rw'){
    Dinv = diag(degsG^(-1))
    L =  diag(n) - Dinv %*% G
  }

  return(L)
}

runLapEmb = function(G, topk, typeLap=c('comb','rw')){

  # Runs Laplace embedding of an adjacency matrix G, given number of clusters and
  # Laplacian type. Returns a list with $vects and $vals.

  typeLap = match.arg(typeLap)

  L = buildLaplacian(G, typeLap)
  ansSpect = computeTopSpectrum(L, typeLap, topk)

  return(ansSpect)

}

runMotifEmb = function(G, motif_name, type, typeLap, numEigs){

  M = motifAdjacency(G, motif_name, type)
  comps = largestComponent(M)
  M = M[comps,comps, drop=FALSE]
  Gcomps = G[comps,comps, drop=FALSE]
  spect = runLapEmb(M,numEigs,typeLap)
  ans = list()
  ans$vects = spect$vects
  ans$vals = spect$vals
  ans$comps = comps
  ans$M = M
  ans$Gcomps = Gcomps

  return(ans)
}


# Plotting
#########################################################################################################################################
theme_diss = function(){

  # TODO more in here

  temp = theme_bw() +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(plot.margin = margin(t=0, r=10, b=0, l=0)) +
    theme(text = element_text(family='CM Roman')) +
    theme(axis.text = element_text(color='black',size=20)) + # size=15 before
    theme(axis.title = element_text(size=20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


  return(temp)
}


# Experiments
#########################################################################################################################################
experimentRunAll = function(){

  t00 = Sys.time()

  # Eigenvector sweep
  experimentEigenvectorSweep()

  # Symmetric DSBM
  ns = c(50,100)
  ps = c(0.3,0.15)
  qs = c(0.2,0.1)
  n_repeats = 20
  set.seed(6552)
  experimentMotif(ns, ps, qs, qs, ps, n_repeats, 'sym')

  # Asymmetric DSBM
  ns = c(100,200)
  ps = c(0.2,0.15)
  q1s = c(0.35,0.25)
  q2s = c(0.05,0.05)
  n_repeats = 20
  set.seed(3193)
  experimentMotif(ns, ps, q1s, q2s, ps, n_repeats, 'assym')

  # US Political blogs network
  experimentPolblogs()

  # US Migration network
  motif_names = c('Ms','M6','M9')
  eigenvectors_for_map = c(2:7)
  n_vects_for_cluster = 7
  n_clusters = 7
  simplify_sweep = 1
  experimentMigration(motif_names, eigenvectors_for_map, n_vects_for_cluster, n_clusters, simplify_sweep)
  us_migration_table(motif_names, eigenvectors_for_map)

  # BSBM
  ns = c(100,200)
  ps = c(0.2,0.1)
  qs = c(0.1,0.06)
  n_repeats = 20
  set.seed(6342)
  experimentBipartite(ns, ps, qs, n_repeats)

  # American Revolution network
  experimentAmericanRevolution()

  # Unicode Languages network
  eigenvectors_for_map = 2:6
  eigenvectors_for_cluster = 2:6
  n_clusters_to_print = 6
  experimentLanguages(eigenvectors_for_map, eigenvectors_for_cluster, n_clusters_to_print)

  # Timings for MAM computation
  experimentTiming()

  cat('All experiments run!\n')
  cat(capture.output(Sys.time()-t00), '\n\n')
  return()
}

experimentMotif = function(ns, p1s, p2s, p3s, p4s, n_repeats, experiment_name){

  if(experiment_name == 'sym'){
    cat('Symmetric DSBM\n')
  }

  if(experiment_name == 'assym'){
    cat('Asymmetric DSBM\n')
  }

  t0 = Sys.time()
  n_experiments = length(ns)
  motif_names = c('Ms','Md','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','M13')
  plot_motif_names = c('M[s]','M[d]','M[1]','M[2]','M[3]','M[4]','M[5]','M[6]','M[7]','M[8]','M[9]','M[10]','M[11]','M[12]','M[13]')

  for(expmt in 1:n_experiments){

    # prepare experiment parameters
    n = ns[expmt]
    p1 = p1s[expmt]
    p2 = p2s[expmt]
    p3 = p3s[expmt]
    p4 = p4s[expmt]
    block_sizes = rep(n,2)
    sparsity_matrix = matrix(c(p1,p2,p3,p4), nrow=2, byrow=TRUE)
    truth = c(rep(1,n),rep(2,n))

    # prepare results table and row counter
    resultsColNames = c('Motif','ARI','CompSize')
    resultsNRows = n_repeats * length(motif_names)
    results = dataFrameNames(resultsColNames, resultsNRows)
    row_number = 1

    for(rep in 1:n_repeats){

      # sample the graph
      G = sampleDSBM(block_sizes, sparsity_matrix)

      for(motif_name in motif_names){

        # run algorithm
        motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 2)
        comps = motifEmb$comps
        vects = motifEmb$vects
        clusts = kmeanspp(vects[,2],k=2)$cluster
        randIndex = adjustedRandIndex(clusts,truth[comps])

        # enter results and increment counter
        results[row_number,1] = motif_name
        results[row_number,2] = randIndex
        results[row_number,3] = length(comps)
        row_number = row_number + 1

        # progress
        cat('\tExperiment ',expmt,'/',n_experiments,', repeat ',rep,'/',n_repeats,', motif ',motif_name,'\t\t\t\r',sep='')
      }
    }

    # generate mean connected component sizes
    meanCompSize = data.frame(meanCompSize = round(c(by(results$CompSize, results$Motif, mean))[rank(motif_names)], 0))

    # generate plot
    pl = ggplot() +
      geom_abline(slope=0, intercept=0, linetype='dashed', color='gray50') +
      geom_abline(slope=0, intercept=1, linetype='dashed', color='gray50') +
      geom_violin(data=results, mapping=aes(x=Motif, y=ARI), scale='width', color='black', fill='gray90', width=0.7) +
      stat_summary(data=results, mapping=aes(x=Motif, y=ARI), fun.y=median, geom='point', size=1.5) +
      stat_summary(data=results, mapping=aes(x=Motif, y=ARI), fun.y=function(z){10}, fun.ymin=function(z){quantile(z,0.25)}, fun.ymax=function(z){quantile(z,0.75)}, size=0.5) +
      scale_x_discrete(limits=motif_names, labels=parse(text=plot_motif_names)) +
      xlab('Motif') +
      ylab('ARI') +
      coord_cartesian(ylim = c(-0.5,1), clip='off') +
      geom_text(data=meanCompSize, mapping=aes(x=1:15, y=1.15, label=meanCompSize), family='CM Roman',size=5)

    # legend on first plot only
    if((experiment_name == 'sym') & (expmt == 1)){
      pl = pl +
        annotate(geom='rect', xmin=10.5, xmax=15.5, ymin=-0.55, ymax=-0.02, fill='white', color='white') +
        geom_segment(data=data.frame(x=11, xend=11, y=-0.07, yend=-0.5), aes(x=x, xend=xend, y=y, yend=yend)) +
        geom_point(data=data.frame(x=11, y=0.5*(-0.07-0.5)), aes(x=x, y=y)) +
        annotate(geom='text', x=13.3, y=-0.1, label='Upper quartile', family='CM Roman', size=6) +
        annotate(geom='text', x=12.4, y=0.5*(-0.07-0.5), label='Median', family='CM Roman',size=6) +
        annotate(geom='text', x=13.3, y=-0.47, label='Lower quartile', family='CM Roman',size=6)
    }

    # pdf name
    if(experiment_name == 'sym'){
      pdf_name = paste('../results/motifsym/motifsym_',expmt,'.pdf',sep='')
    }

    else if(experiment_name == 'assym'){
      pdf_name = paste('../results/motifassym/motifassym_',expmt,'.pdf',sep='')
    }

    pl = pl + (theme_diss() +
                 theme(plot.title = element_text(vjust=7)) +
                 theme(plot.margin = margin(t=20, r=10, b=0, l=0)) +
                 theme(axis.title.y = element_text(margin=margin(r=-12))) +
                 theme(axis.text.x = element_text(size=14))) +
      geom_text(mapping=aes(x=1, y=1.145, label='"|"*italic(C)*"| = "'), parse=TRUE, family='CM Roman',size=5,hjust=1.3)

    # print to pdf
    ggsave(filename=pdf_name, plot=pl, width=7, height=5)
    embed_fonts(pdf_name, outfile=pdf_name)
  }

  cat('\n')
  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')
  return()
}

experimentPolblogs = function(){

  cat('US Political Blogs network\n')

  t0 = Sys.time()

  # read data
  graph_data = readGraphData('polblogs')
  G = graph_data$adj
  truth = graph_data$labels+1
  motif_names = c('Ms','Md','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','M13')
  plot_motif_names = c('M[s]','M[d]','M[1]','M[2]','M[3]','M[4]','M[5]','M[6]','M[7]','M[8]','M[9]','M[10]','M[11]','M[12]','M[13]')

  # run alg so points can be circled
  comps = largestComponent(G)
  cat('\tNumber of liberal blogs:', sum(truth[comps] == 1), '\n')
  cat('\tNumber of conservative blogs:', sum(truth[comps] == 2), '\n')
  cat('\tNumber of edges:', sum(G[comps,comps] > 0), '\n')
  set.seed(9246)
  motifEmb = runMotifEmb(G, 'Ms', 'func', 'rw', 3)
  comps = motifEmb$comps
  vects = motifEmb$vects
  clusts = kmeanspp(vects[,2],k=2)$cluster
  outliers = (1:length(comps))[clusts == 2]

  # plot network
  G_graph = graph_from_adjacency_matrix(G[comps,comps], weighted=TRUE)
  set.seed(407)
  plotcoord = data.frame(layout_with_fr(G_graph, grid='grid'))
  plotcoord$cols = as.character(truth[comps])
  colnames(plotcoord) = c("X1","X2",'cols')
  edgelist = get.edgelist(G_graph)
  edges = data.frame(plotcoord[edgelist[,1],1:2], plotcoord[edgelist[,2],1:2])
  colnames(edges) = c("X1","Y1","X2","Y2")
  degrees = apply(G[comps,comps], 1, sum) + apply(G[comps,comps], 2, sum)
  degrees = 0.7*(degrees)^0.3
  outliers = data.frame(x=plotcoord$X1[outliers], y=plotcoord$X2[outliers])
  dividing_line = data.frame(x=c(0.15,0.37), y=c(-0.85,-0.5))

  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='grey', size=0.1) +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord, size=degrees, alpha=0.4) +
    geom_point(data=outliers, aes(x=x,y=y), shape=1, size=3, stroke=1) +
    geom_line(data=dividing_line, aes(x=x,y=y), linetype='dashed') +
    scale_color_manual(labels=c('Liberal','Conservative'), values=c('dodgerblue','firebrick2')) +
    ggtitle('') +
    xlab('') +
    ylab('') +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank(), plot.margin = margin(t=-10,b=-20,r=10,l=-20),
                          legend.position = c(0.83,0.92), legend.title=element_blank(), legend.text=element_text(size=20),
                          panel.border = element_blank(), panel.grid = element_blank()))

  # print network plot pdf
  pdf_name = '../results/polblogs/polblogs_network.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)


  # prepare results table and row counter
  resultsColNames = c('Motif','ARI','CompSize')
  resultsNRows = length(motif_names)
  results = dataFrameNames(resultsColNames, resultsNRows)
  row_number = 1

  for(motif_name in motif_names){

    # run algorithm
    motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 3)
    comps = motifEmb$comps
    vects = motifEmb$vects
    clusts = kmeanspp(vects[,2],k=2)$cluster
    randIndex = adjustedRandIndex(clusts,truth[comps])

    # enter results and increment counter
    results[row_number,1] = motif_name
    results[row_number,2] = randIndex
    results[row_number,3] = length(comps)
    row_number = row_number + 1

    # progress
    cat('\tMotif ',motif_name,', ARI = ',round(randIndex,2),', |C| = ',length(comps),'\t\t\t\n',sep='')
  }

  # process results for ari/comp plot
  hjusts=c(0.5, 1.1, 1.2, -0.3, 1.2, 0.5, -0.1, -0.1, 0.4, 0.7, 0.2, 0.5, 0.7, 0.5, -0.2)
  vjusts=c(-0.6, 1.5, -0.4, -0.4, -0.4, 1.5, -0.4, -0.5, -0.5, -0.6, -0.6, 1.5, -0.6, 1.5, 1.5)
  results$plot_motif_names = plot_motif_names

  # generate ari/comp plot
  pl = ggplot(data=results, mapping=aes(x=CompSize, y=ARI)) +
    geom_abline(slope=0, intercept=0, linetype='dashed', color='gray50') +
    geom_abline(slope=0, intercept=1, linetype='dashed', color='gray50') +
    geom_point() +
    xlab(as.expression(bquote('Connected component size, |'*italic(C)*'|'))) +
    ylab('ARI') +
    coord_cartesian(xlim=c(0,1222), ylim = c(-0.5,1), clip='off') +
    scale_x_continuous(breaks=c(0,500,1000,1222)) +
    geom_text(data=results, mapping=aes(x=CompSize, y=ARI, label=plot_motif_names), hjust=hjusts, vjust=vjusts, parse=TRUE, family='CM Roman',size=5) +
    (theme_diss() + theme(plot.margin = margin(t=1,r=1)) + theme(axis.title.y = element_text(margin=margin(r=-12))))

  # print ari/comp plot to pdf
  pdf_name = '../results/polblogs/polblogs_ari_conn.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)




  # run M12 again
  cat('\tGenerating eigenvector plots\t\t\n')
  motif_name = 'M12'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 3)
  comps = motifEmb$comps
  vects = motifEmb$vects
  clusts = kmeanspp(vects[,2],k=2)$cluster

  if(clusts[1] == 2){
    clusts = 3 - clusts
  }

  M12_df = as.data.frame(vects)
  M12_df$truth = as.character(truth[comps])
  M12_df$clusts = as.character(clusts)

  # plot eigs for M12 with clusts coloring
  pl = ggplot(data=M12_df, mapping=aes(x=V2, y=V3, color=clusts)) +
    geom_point(alpha=0.4) +
    xlab('Eigenvector 2') +
    ylab('Eigenvector 3') +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    (theme_diss() + theme(plot.margin = margin(b=2,r=1,t=1), legend.position = c(0.135,0.08),
                          legend.title=element_blank(), legend.text=element_text(size=20), axis.title.y=element_text(margin=margin(r=10)))) +
    scale_color_manual(name='Cluster',labels=c('Cluster 1','Cluster 2'),values=c('dodgerblue','firebrick2'))

  # print M12 clusts plot to pdf
  pdf_name = '../results/polblogs/polblogs_M12_clusts.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)

  # plot eigs for M12 with truth coloring
  pl = ggplot(data=M12_df, mapping=aes(x=V2, y=V3, color=truth)) +
    geom_point(alpha = 0.4) +
    xlab('Eigenvector 2') +
    ylab('Eigenvector 3') +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    (theme_diss() + theme(plot.margin = margin(b=2,r=1,t=1), legend.position = c(0.175,0.08),
                          legend.title=element_blank(), legend.text=element_text(size=20), axis.title.y=element_text(margin=margin(r=10)))) +
    scale_color_manual(name='Truth',labels=c('Liberal','Conservative'),values=c('dodgerblue','firebrick2'))

  # print M12 truth plot to pdf
  pdf_name = '../results/polblogs/polblogs_M12_truth.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)

  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')

  return()
}

experimentMigration = function(motif_names, eigenvectors_for_map, n_vects_for_cluster, n_clusters, simplify_sweep){

  cat('US Migration network\n')

  t0 = Sys.time()

  # read data
  cat('\tReading data...\n')
  Ggraph = readGraphData('us_migration')
  counties = Ggraph$counties
  states = Ggraph$states
  G = Ggraph$adj

  cat('\tNumber of counties:', nrow(G), '\n')
  cat('\tNumber of states:', nrow(states), '\n')
  cat('\tNumber of edges:', sum(G>0), '\n')

  # capped at 1e4
  G = G*(G <= 1e4) + 1e4*(G > 1e4)
  G = drop0(G)


  # plot map with states labelled
  cat('\tPlotting reference map\n')
  pl = ggplot() +
    geom_sf(data=states, color='black', size=0.2, aes(geometry=geometry, fill=shading)) +
    scale_fill_manual(values=c('white','gray92','gray84','gray75')) +
    coord_sf(crs = st_crs(2163)) +
    geom_text(data=states, aes(x=X, y=Y, label=abbr_name), family='CM Roman', size=states$label_sizes) +
    geom_segment(data=states, aes(x=label_line_x, xend=label_line_xend, y=label_line_y, yend=label_line_yend), size=0.15, color='black') +
    theme(panel.background = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = margin(t=0, r=0, b=0, l=0),
          legend.position='none'
    )

  # print pdf
  pdf_name = '../results/us_migration/us_migration_map_state_names.pdf'
  suppressWarnings(ggsave(filename=pdf_name, plot=pl, width=6, height=4))
  embed_fonts(file=pdf_name, outfile=pdf_name)



  # get baseline embeddings for calibration
  cat('\tGetting baseline embedding\n')
  k = max(eigenvectors_for_map)
  motifEmbBase = runMotifEmb(G, 'M6', 'func', 'rw', k)
  vectsBase = motifEmbBase$vects
  set.seed(314159)
  clustsBase = kmeanspp(vectsBase[,2:n_vects_for_cluster], n_clusters)$cluster

  for(motif_name in motif_names){

    # run algorithm
    cat('\tMotif ',motif_name,' embedding:\n',sep='')
    motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', k)
    vects = motifEmb$vects
    clusts = kmeanspp(vects[,2:n_vects_for_cluster], n_clusters)$cluster
    clusts = makeClusterLike(clusts, clustsBase)
    counties$clust = as.character(clusts)
    n_cut_value = scorePartition(G, clusts, 'n_cut', directed=TRUE)
    n_cut = data.frame(x=-13e5, y=-205e4, label=paste('Ncut =', format(round(n_cut_value,2), nsmall=2)))




    # plot sweep profiles
    sweep_scores = sweepScores(G, vects[,2], 'n_cut', simplify_sweep)
    sweep_profile = data.frame(x=1:length(sweep_scores), y=sweep_scores)
    sweepClusts = minSweepPartition2(G, vects[,2], sweep_scores, 'n_cut')
    pl = ggplot(data=sweep_profile, aes(x=x,y=y)) +
      geom_line() +
      ggtitle('') +
      xlab(as.expression(bquote('Splitting point'))) +
      ylab('Ncut score') +
      scale_y_continuous(limits=c(0,0.5)) +
      (theme_diss() + theme(plot.margin = margin(t=-5,r=10,b=5)))

    file_name = paste('../results/us_migration/sweepClusts_', motif_name, '.txt', sep='')
    write(sweepClusts, file_name, ncolumns=1)

    pdf_name = paste('../results/us_migration/us_migration_sweep_profile_',motif_name,'.pdf',sep='')
    suppressWarnings(ggsave(filename=pdf_name, plot=pl, width=4, height=4))
    embed_fonts(pdf_name)


    # print vect-coloured maps
    for(i in eigenvectors_for_map){
      cat('\t\tPlotting eigenvector',i,'\n')

      # flip sign if necessary so maps look similar
      if(sum(abs(-vects[,i] - vectsBase[,i])) < sum(abs(vects[,i] - vectsBase[,i]))){
        counties$vect = -vects[,i]
      }
      else{
        counties$vect = vects[,i]
      }

      # make plot
      pl = ggplot() +
        geom_sf(data=counties, size=0.05, aes(geometry=geometry, fill=vect, color=vect)) +
        scale_color_gradientn(colors = cet_pal(100, name = "d13")) +
        scale_fill_gradientn(colors = cet_pal(100, name = "d13")) +
        geom_sf(data=states, color='black', fill=NA, size=0.2, aes(geometry=geometry)) +
        coord_sf(crs = st_crs(2163)) +
        theme(panel.background = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = margin(t=0, r=0, b=0, l=0),
          legend.position='none'
        )

      # print pdf
      pdf_name = paste('../results/us_migration/us_migration_',motif_name,'_',i,'.pdf',sep='')
      suppressWarnings(ggsave(filename=pdf_name, plot=pl, width=6, height=4))
    }



    # print clust-coloured maps
    cat('\t\tPlotting clusters\n')
    # make plot
    pl = ggplot() +
      geom_sf(data=counties, size=0.05, aes(geometry=geometry, fill=clust, color=clust)) +
      scale_color_brewer(palette='Set3') +
      scale_fill_brewer(palette='Set3') +
      geom_sf(data=states, color='black', fill=NA, size=0.2, aes(geometry=geometry)) +
      geom_text(data=n_cut, aes(x=x, y=y, label=label), family='CM Roman', size=10) +
      coord_sf(crs = st_crs(2163)) +
      theme(panel.background = element_blank(),
            axis.line=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin = margin(t=0, r=0, b=0, l=0),
            legend.position='none'
      )

    # print pdf
    pdf_name = paste('../results/us_migration/us_migration_clusts_',motif_name,'.pdf',sep='')
    suppressWarnings(ggsave(filename=pdf_name, plot=pl, width=6, height=4))
    embed_fonts(pdf_name, outfile=pdf_name)

  }

  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')
  return()
}

us_migration_table = function(motif_names, eigenvectors_for_map){

  string = '\\hline & $\\ca{M}_s$ & $\\ca{M}_6$ & $\\ca{M}_9$ \\\\ '

  for(i in eigenvectors_for_map){

    string = paste(string, '\\hline $ \\ x_{', i, '} \\ $ & ', sep='')

    for(motif_name in motif_names){
      string = paste(string, '\\begin{minipage}{0.31\\textwidth} \\centering \\includegraphics[scale=0.3,draft=false]{../../results/us_migration/us_migration_',motif_name,'_',i,'.pdf} \\end{minipage} &', sep='')
    }

    string = paste(substr(string,1,nchar(string)-1), '\\\\')
  }

  string = paste(string, '\\hline \\ $C$ \\ & ')

  for(motif_name in motif_names){
    string = paste(string, '\\begin{minipage}{0.31\\textwidth} \\centering \\includegraphics[scale=0.3,draft=false]{../../results/us_migration/us_migration_clusts_',motif_name,'.pdf} \\end{minipage} &', sep='')
  }

  string = paste(substr(string,1,nchar(string)-1), '\\\\ \\hline')

  file_name = '../results/us_migration/us_migration_table.txt'
  write(string, file_name)
}

us_migration_sweep_ari = function(motif_names){

  for(i in motif_names){
    for(j in motif_names){

      file_name_i = paste('../results/us_migration/sweepClusts_', i, '.txt', sep='')
      clusts_i = scan(file_name_i, quiet=TRUE)

      file_name_j = paste('../results/us_migration/sweepClusts_', j, '.txt', sep='')
      clusts_j = scan(file_name_j, quiet=TRUE)

      randIndex = adjustedRandIndex(clusts_i, clusts_j)

      cat(i, j, randIndex, '\n')
    }
  }

}

experimentBipartite = function(ns, p1s, p2s, n_repeats){

  cat('BSBM\n')

  t0 = Sys.time()

  n_experiments = length(ns)
  motif_names = c('coll','expa','SVD_source','SVD_dest')
  plot_motif_names = c('M[coll]','M[expa]','Co*"-"*clust[italic(S)]','Co*"-"*clust[italic(D)]')

  for(expmt in 1:n_experiments){

    # prepare experiment parameters and results tables
    n = ns[expmt]
    p1 = p1s[expmt]
    p2 = p2s[expmt]
    dest_block_sizes = rep(n,2)
    targ_block_sizes = rep(n,2)
    truth = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n))

    # prepare results table and row counter
    resultsColNames = c('Motif','ARI','CompSize')
    resultsNRows = n_repeats * length(motif_names)
    results = dataFrameNames(resultsColNames, resultsNRows)
    row_number = 1

    for(rep in 1:n_repeats){

      # sample the graph
      G = sampleBSBM(dest_block_sizes, targ_block_sizes, p1, p2)

      for(motif_name in motif_names[1:2]){

        # run algorithm
        motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 3)
        comps = motifEmb$comps
        vects = motifEmb$vects
        clusts = kmeanspp(vects[,2],k=2)$cluster
        randIndex = adjustedRandIndex(clusts,truth[comps])

        # enter results and increment counter
        results[row_number,1] = motif_name
        results[row_number,2] = randIndex
        results[row_number,3] = length(comps)
        row_number = row_number + 1

        # progress
        cat('\tExperiment ',expmt,'/',n_experiments,', repeat ',rep,'/',n_repeats,', motif ',motif_name,'\t\t\t\r',sep='')
      }

      # run SVD embedding
      vect = runSVDEmb(G)
      clustsSource = kmeanspp(vect[1:(2*n)],k=2)$cluster
      clustsDest = kmeanspp(vect[(2*n+1):(4*n)],k=2)$cluster
      randIndexSource = adjustedRandIndex(clustsSource,truth[1:(2*n)])
      randIndexDest = adjustedRandIndex(clustsDest,truth[(2*n+1):(4*n)])

      # enter results and increment counter
      results[row_number,1] = motif_names[3]
      results[row_number,2] = randIndexSource
      results[row_number,3] = 2*n
      row_number = row_number + 1
      results[row_number,1] = motif_names[4]
      results[row_number,2] = randIndexDest
      results[row_number,3] = 2*n
      row_number = row_number + 1
    }

    # generate mean connected component sizes
    meanCompSize = data.frame(meanCompSize = round(c(by(results$CompSize, results$Motif, mean))[rank(motif_names)]))

    # add x-coordinates of violins
    desired_x_coords = c(1,2,3.2,4.2)
    midpoint_x_coords = (desired_x_coords[2] + desired_x_coords[3]) / 2
    results$x_coord_violin = NA
    results$x_coord_violin[results[,1] == motif_names[1]] = desired_x_coords[1]
    results$x_coord_violin[results[,1] == motif_names[2]] = desired_x_coords[3]
    results$x_coord_violin[results[,1] == motif_names[3]] = desired_x_coords[2]
    results$x_coord_violin[results[,1] == motif_names[4]] = desired_x_coords[4]

    # make dividing line
    dividing_line = data.frame(x=midpoint_x_coords, xend=midpoint_x_coords, y=-1, yend=2)

    # set box widths
    results$width = 5



    # generate plot
    pl = ggplot() +
      geom_abline(slope=0, intercept=0, linetype='dashed', color='gray50') +
      geom_abline(slope=0, intercept=1, linetype='dashed', color='gray50') +
      geom_violin(data=results, mapping=aes(x=x_coord_violin, y=ARI, fill=factor(x_coord_violin)), scale='width', color='black',width=0.7) +
      scale_fill_manual(values=rep('gray90',4)) +
      scale_x_continuous(breaks=desired_x_coords, labels=parse(text=plot_motif_names)[c(1,3,2,4)]) +
      stat_summary(data=results, mapping=aes(x=x_coord_violin, y=ARI), fun.y=median, geom='point', size=1.5) +
      stat_summary(data=results, mapping=aes(x=x_coord_violin, y=ARI), fun.y=function(z){10}, fun.ymin=function(z){quantile(z,0.25)}, fun.ymax=function(z){quantile(z,0.75)}, size=0.5) +
      xlab('Method') +
      ylab('ARI') +
      geom_segment(data=dividing_line, aes(x=x, xend=xend, y=y, yend=yend)) +
      coord_cartesian(ylim = c(-0.5,1), clip='on') +

      (theme_diss() +
         theme(plot.title = element_text(vjust=7)) +
         theme(plot.margin = margin(t=20, r=10, b=0, l=0)) +
         theme(legend.position='none') +
         theme(axis.title.y = element_text(margin=margin(r=-12))))

    # print to pdf
    pdf_name = paste('../results/bipartite/bipartite',expmt,'.pdf',sep='')
    ggsave(filename=pdf_name, plot=pl, width=7, height=5)
    embed_fonts(pdf_name, outfile=pdf_name)
  }

  cat('\n')
  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')
  return()
}

experimentAmericanRevolution = function(){

  cat('American Revolution network\n')

  t0 = Sys.time()

  # read data
  graph_data = readGraphData('american_revolution')
  G = graph_data$adj

  cat('\tNumber of source vertices (people):', sum(apply(G,1,sum)>0), '\n')
  cat('\tNumber of destination vertices (organisations):', sum(apply(G,2,sum)>0), '\n')
  cat('\tTotal number of vertices:', nrow(G), '\n')
  cat('\tNumber of edges:', sum(G>0), '\n')


  # cluster source with coll
  motif_name = 'coll'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 5)
  comps = motifEmb$comps
  vects = motifEmb$vects
  set.seed(2642)
  clusts = kmeanspp(vects[,2:5],k=5)$cluster

  # plot source clustering
  G_graph = graph_from_adjacency_matrix(G)
  set.seed(6548)
  plotcoord = data.frame(layout.fruchterman.reingold(G_graph)) # TODO add a seed to this layout?
  plotcoord$cols = c(as.character(clusts),rep('0',5))
  colnames(plotcoord) = c("X1","X2",'cols')
  edgelist = get.edgelist(G_graph)
  edges = data.frame(plotcoord[edgelist[,1],1:2], plotcoord[edgelist[,2],1:2])
  colnames(edges) = c("X1","Y1","X2","Y2")

  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='gray', size=0.5) +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord[1:136,], size=2) +
    scale_color_brewer(palette='Set1', labels=paste('Cluster',1:5)) +
    geom_point(aes(x=X1,y=X2), data=plotcoord[137:141,], color='gray50', size=3, shape=15) +
    xlab('') +
    ylab('') +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                          legend.position = c(0.12,0.83), legend.title=element_blank(), legend.text=element_text(size=20),
                          plot.margin=margin(r=5,l=-20,t=5,b=-20), panel.border=element_blank(), panel.grid=element_blank()))


  # print source clustering pdf
  pdf_name = '../results/american_revolution/american_revolution_source.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)





  # cluster dest with expa
  motif_name = 'expa'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 2)
  comps = motifEmb$comps
  vects = motifEmb$vects
  set.seed(9834)
  clusts = kmeanspp(vects[,2],k=2)$cluster

  # plot dest clustering
  plotcoord$cols = c(rep('0',136),as.character(clusts))

  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='gray', size=0.5) +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord[137:141,], size=3, shape=15) +
    scale_color_brewer(palette='Set1', labels=paste('Cluster',1:2)) +
    geom_point(aes(x=X1,y=X2), data=plotcoord[1:136,], color='gray50', size=2) +
    xlab('') +
    ylab('') +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                          legend.position = c(0.12,0.9), legend.title=element_blank(), legend.text=element_text(size=20),
                          plot.margin=margin(r=5,l=-20,t=5,b=-20), panel.border=element_blank(), panel.grid=element_blank()))


  # print dest clustering pdf
  pdf_name = '../results/american_revolution/american_revolution_dest.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)

  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')
  return()
}

experimentLanguages = function(eigenvectors_for_map, eigenvectors_for_cluster, n_clusters_to_print){

  cat('Unicode Languages network\n')
  cat('\tReading data...\n')

  t0 = Sys.time()

  # number of eigenvectors to plot
  k = max(eigenvectors_for_map)

  # read data
  graph_data = readGraphData('languages')
  G = graph_data$adj
  countries = graph_data$countries
  languages = graph_data$languages
  populations = graph_data$populations
  country_names = graph_data$country_names
  language_names = graph_data$language_names
  language_speakers = graph_data$language_speakers
  n_source = length(countries)
  n_dest = length(languages)


  # cluster source with coll
  motif_name = 'coll'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', k)
  comps = motifEmb$comps
  vects = motifEmb$vects
  allVects = matrix(NA, nrow=n_source, ncol=ncol(vects))
  allVects[comps,] = vects

  # print network summaries
  cat('\tNumber of source vertices (territories):', length(comps), '\n')
  cat('\tNumber of destination vertices (languages):', nrow(G)-length(comps), '\n')
  cat('\tTotal number of vertices:', nrow(G), '\n')
  cat('\tNumber of edges:', sum(G[largestComponent(G),largestComponent(G)]>0), '\n')

  # Kazakhs who speak German
  kazakhstan_index = which(countries=='KZ')
  german_index = which(languages=='de')
  kazakhs_german = G[kazakhstan_index, (n_source+german_index)]
  cat('\tKazakhs who speak German:', kazakhs_german, '\n')

  # Panamanians who speak Chinese
  panama_index = which(countries=='PA')
  chinese_index = which(languages=='zh')
  panamanians_chinese = G[panama_index, (n_source+chinese_index)]
  cat('\tPanamanians who speak Chinese:', panamanians_chinese, '\n')

  # Japanese who speak Korean
  japan_index = which(countries=='JP')
  korean_index = which(languages=='ko')
  japanese_korean = G[japan_index, (n_source+korean_index)]
  cat('\tJapanese who speak Korean:', japanese_korean, '\n')


  # read world map
  world_map = ne_countries(scale = 'small', returnclass = "sf")
  world_map$vects = rep(NA, nrow(world_map))
  world_map$clusts = rep(NA, nrow(world_map))

  world_map$iso_a2[world_map$name == 'Norway'] = 'NO'
  world_map$iso_a2[world_map$name == 'France'] = 'FR'
  world_map$iso_a2[world_map$name == 'Somaliland'] = 'SO'
  world_map$iso_a2[world_map$name == 'N. Cyprus'] = 'CY'
  world_map$iso_a2[world_map$name == 'Baikonur'] = 'KZ'
  world_map = world_map[world_map$iso_a2 != 'TF',]


  # trim map
  world_map = suppressMessages(suppressWarnings(st_crop(world_map, c(xmin=-168.6, xmax=180, ymin=-60, ymax=90))))


  # generate clusters
  set.seed(3596)
  clusts = kmeanspp(vects[,eigenvectors_for_cluster],k=n_clusters_to_print)$cluster
  ordered_clusts = rep(NA, length(clusts))
  for(i in 1:n_clusters_to_print){
    ordered_clusts[clusts == order(table(clusts), decreasing=TRUE)[i]] = i
  }
  allClusts = rep(0, n_source)
  allClusts[comps] = ordered_clusts





  # plot world map with clust colouring
  cat('\tMaking map\n')

  # add clusts values to map
  for(j in 1:n_source){
    if(countries[j] %in% world_map$iso_a2){
      ind = which(countries[j] == world_map$iso_a2)
      world_map$clusts[ind] = as.character(allClusts[j])
    }
  }

  # plot map
  na_color = '#DDDDDD'
  pl = ggplot(data=world_map, aes(geometry=geometry)) +
    geom_sf(aes(fill=clusts), size=0.05, color='black') +
    scale_fill_brewer(palette='Set3', na.value=na_color, labels=c(paste('   Cluster',1:6),'   No cluster')) +
    coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs") +
    (theme_diss() +
    theme(panel.background = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = margin(t=-80, r=100, b=-80, l=-50),
          legend.position=c(1.05,0.65),
          legend.key.size=unit(0.7, "cm"),
          legend.title=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(colour='transparent'),
          legend.text = element_text(size=15)
    ))


  # save map
  pdf_name = '../results/languages/languages_source_map_clusts.pdf'
  ggsave(filename=pdf_name, plot=pl, width=9, height=3.8)
  embed_fonts(pdf_name)

  # record colours used
  clust_cols = c(brewer.pal(6,'Set3'), na_color)
  clust_cols = substr(clust_cols, 2, nchar(clust_cols))










  # format clusters for tex input
  max_length = max(table(allClusts))
  padded_country_names = matrix('', nrow=max_length, ncol=n_clusters_to_print)

  for(i in 1:n_clusters_to_print){
    inds = (allClusts==i)
    padded_country_names[(1:sum(inds)),i] = (country_names[inds])[order(populations[inds], decreasing=TRUE)]
  }

  # only print a few in each cluster
  padded_country_names = padded_country_names[1:20,]

  # add dots and cluster size
  padded_country_names = rbind(padded_country_names, '')
  padded_country_names = rbind(padded_country_names, '')
  for(i in 1:n_clusters_to_print){
    if(padded_country_names[nrow(padded_country_names)-2,i] != ''){
      padded_country_names[nrow(padded_country_names)-1,i] = '$\\cdots$'
    }
    padded_country_names[nrow(padded_country_names),i] = paste('$|\\textrm{Cluster\\ }',i,'|$ =', sum(allClusts == i))
  }


  # write table to data string
  data_string = ''

  # cluster headers
  for(i in 1:n_clusters_to_print){
    data_string = paste(data_string, ' \\cellcolor[HTML]{',clust_cols[i],'} Cluster ', i, ' & ', sep='')
  }
  data_string = paste(substr(data_string, 1, nchar(data_string)-2), '\\\\[0.1cm] \\hline \\rule{0pt}{1.2em}')


  for(row_number in 1:nrow(padded_country_names)){
    for(i in 1:n_clusters_to_print){
      data_string = paste(data_string, padded_country_names[row_number,i], '&')
    }

    data_string = substr(data_string, 1, nchar(data_string)-1)
    data_string = paste(data_string, '\\\\')
  }
  data_string = paste(substr(data_string, 1, nchar(data_string)-3), '')

  filename = '../results/languages/languages_source_clusters.txt'
  write(data_string, file=filename)







  # cluster dest with expa
  motif_name = 'expa'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', k)
  comps = motifEmb$comps - n_source
  vects = motifEmb$vects
  set.seed(6320)
  clusts = kmeanspp(vects[,eigenvectors_for_cluster],k=n_clusters_to_print)$cluster
  ordered_clusts = rep(NA, length(clusts))
  for(i in 1:n_clusters_to_print){
    ordered_clusts[clusts == order(table(clusts), decreasing=TRUE)[i]] = i
  }
  allClusts = rep(0, n_dest)
  allClusts[comps] = ordered_clusts


  # format clusters for tex input
  max_length = max(table(allClusts))
  padded_language_names = matrix('', nrow=max_length, ncol=n_clusters_to_print)

  for(i in 1:n_clusters_to_print){
    inds = (allClusts==i)
    padded_language_names[(1:sum(inds)),i] = (language_names[inds])[order(language_speakers[inds], decreasing=TRUE)]
  }

  # only print a few in each cluster
  padded_language_names = padded_language_names[1:20,]

  # add dots and cluster size
  padded_language_names = rbind(padded_language_names, '')
  padded_language_names = rbind(padded_language_names, '')
  for(i in 1:n_clusters_to_print){
    if(padded_language_names[nrow(padded_language_names)-2,i] != ''){
      padded_language_names[nrow(padded_language_names)-1,i] = '$\\cdots$'
    }
    padded_language_names[nrow(padded_language_names),i] = paste('$|\\textrm{Cluster\\ }',i,'|$ =', sum(allClusts == i))
  }


  # write table to data string
  data_string = ''

  # headers
  for(i in 1:n_clusters_to_print){
    data_string = paste(data_string, 'Cluster', i, '&')
  }
  data_string = paste(substr(data_string, 1, nchar(data_string)-2), '\\\\[0.1cm] \\hline \\rule{0pt}{1.2em}')


  for(row_number in 1:nrow(padded_language_names)){
    for(i in 1:n_clusters_to_print){
      data_string = paste(data_string, padded_language_names[row_number,i], '&')
    }

    data_string = substr(data_string, 1, nchar(data_string)-1)
    data_string = paste(data_string, '\\\\')
  }
  data_string = paste(substr(data_string, 1, nchar(data_string)-3), '')

  filename = '../results/languages/languages_dest_clusters.txt'
  write(data_string, file=filename)


  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')
  return()
}

experimentDirectedCP = function(ns, p1s, p2s, n_repeats){

  t0 = Sys.time()

  n_experiments = length(ns)
  motif_names = c('Ms','M1','M5','M8','M9','M10')
  plot_motif_names = c('M[italic(s)]','M[1]','M[5]','M[8]','M[9]','M[10]')

  for(expmt in 1:n_experiments){

    # prepare experiment parameters and results tables
    n = ns[expmt]
    p1 = p1s[expmt]
    p2 = p2s[expmt]
    block_sizes = rep(n,4)
    truth = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n))

    # prepare results table and row counter
    resultsColNames = c('Motif','ARI','CompSize')
    resultsNRows = n_repeats * length(motif_names)
    results = dataFrameNames(resultsColNames, resultsNRows)
    row_number = 1

    for(rep in 1:n_repeats){

      # sample the graph
      G = sampleDirectedCP(block_sizes, p1, p2)

      for(motif_name in motif_names){

        # run algorithm
        motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 2)
        comps = motifEmb$comps
        vects = motifEmb$vects
        clusts = kmeanspp(vects[,2],k=4)$cluster
        randIndex = adjustedRandIndex(clusts,truth[comps])

        # enter results and increment counter
        results[row_number,1] = motif_name
        results[row_number,2] = randIndex
        results[row_number,3] = length(comps)
        row_number = row_number + 1

        # progress
        motif_name_print = substr(paste(motif_name,',    ',sep=''),1,5)
        randIndex_print = substr(paste(round(randIndex,2),'000000',sep=''),1,4)
        cat('experiment ',expmt,'/',n_experiments,', repeat ',rep,'/',n_repeats,', motif ',motif_name_print,'ARI = ',randIndex_print,'\t\t\t\n',sep='')
      }
    }

    # generate mean connected component sizes
    meanCompSize = data.frame(meanCompSize = round(c(by(results$CompSize, results$Motif, mean))[rank(motif_names)], 1))

    # generate plot
    pl = ggplot() +
      geom_boxplot(data=results, mapping=aes(x=Motif, y=ARI), outlier.size=0.5, color='black') +
      scale_x_discrete(limits=motif_names, labels=parse(text=plot_motif_names)) +
      xlab('Motif') +
      ylab('ARI') +
      coord_cartesian(ylim = c(-0.5,1), clip='off') +
      geom_text(data=meanCompSize, mapping=aes(x=1:6, y=1.15, label=meanCompSize), family='CM Roman') +
      ggtitle(as.expression(bquote(italic(n)*" = "*.(n)*", "*italic(p)*" = "*.(p1)*", "*italic(q)*" = "*.(p2)))) +
      (theme_diss() +
         theme(plot.title = element_text(vjust=7)) +
         theme(plot.margin = margin(t=20, r=10, b=0, l=0))) +
      geom_text(mapping=aes(x=0.2, y=1.14, label='"|"*italic(C)*"| = "'), parse=TRUE, family='CM Roman')

    # print to pdf
    pdf_name = paste('../results/directed_cp/directed_cp',expmt,'.pdf',sep='')
    ggsave(filename=pdf_name, plot=pl, width=7, height=5)
    embed_fonts(pdf_name, outfile=pdf_name)
  }

  print(Sys.time() - t0)
  return()
}

experimentFacultyHiring = function(){

  faculty_ids = c('business','cs','history')

  # loop over business, cs, history
  for(faculty_num in 1:1){

    # read data
    faculty_id = faculty_ids[faculty_num]
    graph_data = readGraphData(paste('faculty_', faculty_id, sep=''))
    G = graph_data$adj
    unis = graph_data$unis

    # run algorithm
    motifEmb = runMotifEmb(G, 'M5', type='func', 'rw', 3)
    vects = motifEmb$vects
    comps = motifEmb$comps
    G = G[comps,comps]
    unis = unis[comps,]

    clusts = kmeanspp(vects[,2],k=4)$cluster
    sortClusts = clustDirectedCP(G, clusts)

    # plot
    df = as.data.frame(vects)
    df$sortClusts = as.character(sortClusts)
    pl = ggplot(data=df) +
      geom_point(aes(x=V2,y=V3,col=sortClusts), alpha=0.6) +
      scale_color_jcolors(name='Cluster',labels=c('P-out','C-in','C-out','P-in'),palette = 'pal3') +
      xlab('Eigenvector 2') +
      ylab('Eigenvector 3') +
      ggtitle('M5, k-means++') +
      (theme_diss() +
         theme(plot.margin = margin(b=2)))

    # print pdf
    pdf_name = paste('../results/faculty_hiring/faculty_hiring_', faculty_id, '.pdf', sep='')
    ggsave(filename=pdf_name, plot=pl, width=7, height=5)
    embed_fonts(pdf_name, outfile=pdf_name)

    # format and print clusters for tex input
    max_length = max(sum(sortClusts==1), sum(sortClusts==2), sum(sortClusts==3), sum(sortClusts==4))
    padded_unis = matrix(' ', nrow=max_length, ncol=4)

    for(i in 1:4){
      padded_unis[(1:sum(sortClusts==i)),i] = sort(as.character(unis$institution[sortClusts==i]))
    }

    data_string = ''

    for(row_number in 1:max_length){
      for(i in 1:4){
        data_string = paste(data_string, padded_unis[row_number,i], '&')
      }

      data_string = substr(data_string, 1, nchar(data_string)-1)
      data_string = paste(data_string, '\\\\')
    }

    data_string = gsub('Texas A&M', 'Texas A and M', data_string)
    data_string = gsub('University','Uni',data_string)
    data_string = gsub('Thunderbird School of Global Management','Thunderbird School',data_string)
    data_string = gsub(', Binghamton','',data_string)
    data_string = gsub('North Carolina','N. Carolina',data_string)
    data_string = gsub('Uni of Illinois, Urbana Champaign','Uni of Illinois, UC',data_string)
    data_string = gsub('Technology','Tech',data_string)
    data_string = gsub('Polytechnic','Poly',data_string)
    data_string = substr(data_string, 1, nchar(data_string)-2)

    filename = paste('../results/faculty_hiring/', faculty_id,'.txt', sep='')
    write(data_string, file=filename)

    # plot network
    graph_plot = graph_from_adjacency_matrix(G, weighted=TRUE)
    pdf_name = '../results/faculty_hiring/faculty_network.pdf'
    pdf(pdf_name)
    plot(graph_plot, vertex.size=3*apply(G,1,sum)^0.2, edge.width=0.1*edge_attr(graph_plot,'weight'), vertex.label=NA, edge.arrow.size=0.5, edge.arrow.width=0.05, vertex.color=c('pink','blue','red','cyan')[sortClusts])
    dev.off()
  }

  return()
}

experimentPollination = function(){

  # read data
  graph_data = readGraphData('clements_pollination')
  G = graph_data$adj
  pollinators = graph_data$pollinators
  plants = graph_data$plants
  n_source = length(pollinators)
  n_dest = length(plants)

  # cluster source with coll
  k = 4

  motif_name = 'coll'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', k)
  comps = motifEmb$comps
  vects = motifEmb$vects
  clusts = kmeanspp(vects[,2:k],k=k)$cluster
  allClusts = rep(0, n_source)
  allClusts[comps] = clusts

  # plot source clustering
  G_graph = graph_from_adjacency_matrix(G, weighted=TRUE)
  plotcoord = data.frame(layout.fruchterman.reingold(G_graph)) # TODO seed?
  plotcoord$cols = c(as.character(allClusts),rep(as.character(k+1), n_dest)) # 0: unclustered, k+1: dest
  colnames(plotcoord) = c("X1","X2",'cols')
  edgelist = get.edgelist(G_graph)
  edges = data.frame(plotcoord[edgelist[,1],1:2], plotcoord[edgelist[,2],1:2])
  colnames(edges) = c("X1","Y1","X2","Y2")

  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='gray', size=0.3) +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord[1:n_source,], size=1) +
    scale_color_jcolors(palette='pal3', name='Clusters', labels=c('No Cluster',paste(rep('Cluster',k),1:k))) +
    geom_point(aes(x=X1,y=X2), data=plotcoord[(n_source+1):(n_source+n_dest),], color='black', size=1) +
    ggtitle(paste('Clustering pollinators into',k,'clusters')) +
    xlab('') +
    ylab('') +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank()))

  # print source clustering pdf
  pdf_name = '../results/clements_pollination/clements_pollination_source.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=4)
  embed_fonts(pdf_name, outfile=pdf_name)

  # print clusters
  for(i in 1:k){
    filename = paste('../results/clements_pollination/pollinators_cluster_',i,'.txt', sep='')
    datatable = pollinators[allClusts==i]
    datatable = gsub(pattern='&', replacement=' and ', datatable)
    write.table(datatable, file=filename, row.names=FALSE, col.names=FALSE, quote=FALSE, eol='\\\\')
  }



  # cluster dest with expa
  k = 4

  motif_name = 'expa'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', k)
  comps = motifEmb$comps - n_source
  vects = motifEmb$vects
  clusts = kmeanspp(vects[,2:k],k=k)$cluster
  allClusts = rep(0, n_dest)
  allClusts[comps] = clusts

  # plot dest clustering
  plotcoord$cols = c(rep(as.character(k+1), nrow(plotcoord)-length(allClusts)),as.character(allClusts))

  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='gray', size=0.3) +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord[(n_source+1):(n_source+n_dest),], size=1) +
    scale_color_jcolors(palette='pal3', name='Clusters', labels=c('No Cluster',paste(rep('Cluster',k),1:k))) +
    geom_point(aes(x=X1,y=X2), data=plotcoord[1:n_source,], color='black', size=1) +
    ggtitle(paste('Clustering plants into',k,'clusters')) +
    xlab('') +
    ylab('') +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank()))

  # print source clustering pdf
  pdf_name = '../results/clements_pollination/clements_pollination_dest.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=4)
  embed_fonts(pdf_name, outfile=pdf_name)

  # print clusters
  for(i in 1:k){
    filename = paste('../results/clements_pollination/plants_cluster_',i,'.txt', sep='')
    datatable = plants[allClusts==i]
    datatable = gsub(pattern='&', replacement=' and ', datatable)
    write.table(datatable, file=filename, row.names=FALSE, col.names=FALSE, quote=FALSE, eol='\\\\')
  }

  return()
}

experimentUCIrvine = function(){

  # TODO automate magic numbers in bipartite sizes

  source_size = 899
  dest_size = 522

  # read data
  graph_data = readGraphData('uc_irvine_forum')
  G = graph_data$adj

  # cluster source with coll
  motif_name = 'coll'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 5)
  comps = motifEmb$comps
  vects = motifEmb$vects
  clusts = kmeanspp(vects[,2:5],k=5)$cluster

  print(length(comps))

  # plot source clustering
  G_graph = graph_from_adjacency_matrix(G)
  plotcoord = data.frame(layout.fruchterman.reingold(G_graph)) # TODO add a seed to this layout?
  print('here')

  plotcoord$cols = c(as.character(clusts),rep('0', dest_size))

  colnames(plotcoord) = c("X1","X2",'cols')
  edgelist = get.edgelist(G_graph)
  edges = data.frame(plotcoord[edgelist[,1],1:2], plotcoord[edgelist[,2],1:2])
  colnames(edges) = c("X1","Y1","X2","Y2")


  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='gray', size=0.5) +
    #geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord[1:136,], size=2) +
    scale_color_jcolors(palette='pal3', name='Clusters', labels=c('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5')) +
    #geom_point(aes(x=X1,y=X2), data=plotcoord[137:141,], color='black', size=2) +
    ggtitle('Clustering people into 5 clusters') +
    xlab('') +
    ylab('') +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank()))

  # print source clustering pdf
  pdf_name = '../results/uc_irvine_forum/uc_irvine_forum_source.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=4)
  embed_fonts(pdf_name, outfile=pdf_name)



  # cluster dest with expa
  motif_name = 'expa'
  motifEmb = runMotifEmb(G, motif_name, 'func', 'rw', 2)
  comps = motifEmb$comps
  vects = motifEmb$vects
  clusts = kmeanspp(vects[,2],k=2)$cluster

  # plot dest clustering
  plotcoord$cols = c(rep('0',136),as.character(clusts))

  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='gray', size=0.5) +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord[137:141,], size=2) +
    scale_color_jcolors(palette='pal3', name='Clusters', labels=c('Cluster 1','Cluster 2')) +
    geom_point(aes(x=X1,y=X2), data=plotcoord[1:136,], color='black', size=2) +
    ggtitle('Clustering organisations into 2 clusters') +
    xlab('') +
    ylab('') +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank()))

  # print dest clustering pdf
  pdf_name = '../results/uc_irvine_forum/uc_irvine_forum_dest.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=4)
  embed_fonts(pdf_name, outfile=pdf_name)

  return()

}

experimentEigenvectorSweep = function(){

  cat('Eigenvector sweep\n')

  t0 = Sys.time()

  n = 5

  edge_list = matrix(c(
    1,3,
    1,2,
    4,3,
    4,2,
    4,5,
    3,5,

    4,6,
    5,6,
    5,7,

    6,8,
    6,9,
    6,7,
    8,9,
    8,7,
    8,10,
    9,10,
    7,10
  ), ncol=2, byrow=TRUE)


  G = adj_from_edge_list(edge_list)
  G = G + t(G)

  set.seed(6561)
  truth = c(rep(1,n),rep(2,n))

  plotcoord = data.frame(
    x = c(0,  0,  1  ,1,  2,       2,  3,  3,  2.5,4  ),
    y = c(3.5,2,  4  ,2.5,3,       1.5,2,  1,  0,  0.5)
  )
  plotcoord$cols = as.character(truth)
  colnames(plotcoord) = c("X1","X2",'cols')
  edges = data.frame(plotcoord[edge_list[,1],1:2], plotcoord[edge_list[,2],1:2])
  colnames(edges) = c("X1","Y1","X2","Y2")

  # run algorithm
  LapEmb = runLapEmb(G, 2, 'rw')
  vect = LapEmb$vects[,2]
  scores = sweepScores(G, vect, 'n_cut', verbose=FALSE)
  plotcoord$index = 1:nrow(plotcoord)


  # make network plot
  dividing_line = data.frame(x=c(0.5,3), y=c(0.5,3.4))
  pl = ggplot() +
    geom_segment(data=edges, aes(x=X1, y=Y1, xend = X2, yend = Y2), color='grey', size=0.5) +
    geom_line(data=dividing_line, aes(x=x,y=y), linetype='dashed') +
    geom_point(aes(x=X1,y=X2,col=cols), data=plotcoord, alpha=1, size=15) +
    scale_color_manual(labels=c(' Cluster 1',' Cluster 2'), values=brewer.pal(3,'Set3')[c(1,3)]) +
    geom_text(data=plotcoord, aes(x=X1,y=X2,label=index), family='CM Roman', size=7) +
    ggtitle('') +
    xlab('') +
    ylab('') +
    coord_cartesian(xlim=c(-0.5,4.5), ylim=c(-0.5,4.5)) +
    guides(colour = guide_legend(override.aes = list(size=13))) +
    (theme_diss() + theme(axis.ticks=element_blank(), axis.text=element_blank(), plot.margin = margin(l=-20,t=-5,b=-20,r=20),
                          legend.position = c(0.842,0.82), legend.title=element_blank(), legend.text=element_text(size=20),
                          panel.grid = element_blank(), panel.border=element_blank()))

  # print network plot pdf
  pdf_name = '../results/eigenvector_sweep/eigenvector_sweep_network.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)



  minimiser = data.frame(x=which.min(scores), y=min(scores))
  cat('\tSplitting point is i =', which.min(scores), '\n')
  cat('\tP1 =', (1:nrow(G))[truth==1], '\n')
  cat('\tP2 =', (1:nrow(G))[truth==2], '\n')

  # make sweep plot
  df = data.frame(x=1:length(scores), y=rev(scores))
  pl = ggplot(data=df, mapping=aes(x=x,y=y)) +
    geom_point(size=3) +
    geom_line() +
    geom_point(data=minimiser, aes(x=x,y=y), size=10, shape=1) +
    ggtitle('') +
    xlab(as.expression(bquote('Splitting point, '*italic(i)))) +
    ylab('Ncut score') +
    coord_cartesian(ylim=c(0,0.6)) +
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6)) +
    scale_x_continuous(breaks=1:9) +
    (theme_diss() + theme(plot.margin = margin(t=-5,r=1,b=5)))

  # print sweep plot to pdf
  pdf_name = '../results/eigenvector_sweep/eigenvector_sweep_scores.pdf'
  ggsave(filename=pdf_name, plot=pl, width=7, height=5)
  embed_fonts(pdf_name, outfile=pdf_name)

  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')

  return()
}

experimentTiming = function(){

  t0 = Sys.time()

  ns = c(100,1000,10000)
  ps = c(0.0001,0.001,0.01,0.1)
  motif_names = c('Ms','Md','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','M13')

  column_names = c('n','p','motif','t')
  results = dataFrameNames(column_names, n_rows = length(ns)*length(ps)*length(motif_names))
  row_number = 1

  for(n in ns){

    for(p in ps){

      # make graph
      G = sampleDSBM(n, matrix(p))

      for(motif_name in motif_names){

        #time MAM
        start_time = Sys.time()
        M = motifAdjacency(G, motif_name, 'func')
        end_time = Sys.time()
        time_diff = difftime(end_time, start_time, units="secs")
        time_diff = sigfig(time_diff, 2)

        # enter results
        results$n[row_number] = n
        results$p[row_number] = p
        results$motif[row_number] = motif_name
        results$t[row_number] = time_diff

        # progress
        cat('n = ',n,', p = ',p,', motif = ',motif_name,', time = ',time_diff,'\n', sep='')
        row_number = row_number + 1
      }
    }
  }
  print(results)

  # write results
  for(i in 1:length(ns)){

    n = ns[i]
    data_string = '\\cellcolor[HTML]{E9E9E9} \\smash{\\raisebox{0.7pt}{$p$}} & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_\\mathrm{s}$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_\\mathrm{d}$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_1$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_2$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_3$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_4$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_5$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_6$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_7$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_8$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_9$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_{10}$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_{11}$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_{12}$ & \\cellcolor[HTML]{E9E9E9} $\\ca{M}_{13}$ \\\\ \\hline'

    for(j in 1:length(ps)){

      p = ps[j]
      p_print = formatC(signif(p,digits=1), digits=1,format="fg", flag="#")
      data_string = paste(data_string, '\\cellcolor[HTML]{E9E9E9}', p_print)

      for(k in 1:length(motif_names)){

        motif_name = motif_names[k]
        data_string = paste(data_string, '& ',results$t[(results$n ==n ) & (results$p == p) & (results$motif == motif_name)])
      }

      data_string = paste(data_string, '\\\\ \\hline')
    }

    file_name = paste( '../results/timing/timing_n_',n,'.txt', sep='')
    write(data_string, file_name)
  }

  # overall time
  cat('\t',capture.output(Sys.time()-t0), sep='')
  cat('\n\n')
  return()
}

# Graph Data
#########################################################################################################################################
readGraphData = function(graphname){

  # Reads the sample data available for this project, as a list containing (weighted) adjacency matrix and other attributes
  # TODO use fread

  if(graphname == 'us_canada_airlines'){ans = read_us_canada_airlines()}
  else if(graphname == 'karate_club'){ans = read_karate_club()}
  else if(graphname == 'us_airports'){ans = read_us_airports()}
  else if(graphname == 'polblogs'){ans = read_polblogs()}
  else if(graphname == 'hessen_roads'){ans = read_hessen_roads()}
  else if(graphname == 'uk_roads'){ans = read_uk_roads()}
  else if(graphname == 'us_migration'){ans = read_us_migration()}
  else if(graphname == 'oldenburg'){ans = read_oldenburg()}
  else if(graphname == 'american_revolution'){ans = read_american_revolution()}
  else if(graphname == 'faculty_business'){ans = read_faculty_business()}
  else if(graphname == 'faculty_cs'){ans = read_faculty_cs()}
  else if(graphname == 'faculty_history'){ans = read_faculty_history()}
  else if(graphname == 'clements_pollination'){ans = read_clements_pollination()}
  else if(graphname == 'uc_irvine_forum'){ans = read_uc_irvine_forum()}
  else if(graphname == 'languages'){ans = read_languages()}

  else(stop('Invalid graph name'))
  return(ans)

}

read_us_canada_airlines = function(){

  df = read.csv('../datasets/us_canada_airlines/us_canada_airlines.csv',sep = ' ')
  cities = unname(unlist(fread('../datasets/us_canada_airlines/us_canada_airlines_names.csv',sep = '\n', header=FALSE,data.table=FALSE)))
  gr = graph.data.frame(df)

  ans = list()
  ans$adj = as_adj(gr, attr='similarity')
  ans$cities = cities
  return(ans)
}

read_karate_club = function(){
  gr = read_graph('../datasets/karate_club/karate_club.gml', format = 'gml')

  ans = list()
  ans$adj = as_adj(gr)
  return(ans)
}

read_us_airports = function(){

  df = read.csv('../datasets/us_airports/us_airports.csv',sep = ' ')
  all_airports = read.csv('../datasets/us_airports/us_airports_names.csv',sep = ' ')
  gr = graph.data.frame(df)
  vertex_list = sort(unique(as.numeric(vertex_attr(gr,'name'))))
  used_airports = as.character((all_airports[vertex_list,])$airport_code)

  ans = list()
  ans$adj = as_adj(gr, attr='similarity')
  ans$airports = used_airports
  return(ans)
}

read_polblogs = function(){

  gr = read_graph('../datasets/polblogs/polblogs.gml',format='gml')

  ans = list()
  ans$adj = as_adj(gr)
  ans$labels = vertex_attr(gr,'value')
  ans$source = vertex_attr(gr,'source')
  return(ans)
}

read_hessen_roads = function(){

  df = read.csv('../datasets/hessen_roads/hessen_roads.csv',sep = '\t')
  gr = graph.data.frame(df)

  ans = list()
  ans$adj = as_adj(gr, attr='Capacity')
  ans$length = edge_attr(gr,'Length')
  return(ans)
}

read_uk_roads = function(){

  # TODO
  df = fread('../datasets/uk_roads/uk_roads.csv',data.table=FALSE,nrows = 100000)
  df = df[df$`Region Name (GO)` == 'Wales',]
  print(df[1,])

  gr = graph.data.frame(df[])

  ans = list()
  ans = df
  return(ans)
}

read_us_migration = function(){

  df = fread('../datasets/us_migration/us_migration.txt',data.table=FALSE,sep = ' ',colClasses=c('character','character','integer'))
  colnames(df) = c('county_1','county_2','flow')

  # remove AK, HI
  df = df[!(substr(df$county_1,1,2) %in% c('02','15')),]
  df = df[!(substr(df$county_2,1,2) %in% c('02','15')),]
  gr = graph.data.frame(df)
  fips = vertex_attr(gr,'name')

  # format counties
  counties = us_counties('2000-04-01', resolution='low')
  counties = ms_simplify(counties, keep=0.4, drop_null_geometries=FALSE, keep_shapes=TRUE)
  # remove AK, HI, PR
  counties = counties[!(counties$state_abbr %in% c('AK','HI','PR')),]
  # remove merged counties
  counties = counties[!(counties$fips %in% setdiff(counties$fips,fips)),]
  counties = counties[order(counties$fips),]

  # format states
  states = us_states('2000-04-01', resolution='low')
  states = ms_simplify(states, keep=0.5, drop_null_geometries=FALSE, keep_shapes=TRUE)
  states = states[!(states$state_abbr %in% c('AK','HI')),]

  # add centroids
  state_label_coords = st_as_sf(gCentroid(as(states, Class='Spatial'), byid=TRUE))
  state_label_coords = st_transform(state_label_coords, crs=2163)
  state_label_coords = st_coordinates(state_label_coords)

  # adjust centroid positions
  state_label_coords[states$abbr_name == 'ID',2] = state_label_coords[states$abbr_name == 'ID',2] - 10e4
  state_label_coords[states$abbr_name == 'CA',2] = state_label_coords[states$abbr_name == 'CA',2] - 8e4
  state_label_coords[states$abbr_name == 'LA',1] = state_label_coords[states$abbr_name == 'LA',1] - 6e4
  state_label_coords[states$abbr_name == 'LA',2] = state_label_coords[states$abbr_name == 'LA',2] + 8e4
  state_label_coords[states$abbr_name == 'KY',1] = state_label_coords[states$abbr_name == 'KY',1] + 5e4
  state_label_coords[states$abbr_name == 'MI',1] = state_label_coords[states$abbr_name == 'MI',1] + 7e4
  state_label_coords[states$abbr_name == 'MI',2] = state_label_coords[states$abbr_name == 'MI',2] - 10e4
  state_label_coords[states$abbr_name == 'FL',1] = state_label_coords[states$abbr_name == 'FL',1] + 7e4
  state_label_coords[states$abbr_name == 'VT',2] = state_label_coords[states$abbr_name == 'VT',2] + 5e4
  state_label_coords[states$abbr_name == 'VT',1] = state_label_coords[states$abbr_name == 'VT',1] - 1.5e4
  state_label_coords[states$abbr_name == 'NH',1] = state_label_coords[states$abbr_name == 'NH',1] + 1e4
  state_label_coords[states$abbr_name == 'NH',2] = state_label_coords[states$abbr_name == 'NH',2] - 5.5e4
  state_label_coords[states$abbr_name == 'MA',1] = state_label_coords[states$abbr_name == 'MA',1] + 20e4
  state_label_coords[states$abbr_name == 'MA',2] = state_label_coords[states$abbr_name == 'MA',2] + 11e4
  state_label_coords[states$abbr_name == 'RI',1] = state_label_coords[states$abbr_name == 'RI',1] + 15e4
  state_label_coords[states$abbr_name == 'RI',2] = state_label_coords[states$abbr_name == 'RI',2] - 10e4
  state_label_coords[states$abbr_name == 'NJ',1] = state_label_coords[states$abbr_name == 'NJ',1] + 20e4
  state_label_coords[states$abbr_name == 'CT',2] = state_label_coords[states$abbr_name == 'CT',2] + 0.7e4
  state_label_coords[states$abbr_name == 'DE',1] = state_label_coords[states$abbr_name == 'DE',1] + 20e4
  state_label_coords[states$abbr_name == 'DE',2] = state_label_coords[states$abbr_name == 'DE',2] + 3e4
  state_label_coords[states$abbr_name == 'MD',1] = state_label_coords[states$abbr_name == 'MD',1] + 30e4
  state_label_coords[states$abbr_name == 'MD',2] = state_label_coords[states$abbr_name == 'MD',2] - 5e4
  state_label_coords[states$abbr_name == 'DC',1] = state_label_coords[states$abbr_name == 'DC',1] + 27e4
  state_label_coords[states$abbr_name == 'DC',2] = state_label_coords[states$abbr_name == 'DC',2] - 15e4


  # set centroid label sizes
  state_label_sizes = rep(2, nrow(states))
  state_label_sizes[states$abbr_name == 'VT'] = 1.8
  state_label_sizes[states$abbr_name == 'NH'] = 1.8
  state_label_sizes[states$abbr_name == 'CT'] = 1.8



  # add lines for state labels
  state_label_line_x = rep(NA, nrow(states))
  state_label_line_xend = rep(NA, nrow(states))
  state_label_line_y = rep(NA, nrow(states))
  state_label_line_yend = rep(NA, nrow(states))

  state_label_line_x[states$abbr_name == 'RI'] = state_label_coords[states$abbr_name == 'RI',1] - 5e4
  state_label_line_xend[states$abbr_name == 'RI'] = state_label_coords[states$abbr_name == 'RI',1] - 15e4
  state_label_line_y[states$abbr_name == 'RI'] = state_label_coords[states$abbr_name == 'RI',2]
  state_label_line_yend[states$abbr_name == 'RI'] = state_label_coords[states$abbr_name == 'RI',2] + 9e4

  state_label_line_x[states$abbr_name == 'NJ'] = state_label_coords[states$abbr_name == 'NJ',1] - 5e4
  state_label_line_xend[states$abbr_name == 'NJ'] = state_label_coords[states$abbr_name == 'NJ',1] - 18e4
  state_label_line_y[states$abbr_name == 'NJ'] = state_label_coords[states$abbr_name == 'NJ',2]
  state_label_line_yend[states$abbr_name == 'NJ'] = state_label_coords[states$abbr_name == 'NJ',2]

  state_label_line_x[states$abbr_name == 'DE'] = state_label_coords[states$abbr_name == 'DE',1] - 6e4
  state_label_line_xend[states$abbr_name == 'DE'] = state_label_coords[states$abbr_name == 'DE',1] - 20e4
  state_label_line_y[states$abbr_name == 'DE'] = state_label_coords[states$abbr_name == 'DE',2]
  state_label_line_yend[states$abbr_name == 'DE'] = state_label_coords[states$abbr_name == 'DE',2] - 5e4

  state_label_line_x[states$abbr_name == 'MD'] = state_label_coords[states$abbr_name == 'MD',1] - 7e4
  state_label_line_xend[states$abbr_name == 'MD'] = state_label_coords[states$abbr_name == 'MD',1] - 17e4
  state_label_line_y[states$abbr_name == 'MD'] = state_label_coords[states$abbr_name == 'MD',2]
  state_label_line_yend[states$abbr_name == 'MD'] = state_label_coords[states$abbr_name == 'MD',2]

  state_label_line_x[states$abbr_name == 'DC'] = state_label_coords[states$abbr_name == 'DC',1] - 6e4
  state_label_line_xend[states$abbr_name == 'DC'] = state_label_coords[states$abbr_name == 'DC',1] - 27e4
  state_label_line_y[states$abbr_name == 'DC'] = state_label_coords[states$abbr_name == 'DC',2]
  state_label_line_yend[states$abbr_name == 'DC'] = state_label_coords[states$abbr_name == 'DC',2] + 15e4

  state_label_line_x[states$abbr_name == 'MA'] = state_label_coords[states$abbr_name == 'MA',1] - 7e4
  state_label_line_xend[states$abbr_name == 'MA'] = state_label_coords[states$abbr_name == 'MA',1] - 20e4
  state_label_line_y[states$abbr_name == 'MA'] = state_label_coords[states$abbr_name == 'MA',2]
  state_label_line_yend[states$abbr_name == 'MA'] = state_label_coords[states$abbr_name == 'MA',2] - 11e4

  # make grey/white patchwork
  state_shading = rep('g', nrow(states))
  state_shading[states$abbr_name %in% c('AZ','WY','AL','KY','NC','MD')] = 'r'
  state_shading[states$abbr_name %in% c('WA','NV','NM','ND','KS','IA','AR','MI','FL','VA','SC','PA','CT','ME','VT')] = 'o'
  state_shading[states$abbr_name %in% c('OR','UT','MT','NE','OK','IL','MN','TN','LA','OH','NJ','MA')] = 'b'

  # add label data to states object
  states = cbind(states, state_label_coords)
  states$label_sizes = state_label_sizes
  states$label_line_x = state_label_line_x
  states$label_line_xend = state_label_line_xend
  states$label_line_y = state_label_line_y
  states$label_line_yend = state_label_line_yend
  states$shading = state_shading

  ans = list()
  ans$adj = unname(as_adj(gr, attr='flow'))
  ans$counties = counties
  ans$states = states
  return(ans)
}

read_oldenburg = function(){

  df = fread('../datasets/oldenburg/oldenburg_edges.txt',data.table=FALSE,sep = ' ')
  df = df[,2:4]
  gr = graph.data.frame(df)

  ans = list()
  ans$adj = unname(as_adj(gr, attr='length'))
  return(ans)
}

read_american_revolution = function(){

  df = fread('../datasets/american_revolution/american_revolution.csv',data.table=FALSE,sep = ',')
  gr = graph.data.frame(df)

  ans = list()
  ans$adj = unname(as_adj(gr))
  return(ans)
}

read_faculty_business = function(){

  # read edge list
  edge_list = fread('../datasets/faculty_hiring/business/Business_edgelist.txt',data.table=FALSE,sep = '\t')
  edge_list = edge_list[,1:2]
  unis = read.csv('../datasets/faculty_hiring/business/Business_vertexlist.txt', sep='\t')

  # adj
  adj = adj_from_edge_list(edge_list)

  ans = list()
  ans$adj = adj
  ans$unis = unis
  return(ans)

}

read_faculty_cs = function(){

  df = fread('../datasets/faculty_hiring/cs/ComputerScience_edgelist.txt',data.table=FALSE,sep = '\t')
  unis = read.csv('../datasets/faculty_hiring/cs/ComputerScience_vertexlist.txt', sep='\t')
  gr = graph.data.frame(df,vertices = unis$u)
  ord = order(as.numeric(vertex_attr(gr,'name')))
  named = 1:(length(ord)-1)

  ans = list()
  ans$adj = drop0_killdiag(unname(as_adj(gr))[ord,ord])[named,named]
  ans$unis = unis[named,]
  return(ans)

}

read_faculty_history = function(){

  df = fread('../datasets/faculty_hiring/history/History_edgelist.txt',data.table=FALSE,sep = '\t')
  unis = read.csv('../datasets/faculty_hiring/history/History_vertexlist.txt', sep='\t')
  gr = graph.data.frame(df,vertices = unis$u)
  ord = order(as.numeric(vertex_attr(gr,'name')))
  named = 1:(length(ord)-1)

  ans = list()
  ans$adj = drop0_killdiag(unname(as_adj(gr))[ord,ord])[named,named]
  ans$unis = unis[named,]
  return(ans)

}

read_clements_pollination = function(){

  df = fread('../datasets/clements_pollination/clements_1923.csv',data.table=FALSE,sep = ',')

  Pl_ge = df[2,5:100]
  Pl_sp = df[3,5:100]
  Pol_ge = df[5:279,2]
  Pol_sp = df[5:279,3]

  pollinators = paste(Pol_ge, Pol_sp)
  plants = paste(Pl_ge, Pl_sp)

  small_adj = matrix(as.numeric(as.matrix(df[5:279,5:100])), nrow=275)
  adj = matrix(0, nrow=371, ncol=371)
  adj[1:275,276:371] = small_adj

  ans = list()
  ans$adj = drop0(adj)
  ans$pollinators = pollinators
  ans$plants = plants

  return(ans)
}

read_uc_irvine_forum = function(){

  df = fread('../datasets/uc_irvine_forum/uc_irvine_forum.txt',data.table=FALSE,sep = ' ')
  gr = graph.data.frame(df)

  ans = list()
  ans$adj = drop0_killdiag(unname(as_adj(gr, attr='V3')))
  return(ans)
}

read_languages = function(){

  # the master data are kept in adj, countries and languages.
  # all others are formatted to match these.
  # big populations and connected component are extracted at the end.

  # read edge list
  edge_list = fread('../datasets/languages/unicodelang/out.unicodelang',data.table=FALSE,sep = '\t',header=FALSE)
  edge_list$V2 = edge_list$V2 + 254

  # adj
  adj = adj_from_edge_list(edge_list)

  # countries
  countries = unname(unlist(fread('../datasets/languages/unicodelang/ent.unicodelang.coutry.code',data.table=FALSE,sep = '\t',header=FALSE)))
  countries[156] = 'NA' # Nigeria formats as NA not 'NA'

  # languages
  languages = unname(unlist(fread('../datasets/languages/unicodelang/ent.unicodelang.language.code',data.table=FALSE,sep = '\t',header=FALSE)))






  # country names
  country_names = fread('../datasets/languages/languages_country_names.txt',data.table=FALSE,sep = ',',header=TRUE)
  country_names[155,2] = 'NA' # Nigeria

  new_country_names = c()
  for(i in 1:length(countries)){
    country = countries[i]
    if((country %in% country_names$Code)){
      new_country_names = c(new_country_names, country_names$Name[country_names$Code == country])
    }
    else{
      new_country_names = c(new_country_names, country)
    }
  }
  country_names = new_country_names


  # populations
  populations = fread('../datasets/languages/languages_populations.txt',data.table=FALSE,sep = ',',header=FALSE)
  colnames(populations) = c('country','population')
  populations$country[161] = 'NA' # Nigeria

  new_populations = c()
  for(i in 1:length(countries)){
    country = countries[i]
    if((country %in% populations$country)){
      new_populations = c(new_populations, populations$population[populations$country == country])
    }
    else{
      new_populations = c(new_populations, 0)
    }
  }
  populations = new_populations


  # language names
  language_names = fread('../datasets/languages/languages_language_names2.txt',data.table=FALSE,sep = '\t',header=TRUE)

  new_language_names = c()
  for(i in 1:length(languages)){
    language = languages[i]
    if((language %in% language_names$Id)){
      new_language_names = c(new_language_names, language_names$Ref_Name[language_names$Id == language])
    }
    else if((language %in% language_names$Part1)){
      new_language_names = c(new_language_names, language_names$Ref_Name[language_names$Part1 == language])
    }
    else{
      new_language_names = c(new_language_names, language)
    }
  }
  language_names = new_language_names





  # multiply adj by populations
  for(i in 1:length(countries)){
    adj[i,] = adj[i,] * new_populations[i]
  }

  # number of language speakers
  language_speakers = apply(adj[1:length(countries),(length(countries)+1):(length(countries)+length(languages))],2,sum)





  # find countries and languages with large populations
  big_countries = (1:length(countries))[new_populations > 1e6]
  big_languages = (1:length(languages))[language_speakers > 1e6]
  big_countries_and_languages = c(big_countries,(big_languages+length(countries)))

  # include only countries and languages with large populations
  adj = adj[big_countries_and_languages,big_countries_and_languages]
  countries = countries[big_countries]
  languages = languages[big_languages]
  populations = populations[big_countries]
  country_names = country_names[big_countries]
  language_names = language_names[big_languages]
  language_speakers = language_speakers[big_languages]





  # get largest connected component
  comps_countries_and_languages = largestComponent(adj)
  comps_countries = comps_countries_and_languages[comps_countries_and_languages <= length(countries)]
  comps_languages = comps_countries_and_languages[comps_countries_and_languages >  length(countries)] - length(countries)

  # print countries and languages not on largest connected component
  not_comps_countries_and_languages = setdiff(1:nrow(adj), comps_countries_and_languages)
  not_comps_countries = not_comps_countries_and_languages[not_comps_countries_and_languages <= length(countries)]
  not_comps_languages = not_comps_countries_and_languages[not_comps_countries_and_languages >  length(countries)] - length(countries)
  cat('\tTerritories not on largest component: ', paste(c(country_names[not_comps_countries]), ', ', sep=''), '\n', sep='')
  cat('\tLanguages not on largest component: ', paste(language_names[not_comps_languages], ', ',sep=''), '\n', sep='')

  # include only countries and languages on largest component
  adj = adj[comps_countries_and_languages,comps_countries_and_languages]
  countries = countries[comps_countries]
  languages = languages[comps_languages]
  populations = populations[comps_countries]
  country_names = country_names[comps_countries]
  language_names = language_names[comps_languages]
  language_speakers = language_speakers[comps_languages]






  # return ans
  ans = list()
  ans$adj = adj
  ans$countries = countries
  ans$languages = languages
  ans$populations = populations
  ans$country_names = country_names
  ans$language_names = language_names
  ans$language_speakers = language_speakers
  return(ans)
}
