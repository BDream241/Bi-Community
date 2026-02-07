This  provides a C++ implementation for bipartite graph community search.

<h3>Overview</h3>
 
utility.h is the base file

bigraph.h is used for reading and decomposing bipartite graph files

abIndex.h is the index build file

The command format for decomposing bipartite graphs: decompose data/WC

.. /data/WC is the address of the dataset

The command format of community search: community-search data/WC 1 2 3

1, 2, 3 respectively represent the query vertex, alpha and beta

The command format of bi-core search: bi-core search data/WC 2 3

2 and 3 respectively represent alpha and beta
