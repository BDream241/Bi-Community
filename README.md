This  provides a C++ implementation for bipartite graph community search.

<h2>Overview</h2>
We implement a set of tools for processing bipartite graphs, including graph decomposition, index construction, bi-core search, and community search.
The system supports querying  searching (\alpha, \beta)-communities containing a specified query vertex with given parameters.

<h2>Code Structure</h2>

<p>The main components of the code are organized as follows:</p>

<h3>utility.h</h3>
<p>
The base file of the project.<br>
It contains common utility functions, shared data structures, and basic definitions used throughout the codebase.
</p>

<h3>bigraph.h</h3>
<p>
This file is responsible for reading bipartite graph data and performing bipartite graph decomposition.<br>
It handles graph construction from input files and prepares the graph for subsequent index building and query operations.
</p>

<h3>abIndex.h</h3>
<p>
The index construction module.<br>
This file implements the index used to support efficient bi-core search and community search on bipartite graphs.
</p>
 
The command format for decomposing bipartite graphs: decompose data/WC

.. /data/WC is the address of the dataset

The command format of community search: community-search data/WC 1 2 3

1, 2, 3 respectively represent the query vertex, alpha and beta

The command format of bi-core search: bi-core search data/WC 2 3

2 and 3 respectively represent alpha and beta
