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

<h2>Real-world Datasets</h2>

<p>
The real-world datasets used in this project are collected from the KONECT network repository:
</p>

<p>
<a href="http://konect.cc/networks/" target="_blank">http://konect.cc/networks/</a>
</p>

<p>
KONECT is a publicly available repository providing a large collection of real-world networks, 
including bipartite graphs from various application domains.
</p>
 
<h2>Usage</h2>

<h3>1. Bipartite Graph Decomposition</h3>

<p>
To decompose a bipartite graph, use the following command:
</p>

<pre><code>decompose data/WC</code></pre>

<p>
<code>data/WC</code> specifies the path to the dataset.<br>
The graph will be read and decomposed according to the implemented bipartite core definition.
</p>

<h3>2. Community Search</h3>

<p>
The command format for community search is:
</p>

<pre><code>community-search data/WC 1 2 3</code></pre>

<p>Where:</p>
<ul>
  <li><code>data/WC</code> : path to the dataset</li>
  <li><code>1</code> : query vertex</li>
  <li><code>2</code> : parameter α (alpha)</li>
  <li><code>3</code> : parameter β (beta)</li>
</ul>

<p>
This command searches for the community that contains the given query vertex under the specified <code>(α, β)</code> parameters.
</p>

