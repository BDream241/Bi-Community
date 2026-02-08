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

<p>
Each dataset, after preprocessing, is packaged in the same format as <code>data-AR.zip</code>,
which contains three files:
</p>

<ul>
  <li>
    <code>graph.txt</code>: stores the basic statistics of the bipartite graph, including
    the number of upper-layer vertices, the number of lower-layer vertices, and the total
    number of edges.
  </li>
  <li>
    <code>edge.txt</code>: stores all edges of the bipartite graph.
  </li>
  <li>
    <code>query.txt</code>: stores all query instances. Each query consists of four integers,
    e.g., <code>42 1 13 13</code>, representing the query vertex <code>q</code>, whether
    <code>q</code> belongs to the upper layer (0 for no, 1 for yes), and the
    <code>&alpha;</code> and <code>&beta;</code> parameters of the target community,
    respectively.
  </li>
</ul>
 
<h2>Usage</h2>

<h3>Bi-Community Search</h3>

<p>
<code>BiCommunitySearch.exe</code> is the executable file of the proposed method.
After launching the program, users need to provide the dataset path and the query file path as input.
</p>

<p>
The command format for bi-community search is:
</p>

<pre><code>BiCommunitySearch.exe AR AR/query.txt</code></pre>

<p>Where:</p>
<ul>
  <li><code>AR</code> : path to the dataset</li>
  <li><code>AR/query.txt</code> : path to the query file</li>
</ul>

<p>
The program processes all queries listed in <code>query.txt</code> sequentially and reports
the index construction time as well as the search time for each query.
</p>

<p>
An example output is shown below:
</p>

<pre><code>Index construction time: 48.4028 s
community-0 search time: 0.15 ms
community-1 search time: 0.1091 ms
community-2 search time: 0.1083 ms
community-3 search time: 0.1194 ms
community-4 search time: 0.098 ms
</code></pre>

<p>
Here, <code>community-i</code> corresponds to the <em>i</em>-th query in the query file.
</p>


