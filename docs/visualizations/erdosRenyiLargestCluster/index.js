// Function to generate Erdős–Rényi network (directed)
function generateErdosRenyi(n = 100, p = 0.02) { // Lowered p for less connectivity
  const nodes = [];
  for (let i = 0; i < n; i++) {
    nodes.push({ id: `node_${i}`, name: `Node ${i}`, concentration: Math.random() * 10 + 1 }); // Random equilibrium concentration
  }
  const edges = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i !== j && Math.random() < p) {
        const rate = Math.random() * 5 + 0.1; // Random reaction rate
        edges.push({ source: nodes[i].id, target: nodes[j].id, rate });
      }
    }
  }
  return { nodes, edges };
}

// Function to find largest connected component (undirected)
function getLargestComponent(nodes, edges) {
  const adj = {};
  nodes.forEach(n => adj[n.id] = []);
  edges.forEach(e => {
    adj[e.source].push(e.target);
    adj[e.target].push(e.source); // Undirected for connectivity
  });
  const visited = new Set();
  let maxComponent = [];
  nodes.forEach(node => {
    if (!visited.has(node.id)) {
      const component = [];
      const stack = [node.id];
      while (stack.length) {
        const curr = stack.pop();
        if (!visited.has(curr)) {
          visited.add(curr);
          component.push(curr);
          adj[curr].forEach(neigh => {
            if (!visited.has(neigh)) stack.push(neigh);
          });
        }
      }
      if (component.length > maxComponent.length) maxComponent = component;
    }
  });
  const componentSet = new Set(maxComponent);
  const filteredNodes = nodes.filter(n => componentSet.has(n.id));
  const filteredEdges = edges.filter(e => componentSet.has(e.source) && componentSet.has(e.target));
  return { nodes: filteredNodes, edges: filteredEdges };
}

// Generate random graph instead of loading SDF
console.log("Generating Erdős–Rényi network...");
let { nodes: molecules, edges } = generateErdosRenyi(100, 0.01);
console.log(`Generated ${molecules.length} nodes and ${edges.length} edges.`);

// Filter to largest component
const filtered = getLargestComponent(molecules, edges);
molecules = filtered.nodes;
edges = filtered.edges;
console.log(`After filtering: ${molecules.length} nodes and ${edges.length} edges.`);

// Set up SVG
const width = window.innerWidth;
const height = window.innerHeight;
const svg = d3.select('body').append('svg')
  .attr('width', '100%')
  .attr('height', '100%')
  .style('position', 'fixed')
  .style('top', 0)
  .style('left', 0)
  .style('width', '100vw')
  .style('height', '100vh');

// Add arrow marker for directed edges
svg.append('defs').append('marker')
  .attr('id', 'arrow')
  .attr('viewBox', '0 -5 10 10')
  .attr('refX', 13)
  .attr('refY', 0)
  .attr('markerWidth', 6)
  .attr('markerHeight', 6)
  .attr('orient', 'auto')
  .append('path')
  .attr('d', 'M0,-5L10,0L0,5')
  .attr('fill', '#999');

// Create simulation
const simulation = d3.forceSimulation(molecules)
  .force('link', d3.forceLink(edges).id(d => d.id).distance(d => 100 / d.rate)) // Shorter distance for higher rates
  .force('charge', d3.forceManyBody())
  .force('center', d3.forceCenter(width / 2, height / 2))
  .alphaDecay(0) // Slower decay for longer animation
  .alphaMin(0.05); // Higher min alpha to keep constant movement

// Draw links (with arrows)
const link = svg.append('g')
  .selectAll('line')
  .data(edges)
  .enter().append('line')
  .attr('stroke', '#999')
  .attr('stroke-opacity', d => Math.min(d.rate / 2, 1)) // Opacity proportional to rate
  .attr('stroke-width', d => 2 ) // Increased thickness
  .attr('marker-end', 'url(#arrow)');

// Draw nodes (size based on concentration)
const node = svg.append('g')
  .selectAll('circle')
  .data(molecules)
  .enter().append('circle')
  .attr('r', d => Math.sqrt(d.concentration) * 3) // Increased size
  .attr('fill', '#69b3a2')
  .call(d3.drag()
    .on('start', event => {
      if (!event.active) simulation.alphaTarget(0.3).restart();
      event.subject.fx = event.subject.x;
      event.subject.fy = event.subject.y;
    })
    .on('drag', event => {
      event.subject.fx = event.x;
      event.subject.fy = event.y;
    })
    .on('end', event => {
      if (!event.active) simulation.alphaTarget(0);
      event.subject.fx = null;
      event.subject.fy = null;
    }));

// Update positions
simulation.on('tick', () => {
  link
    .attr('x1', d => d.source.x)
    .attr('y1', d => d.source.y)
    .attr('x2', d => d.target.x)
    .attr('y2', d => d.target.y);

  node
    .attr('cx', d => d.x = Math.max(10, Math.min(width - 10, d.x)))
    .attr('cy', d => d.y = Math.max(10, Math.min(height - 10, d.y)));
});
