
use petgraph::prelude::*;
use petgraph::visit::NodeIndexable;
use std::collections::HashSet;
fn bron_kerbosch(
    graph: &Graph<(), (), Undirected>,
    mut r: HashSet<NodeIndex>,
    mut p: HashSet<NodeIndex>,
    mut x: HashSet<NodeIndex>,
    cliques: &mut Vec<HashSet<NodeIndex>>,
) {
    if p.is_empty() && x.is_empty() {
        cliques.push(r.clone());
        return;
    }

    let mut p_candidates = p.clone();
    p_candidates.retain(|&v| {
        let mut neighbors = graph.neighbors(v);
        neighbors.all(|n| p.contains(&n))
    });

    let p_cloned = p.clone();
    for v in p_cloned {
        let mut neighbors = graph.neighbors(v);
        let mut r_new = r.clone();
        r_new.insert(v);

        let p_new = p.intersection(&neighbors.clone().collect::<HashSet<_>>()).cloned().collect();
        let x_new = x.intersection(&neighbors.collect::<HashSet<_>>()).cloned().collect();

        bron_kerbosch(graph, r_new, p_new, x_new, cliques);

        p.remove(&v);
        x.insert(v);
    }
}

fn find_cliques(graph: &Graph<(), (), Undirected>) -> Vec<HashSet<NodeIndex>> {
    let mut cliques = Vec::new();
    let n = graph.node_count();
    let mut all_nodes = (0..n).map(NodeIndex::new).collect::<HashSet<NodeIndex>>();
    let mut r = HashSet::new();

    let mut x = HashSet::new();

    bron_kerbosch(graph, r, all_nodes, x, &mut cliques);
    cliques
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_find_cliques() {
        // Create an undirected graph
        let mut graph = Graph::<(), (), Undirected>::new_undirected();
        let a = graph.add_node(());
        let b = graph.add_node(());
        let c = graph.add_node(());
        let d = graph.add_node(());
        let e = graph.add_node(());

        graph.add_edge(a, b, ());
        graph.add_edge(a, c, ());
        graph.add_edge(b, c, ());
        graph.add_edge(b, d, ());
        graph.add_edge(c, d, ());
        graph.add_edge(c, e, ());
 
        // Find cliques in the graph
        let cliques = find_cliques(&graph);

        // Print the cliques
   
        for clique in cliques {
            let clique_nodes = clique.iter().collect::<Vec<_>>();
            println!("{:?}", clique_nodes);
        }
    }
}