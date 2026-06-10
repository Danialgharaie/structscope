# Interactions and Contact Graphs

`structscope` computes advanced chemical and geometric interactions directly from coordinates, and integrates them into its contact graph representations.

## 1. Interaction Detection Logic

All interaction detectors reside in `structscope-features` under `interactions.rs` and are computed based on atom names and distance/angle rules:

- **Disulfide Bonds**: Formed between two cysteine sidechain sulfur atoms (`CYS` `SG`) within $2.5\text{ Å}$ distance.
- **Salt Bridges**: Formed between sidechain acidic oxygens (on `ASP`/`GLU`) and sidechain basic nitrogens (on `LYS`/`ARG`/`HIS`) within $4.0\text{ Å}$ distance.
- **Hydrogen Bonds (Polar Contacts)**: Formed between a nitrogen/oxygen donor-acceptor pair on different residues within $2.4\text{ to }3.5\text{ Å}$.
- **Cation-Pi Interactions**: Formed between the centroid of an aromatic ring (`PHE`/`TYR`/`TRP`) and a basic sidechain nitrogen within $6.0\text{ Å}$.
- **Pi-Pi Stacking**: Formed between the centroids of two aromatic rings within $5.5\text{ Å}$:
  - **Parallel Stacking**: Ring normal angle $< 30^\circ$ or $> 150^\circ$.
  - **Perpendicular Stacking**: Ring normal angle between $60^\circ$ and $90^\circ$.
- **Hydrophobic Contacts**: Formed between sidechain aliphatic/aromatic carbons on different residues within $4.5\text{ Å}$.

---

## 2. Contact Graph Representation & Prioritization

Contact graphs are built in `structscope-graphs` using `petgraph`. When generating residue or interface contact graphs, the builder can ingest chemical interactions and merge them as prioritized edges.

### Decoupled Data Flow
To avoid a circular dependency between `structscope-features` (which needs graph definitions for feature metrics) and `structscope-graphs`, `structscope-graphs` defines a simple intermediate struct:

```rust
pub struct ChemicalInteraction {
    pub kind: String,
    pub res_id_a: String,
    pub res_id_b: String,
    pub distance: f64,
}
```

The CLI orchestrator parses the structures, computes features/interactions, maps them to `ChemicalInteraction` using residue ID lookups, and passes them to the graph builders.

### Prioritization Rules
When multiple potential interaction types overlap between two residues, they are merged into a single edge according to strict precedence rules. 
- **Backbone Covalent Adjacency (`covalent_adjacency`)**: Always preserved; never overwritten by a chemical contact.
- **Other Contacts**: Overwritten if a higher-priority interaction type is found, or if an interaction of the same type has a shorter distance:

| Priority | Edge Kind |
| :--- | :--- |
| **7** (Highest) | `disulfide` |
| **6** | `salt_bridge` |
| **5** | `hydrogen_bond` |
| **4** | `cation_pi` |
| **3** | `pi_pi_parallel` / `pi_pi_perpendicular` |
| **2** | `hydrophobic` |
| **1** (Lowest) | `distance_contact` / `interface_contact` |

---

## 3. Export Formats

`structscope` supports three serialization formats for contact graphs:

### GraphML
Standard XML representation of graphs.
- **Nodes**: Contain node metadata (e.g., `residue_name`, `seq_number`).
- **Edges**: Contain edge details (`kind`, `distance`).

### GML (Graph Modeling Language)
An easy-to-parse ascii representation of graphs:
```gml
graph [
  directed 0
  node [
    id 0
    label "1nkd:A:1:_"
    residue_name "MET"
    seq_number 1
  ]
  edge [
    source 0
    target 1
    distance 5.291
    kind "covalent_adjacency"
  ]
]
```

### JSON (Node-Link Format)
A JSON format compatible with web-based visualization tools (e.g., D3.js, cytoscape.js):
```json
{
  "nodes": [
    {
      "id": 0,
      "residue_id": "1nkd:A:1:_",
      "chain_id": "1nkd:A",
      "residue_name": "MET",
      "seq_number": 1
    }
  ],
  "links": [
    {
      "source": 0,
      "target": 1,
      "distance": 5.291,
      "kind": "covalent_adjacency"
    }
  ]
}
```
