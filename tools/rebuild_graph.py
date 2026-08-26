"""Rebuild the graphify graph over first-party source only.

`graphify update .` rescans everything and pulls in the vendored HashRF/tqDist C++
build trees plus mkdocs' bundled JS under site/, taking the graph from ~160 useful
nodes to 3600 mostly-irrelevant ones. This filters the detect output first.
"""

import json
import os
import sys
from pathlib import Path

from graphify.analyze import god_nodes, suggest_questions, surprising_connections
from graphify.build import build_from_json
from graphify.cluster import cluster, score_all
from graphify.detect import detect
from graphify.export import to_json
from graphify.extract import collect_files, extract
from graphify.report import generate

EXCLUDE = (
    "pear_ebi/calculate_distances/linux_bin/",
    "pear_ebi/calculate_distances/mac_bin/",
    "site/",
    ".venv/",
    ".venv.stale-pre-poetry/",
    "dist/",
    "graphify-out/",
    "__pycache__/",
    ".ipynb_checkpoints/",
    "pear_ebi.egg-info/",
    "docs_assets/",
    "logos/",
)

out = Path("graphify-out")
out.mkdir(exist_ok=True)
result = detect(Path("."))
root = result["scan_root"]


def keep(f):
    rel = os.path.relpath(f, root)
    return not any(rel.startswith(e) for e in EXCLUDE)


code = [f for f in result.get("files", {}).get("code", []) if keep(f)]
result["files"] = {"code": code, "document": [], "paper": [], "image": [], "video": []}
result["total_files"] = len(code)
(out / ".graphify_detect.json").write_text(
    json.dumps(result, ensure_ascii=False), encoding="utf-8"
)
print(f"corpus: {len(code)} first-party code files")

files = []
for f in code:
    files.extend(collect_files(Path(f)) if Path(f).is_dir() else [Path(f)])
ast = extract(files, cache_root=Path("."))
(out / ".graphify_ast.json").write_text(
    json.dumps(ast, indent=2, ensure_ascii=False), encoding="utf-8"
)
(out / ".graphify_semantic.json").write_text(
    json.dumps(
        {
            "nodes": [],
            "edges": [],
            "hyperedges": [],
            "input_tokens": 0,
            "output_tokens": 0,
        }
    ),
    encoding="utf-8",
)
merged = {
    "nodes": ast["nodes"],
    "edges": ast["edges"],
    "hyperedges": [],
    "input_tokens": 0,
    "output_tokens": 0,
}
(out / ".graphify_extract.json").write_text(
    json.dumps(merged, indent=2, ensure_ascii=False), encoding="utf-8"
)

G = build_from_json(merged, root=".", directed=False)
if G.number_of_nodes() == 0:
    sys.exit("ERROR: empty graph")
communities = cluster(G)
cohesion = score_all(G, communities)
gods = god_nodes(G)
surprises = surprising_connections(G, communities)
labels = {cid: f"Community {cid}" for cid in communities}
questions = suggest_questions(G, communities, labels)
# graph.json is removed by the caller so the #479 shrink-guard does not block us
to_json(G, communities, "graphify-out/graph.json")
(out / "GRAPH_REPORT.md").write_text(
    generate(
        G,
        communities,
        cohesion,
        labels,
        gods,
        surprises,
        result,
        {"input": 0, "output": 0},
        ".",
        suggested_questions=questions,
    ),
    encoding="utf-8",
)
(out / ".graphify_analysis.json").write_text(
    json.dumps(
        {
            "communities": {str(k): v for k, v in communities.items()},
            "cohesion": {str(k): v for k, v in cohesion.items()},
            "gods": gods,
            "surprises": surprises,
            "questions": questions,
        },
        indent=2,
        ensure_ascii=False,
    ),
    encoding="utf-8",
)
print(
    f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, {len(communities)} communities"
)
