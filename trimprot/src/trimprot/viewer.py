"""Minimal self-contained viewer: the trimmed target with hotspot residues highlighted.

One 3Dmol.js panel, the trimmed structure INLINED as a JS string (works on file://
with no sibling fetches). Everything is grey cartoon except the hotspot residues,
which are drawn green (cartoon + sticks). Nothing else is annotated.
"""
from __future__ import annotations

import json

VIEWER_JS = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js"
HOTSPOT_COLOR = "0x33cc33"


def build_viewer_html(*, target: str, pdb_id: str, pdb_text: str,
                      hotspot_residues: list[int]) -> str:
    payload = {
        "pdb": pdb_text,
        "hotspots": sorted({int(n) for n in hotspot_residues}),
    }
    data_json = json.dumps(payload)
    n_hot = len(payload["hotspots"])

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{target} ({pdb_id.upper()}) — hotspots</title>
<script src="{VIEWER_JS}"></script>
<style>
  html, body {{ margin:0; height:100%; background:#101014; color:#eee;
                font-family:-apple-system,Segoe UI,Roboto,sans-serif; }}
  #bar {{ padding:8px 14px; font-size:14px; background:#1d1d22; }}
  #bar b {{ color:#7CFC7C; }}
  #viewer {{ position:absolute; top:38px; left:0; right:0; bottom:0; }}
</style>
</head>
<body>
<div id="bar">{target} · PDB {pdb_id.upper()} · trimmed target —
  <b>{n_hot} hotspot residues highlighted in green</b></div>
<div id="viewer"></div>
<script>
const DATA = {data_json};
window.addEventListener('load', function () {{
  const viewer = $3Dmol.createViewer(document.getElementById('viewer'),
                                     {{ backgroundColor: '0x101014' }});
  viewer.addModel(DATA.pdb, 'pdb');
  // base: everything grey cartoon
  viewer.setStyle({{}}, {{ cartoon: {{ color: '0x8899aa' }} }});
  // hotspots: green cartoon + sticks (highlighted on every kept chain)
  if (DATA.hotspots.length) {{
    viewer.setStyle({{ resi: DATA.hotspots }},
                    {{ cartoon: {{ color: '{HOTSPOT_COLOR}' }},
                      stick: {{ color: '{HOTSPOT_COLOR}' }} }});
  }}
  viewer.zoomTo();
  viewer.render();
}});
</script>
</body>
</html>
"""
