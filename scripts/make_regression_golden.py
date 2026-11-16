"""Write tests/golden/regression_tile.json (aggregate numbers only) from the
regression tile and the paper's frozen Stage 2 forests. Re-run only when a
change in these numbers is intended, and say why in the CHANGELOG."""
import json
from pathlib import Path

from coralremote import io, paths, regression
from coralremote.forest import load_treebagger

dates, X = regression.load_tile()
ref = io.load_stage2_output()
gate = load_treebagger(paths.STAGE2_OUTPUT, "rfGate")
stress = load_treebagger(paths.STAGE2_OUTPUT, "rfStress")
out = regression.tile_summary(dates, X, gate, stress, (ref["muHealthy"], ref["muBleached"]))
p = Path(__file__).resolve().parent.parent / "tests" / "golden" / "regression_tile.json"
p.write_text(json.dumps({"n_pixels": regression.N_PIXELS, "dates": out}, indent=1))
print("wrote", p)
