"""Identify and replay the archived model, independently of revised primary."""
import os,sys,json,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import numpy as np
from coralremote import io,paths,stage2
from coralremote.forest import load_treebagger
from coralremote.evaluate import compare_date_summary,parse_report_confusion
ap=argparse.ArgumentParser();ap.add_argument('--json',required=True);a=ap.parse_args()
out=paths.assert_writable(a.json)  # refuse writes inside the protected source folders
g=load_treebagger(paths.STAGE2_OUTPUT,'rfGate');s=load_treebagger(paths.STAGE2_OUTPUT,'rfStress')
z=stage2.prepare_frozen_inputs(stage2.Stage2Options(mode='compat'));result=stage2.run_frozen(g,s,z)
n,bad,maxdiff,_=compare_date_summary(result['date_summary'],io.read_date_summary_csv(),tol=1e-9)
exact=np.array_equal(result['metrics']['C'],parse_report_confusion(paths.STAGE2_REPORT.read_text()))
assert bad==0 and exact
out.write_text(json.dumps(dict(kind='legacy replay',fixed_Stage1_pair=[.82,.8],Stage2_mode='compat',original_mask_px=int(z['support'].sum()),date_cells=n,max_difference_pp=maxdiff,confusion_exact=bool(exact),revised_primary=False),indent=2))
print('LEGACY REPLAY: archived fixed coefficients and stored-unit Stage2; revised primary is a different entry point.')
