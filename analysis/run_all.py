"""Run the whole analysis in order (about 15 min): PK fit, PD scenarios, alternatives, parameter intervals, figures."""
import subprocess, sys, pathlib
here = pathlib.Path(__file__).resolve().parent
for s in ('pk_stage.py', 'pd_stage.py', 'pd_stage3.py', 'pd_stage4.py', 'fig_schematic.py', 'figs_split.py'):
    print('running', s, flush=True); subprocess.run([sys.executable, str(here / s)], check=True)
