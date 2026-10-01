import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import yaml

ROOT=Path(__file__).resolve().parents[1]
SKILL=ROOT/'skills/tsymperturb'

def test_skill_metadata_and_reference_links():
    text=(SKILL/'SKILL.md').read_text()
    front=yaml.safe_load(text.split('---',2)[1])
    assert front['name']=='tsymperturb'
    assert 'longitudinal' in front['description']
    for target in re.findall(r'\]\(([^)]+)\)',text):
        if '://' not in target:
            assert (SKILL/target).is_file(),target

def test_isolated_skill_cli(tmp_path):
    installed=tmp_path/'tsymperturb'
    shutil.copytree(SKILL,installed)
    out=tmp_path/'result.json'
    proc=subprocess.run([sys.executable,str(installed/'scripts/tsymperturb.py'),
                         str(ROOT/'examples/synthetic_model.json'),'--output',str(out)],
                         cwd=tmp_path,capture_output=True,text=True)
    assert proc.returncode==0,proc.stderr
    result=json.loads(out.read_text())
    assert len(result['tvpps'])==4 and len(result['dose_results'])==5

def test_cli_failure_does_not_create_output(tmp_path):
    inp=tmp_path/'bad.json'; inp.write_text('{}')
    out=tmp_path/'not-created.json'
    proc=subprocess.run([sys.executable,str(SKILL/'scripts/tsymperturb.py'),str(inp),
                         '--output',str(out)],capture_output=True,text=True)
    assert proc.returncode==2 and not out.exists()
