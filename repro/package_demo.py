"""Create small public summaries and on-demand trace shards from checked outputs."""
import argparse, gzip, hashlib, json, shutil
from pathlib import Path
import numpy as np

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();source=Path(__file__).parent;out=Path(args.output)
    (out/'derived').mkdir(exist_ok=True);(out/'repro').mkdir(exist_ok=True)
    records=[];provenance=[]
    configs=[('fish1-analysis.json','blood_glucose_fish1',7,1,'Blood D-glucose · fish 1'),
             ('fish2-analysis.json','blood_glucose_fish2',4,1,'Blood D-glucose · fish 2'),
             ('gut3-analysis.json','glucose_fish3',2,3,'Gut D-glucose · fish 3')]
    for result_name,recording,target,reference,label in configs:
        path=source/result_name
        if not path.exists():continue
        result=json.loads(path.read_text());manifest=json.loads((source/(recording+'_manifest.json')).read_text())
        x=np.load(source/(recording+'_subset.npz'))
        result.update(label=label,target_position=target,reference_position=reference,
                      visual_available=bool(np.ptp(x['visual'])>0),
                      source_n_cells=manifest['source_n_cells'],
                      mapping_status='Unverified position-to-organ mapping. Codes are literal acquisition metadata, not biological labels.')
        for start in range(0,len(result['cells']),32):
            shard={};name='derived/'+recording+'-traces-'+str(start//32).zfill(2)+'.json.gz'
            for cell in result['cells'][start:start+32]:
                shard[str(cell['id'])]={'trace':cell.pop('trace'),'predictions':cell.pop('predictions')}
                cell['shard']=name
            (out/name).write_bytes(gzip.compress(json.dumps(shard,separators=(',',':'),allow_nan=False).encode(),mtime=0))
        summary='derived/'+recording+'.json'
        (out/summary).write_text(json.dumps(result,separators=(',',':'),allow_nan=False))
        records.append(summary)
        keys=['time','traces','coords','cell_ids','baseline_min','motor','visual','gut_event_times',
              'nutrient_event_times','control_event_times','other_event_times','all_event_times','all_event_positions','selected_target_position']
        minimal=out/'repro'/(recording+'_subset.npz')
        np.savez_compressed(minimal,**{k:x[k] for k in keys if k in x.files})
        # Preserve the full source/transform manifest but no local full-data paths.
        manifest['public_minimal_subset_sha256']=hashlib.sha256(minimal.read_bytes()).hexdigest()
        manifest['selected_position_code']=target;manifest['reference_position_code']=reference
        (out/'repro'/(recording+'_manifest.json')).write_text(json.dumps(manifest,indent=2))
        provenance.append({'recording':recording,'label':label,'sampled_cells':len(result['cells']),
                           'source_n_cells':manifest['source_n_cells'],'manifest':'repro/'+recording+'_manifest.json',
                           'minimal_subset':'repro/'+recording+'_subset.npz','sha256':manifest['public_minimal_subset_sha256'],
                           'position_mapping':'unverified','visual_modulation_present':result['visual_available']})
    # The gut record is the default, because visual modulation is actually present.
    records.sort(key=lambda s:(0 if '/glucose_fish3.' in s else 1,s))
    (out/'data.js').write_text('window.DEMO_RECORDS = '+json.dumps(records)+';\n')
    (out/'provenance.json').write_text(json.dumps({'version':2,'date':'2026-09-29',
      'dataset':'https://doi.org/10.25378/janelia.31967640.v1','data_license':'CC BY 4.0',
      'publication':'https://www.nature.com/articles/s41467-026-76242-8',
      'code':'https://github.com/benJames1212/gutBrainPipeline',
      'commit':'4c62612ecd540ec786c43fe9ce67cc78259a6ec1','source_code_license':'MIT',
      'records':provenance,'ownership':'External Janelia dataset; not Hernandez-Nunez Lab data',
      'limitations':['Stimulus-code organ mapping is not deposited; no nutrient-specific biological classification.',
      'Discovery ranks are frozen before testing; interactive inspection does not constitute prospective validation.',
      'Native voxel-space ROI centroids; z display expanded 4x, no registered brain atlas or inferred connections.']},indent=2))
    for name in ['model_audit.py','test_model_audit.py','run_model_audit.py','extract_subset.py','range_access.py','figshare-31967640.json']:
        if (source/name).exists():shutil.copy2(source/name,out/'repro'/name)
    shutil.copy2(source.parent/'gutBrainPipeline'/'LICENSE',out/'repro'/'UPSTREAM_LICENSE')
    shutil.copy2(Path(__file__),out/'repro'/'package_demo.py')
    print(json.dumps({'records':records,'files':len(list((out/'derived').iterdir())),'data_bytes':sum(p.stat().st_size for p in (out/'derived').iterdir())}))

if __name__=='__main__':main()
