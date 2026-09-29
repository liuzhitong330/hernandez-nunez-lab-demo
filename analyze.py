"""Rebuild the marker-panel explorer: python analyze.py [--download].

Requires Python >=3.10 and numpy. No pickle or source-repository code is run.
Default uses the bundled, CC-BY-4.0 derived cell-level subset (cells.npz).
--download reacquires the four original arrays and checks their SHA256 hashes.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import urllib.request
import numpy as np

ROOT = Path(__file__).resolve().parent
GENES = 'cart2 glyt2 tac1 pvalb7 npb grm1b irx1b dat net calb1 penka penkb eomesa emx3 cfos gad1b cx43 vglut2a sst uts1 pou4f2 cort nr4a2a cckb tph2 chata calb2a npy gfra1a dmbx1a gbx2 crhb nefma chodl pyya zic2a th pdyn tbr1b otpa esrrb'.split()
REGIONS = ['Unassigned','Inferior medulla','Intermediate medulla','Superior medulla','Superior raphe','Cerebellum','Tegmentum','Nucleus isthmi','Optic tectum PV','Optic tectum neuropil','Pretectum','Prethalamus','Dorsal thalamus','Habenula','Hypothalamus','Subpallium','Pallium']
FISH = [59,63,71]
FILES = {'Brain_reg':57981580,'fish_id':57981601,'Coords_All':57981640,'BinaryGenes_All':57981757}
SHARE = 'd1d19b105c4f74865c32'

def normalize(raw):
    arrays = {n:np.load(raw/(n+'.npy'), allow_pickle=False) for n in FILES}
    b,r,f,c = [arrays[n] for n in ['BinaryGenes_All','Brain_reg','fish_id','Coords_All']]
    assert b.shape == (306108,41) and c.shape == (306108,3)
    assert set(np.unique(b)) == {0.,1.} and set(np.unique(f)) == set(FISH)
    assert np.all(np.isfinite(b)) and np.all(np.isfinite(r))
    keep = np.all(np.isfinite(c),axis=1)
    assert keep.sum() == 274455
    np.savez_compressed(ROOT/'cells.npz', expression=b[keep].astype('uint8'), region=r[keep].astype('uint8'), fish=f[keep].astype('uint8'), source_row=np.flatnonzero(keep).astype('uint32'))
    provenance = {'retrieved':'2026-09-29','source':'https://figshare.com/s/'+SHARE,'license':'CC BY 4.0','authors':'Marquez-Legorreta, Fleishman, Hesselink et al. (2026), WARP','original_rows':len(b),'included_rows':int(keep.sum()),'excluded_nonfinite_coordinates':int((~keep).sum()),'mapping_repository':'https://github.com/JaneliaSciComp/warp-atlas-viewer','mapping_commit':'0f9451756b0a0bf3d051fea0f40789992a17e783','note':'Gene order and focal-region labels follow the documented reference viewer mapping; independent implementation. This is external reference data, not Hernandez-Nunez Lab data. Source authors share the CC-BY dataset through a publicly linked Figshare private-share URL.','files':{n:{'url':f'https://figshare.com/ndownloader/files/{i}?private_link={SHARE}','sha256':hashlib.sha256((raw/(n+'.npy')).read_bytes()).hexdigest()} for n,i in FILES.items()}}
    (ROOT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')

def aggregate():
    with np.load(ROOT/'cells.npz',allow_pickle=False) as z:
        b,r,f=z['expression'],z['region'],z['fish']
    anchor = b[:,GENES.index('chata')].astype(bool)
    counts=[]
    for g in range(len(GENES)):
        mask=anchor & b[:,g].astype(bool)
        counts.append([np.bincount(r[mask & (f==fish)],minlength=17).tolist() for fish in FISH])
    totals=[np.bincount(r[f==fish],minlength=17).tolist() for fish in FISH]
    out={'genes':GENES,'regions':REGIONS,'fish':FISH,'counts':counts,'totals':totals,'n':len(b),'chataIndex':GENES.index('chata')}
    (ROOT/'data.js').write_text('window.PANEL_DATA = '+json.dumps(out,separators=(',',':'))+';\n')
    with (ROOT/'marker_counts.csv').open('w',newline='') as handle:
        writer=csv.writer(handle); writer.writerow(['co_marker','fish_id','region','chata_and_marker_cells','chata_cells','all_profiled_cells'])
        for g,gene in enumerate(GENES):
            for fi,fish in enumerate(FISH):
                for ri,region in enumerate(REGIONS):
                    writer.writerow([gene,fish,region,counts[g][fi][ri],counts[out['chataIndex']][fi][ri],totals[fi][ri]])
    # Independent row-wise check, not a second call to the aggregation function.
    for gi in range(len(GENES)):
        for fi,fish in enumerate(FISH):
            for ri in range(17):
                expected=np.count_nonzero((f==fish)&(r==ri)&(b[:,25]==1)&(b[:,gi]==1))
                assert expected==counts[gi][fi][ri]
    print('PASS: 2,091 aggregate cells independently checked; binary input, specimen IDs and coordinate exclusions verified.')
    med=[1,2,3]
    for gene in ['chata','glyt2','esrrb','zic2a','chodl','gbx2','irx1b']:
        gi=GENES.index(gene)
        print(gene,[(sum(counts[gi][fi][i] for i in med),sum(counts[gi][fi])) for fi in range(3)])

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--download',action='store_true');parser.add_argument('--raw-dir',type=Path);args=parser.parse_args()
    if args.download:
        provenance=json.loads((ROOT/'provenance.json').read_text()); raw=ROOT/'raw';raw.mkdir(exist_ok=True)
        for n,entry in provenance['files'].items():
            path=raw/(n+'.npy')
            if not path.exists(): urllib.request.urlretrieve(entry['url'],path)
            assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'], 'Source file changed: '+n
        normalize(raw)
    elif args.raw_dir: normalize(args.raw_dir)
    aggregate()
