"""Remove plot clipping groups; each line is a standalone editable SVG object."""
from copy import deepcopy
from pathlib import Path
import json
import re
import xml.etree.ElementTree as ET

REPO=Path(__file__).resolve().parents[1]
FOLDER=Path('/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples')
SOURCE=FOLDER/'csChrimson_all84_onset_aligned_rotation_mean_SEM_red_gray.svg'
DEST=FOLDER/'csChrimson_all84_onset_aligned_rotation_independent_lines.svg'
NS='http://www.w3.org/2000/svg'
ET.register_namespace('',NS)

def styles(text):
    return {key.strip():value.strip() for part in text.split(';') if ':' in part for key,value in [part.split(':',1)]}

def main():
    original=ET.parse(SOURCE).getroot()
    root=ET.Element(original.tag,{k:v for k,v in original.attrib.items() if k!='class'})
    assert all(set(g.attrib)<= {'clip-path'} for g in original.iter('{'+NS+'}g'))
    records=json.loads((REPO/'outputs/chrimson_onset_traces_20261103/extraction_results.json').read_text())
    order=['ves041_1_Chrim','ves041_2_Chrim','ves041_1_gfp','ves041_2_gfp','empty_Chrim']
    ids={g:sorted(r['fly'] for r in records if r['genotype']==g) for g in order}
    clips={}
    for c in original.iter('{'+NS+'}clipPath'):
        rect=c.find('{'+NS+'}rect')
        clips[c.get('id')]=tuple(float(rect.get(a)) for a in ['x','y','width','height'])
    panels=sorted({v for v in clips.values() if 200<v[2]<300 and v[3]>150},key=lambda v:v[0])
    assert len(panels)==5
    panel_to_genotype=dict(zip(panels,order));counter={g:0 for g in order}
    nfly=nmean=nband=0
    defaults=dict(fill='none',stroke='#000000',**{'stroke-linecap':'round','stroke-linejoin':'round','stroke-miterlimit':'10.00'})
    def walk(parent,clip=None):
        nonlocal nfly,nmean,nband
        for e in parent:
            tag=e.tag.split('}')[-1]
            if tag=='defs':continue
            if tag=='g':
                ref=re.fullmatch(r'url\(#([^)]*)\)',e.get('clip-path',''))
                walk(e,clips[ref.group(1)] if ref else clip);continue
            node=deepcopy(e)
            st=styles(node.get('style',''))
            if tag in ['line','polyline','polygon','path','rect','circle']:
                st={**defaults,**st}
            elif tag=='text':st={'white-space':'pre',**st}
            node.set('style','; '.join(f'{k}: {v}' for k,v in st.items())+';')
            genotype=panel_to_genotype.get(clip)
            identifier=f'{tag}_{len(root)+1:03d}'
            points=node.get('points','').split()
            if tag=='polyline' and len(points)==300 and st.get('stroke-opacity')=='0.28':
                assert genotype is not None
                index=counter[genotype];fly=ids[genotype][index];counter[genotype]+=1
                identifier=f'{genotype}_fly{fly:03d}_rotation_trace';nfly+=1
            elif tag=='polyline' and len(points)==300 and st.get('stroke-width')=='2.45':
                identifier=f'{genotype}_genotype_mean';nmean+=1
            elif tag=='polygon' and st.get('fill-opacity')=='0.23':
                identifier=f'{genotype}_SEM_band';nband+=1
            elif tag=='text':
                identifier='text_'+str(len(root)+1)+'_'+re.sub('[^a-zA-Z0-9]+','_', ''.join(node.itertext()))[:55]
            node.set('id',identifier)
            node.attrib.pop('clip-path',None)
            root.append(node)
    walk(original)
    assert (nfly,nmean,nband)==(84,5,5)
    assert all(counter[g]==len(ids[g]) for g in order)
    assert not list(root.iter('{'+NS+'}g')) and not list(root.iter('{'+NS+'}clipPath'))
    assert not any('clip-path' in e.attrib for e in root.iter())
    # Every data-trace coordinate and its original per-object visual settings
    # remain intact. Only shared structural groups/masks/CSS are removed.
    old=[e.get('points') for e in original.iter('{'+NS+'}polyline') if len(e.get('points','').split())==300]
    new=[e.get('points') for e in root.iter('{'+NS+'}polyline') if len(e.get('points','').split())==300]
    assert old==new
    assert not DEST.exists(),f'Refusing to overwrite: {DEST}'
    ET.ElementTree(root).write(DEST,encoding='utf-8',xml_declaration=True)
    print(f'Saved {DEST}\n84 independent fly traces, 5 independent mean lines, 5 independent SEM bands. No groups, clipping masks, raster images or linked instances.')

if __name__=='__main__':main()
