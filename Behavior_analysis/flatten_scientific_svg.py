"""Make svglite scientific figure objects independent, preserving coordinates."""
from copy import deepcopy
from pathlib import Path
import sys
import xml.etree.ElementTree as E

def main(source,destination):
    ns='http://www.w3.org/2000/svg';E.register_namespace('',ns)
    old=E.parse(source).getroot();root=E.Element(old.tag,{k:v for k,v in old.attrib.items() if k!='class'})
    assert all(set(g.attrib)<={'clip-path'} for g in old.iter('{'+ns+'}g'))
    defaults={'fill':'none','stroke':'#000000','stroke-linecap':'round','stroke-linejoin':'round','stroke-miterlimit':'10.00'}
    def walk(parent):
        for item in parent:
            tag=item.tag.split('}')[-1]
            if tag=='defs':continue
            if tag=='g':walk(item);continue
            node=deepcopy(item);style={k.strip():v.strip() for s in node.get('style','').split(';') if ':' in s for k,v in [s.split(':',1)]}
            if tag in ['line','polyline','polygon','path','rect','circle']:style={**defaults,**style}
            elif tag=='text':style={'white-space':'pre',**style}
            node.set('style','; '.join(f'{k}: {v}' for k,v in style.items())+';')
            node.set('id',f'{tag}_{len(root)+1:04d}');node.attrib.pop('clip-path',None);root.append(node)
    walk(old)
    assert not any(e.tag.split('}')[-1] in ['g','clipPath','image','use'] for e in root.iter())
    for attr in ['points','d','x1','x2','y1','y2']:
        assert [e.get(attr) for e in old.iter() if e.get(attr) is not None]==[e.get(attr) for e in root.iter() if e.get(attr) is not None]
    assert not Path(destination).exists()
    E.ElementTree(root).write(destination,encoding='utf-8',xml_declaration=True)
    print('Saved ungrouped SVG:',destination)

if __name__=='__main__':main(*sys.argv[1:])
