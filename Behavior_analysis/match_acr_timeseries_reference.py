"""Match the user's edited control SVG coordinate scales exactly."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/ves041_original_acr_pre_light_example_20261104'
REFERENCE=Path('/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples/Control_inhibition_acr_timeseries.svg')
NS='http://www.w3.org/2000/svg'

def main():
 m=json.loads((OUT/'selected_example.json').read_text())
 raw=np.genfromtxt(OUT/'ACR_Fly22C_pulse31_raw_velocity.csv',delimiter=',',names=True)
 dt=float(np.median(np.diff(raw['time_s'])));sigma=.010;radius=int(np.ceil(4*sigma/dt));offset=np.arange(-radius,radius+1)*dt
 kernel=np.exp(-.5*(offset/sigma)**2);kernel/=kernel.sum()
 forward=np.convolve(np.pad(raw['forward_mm_s'],radius,mode='edge'),kernel,mode='valid')
 yaw=np.convolve(np.pad(raw['yaw_deg_s'],radius,mode='edge'),kernel,mode='valid')
 lo=m['window_zoom'][0]-1;hi=lo+10;sel=(raw['time_s']>=lo)&(raw['time_s']<=hi)
 tree=ET.parse(REFERENCE);root=tree.getroot();ET.register_namespace('',NS)
 # These scale constants are read from the supplied reference's axis ticks
 # and 1 s scale bar, rather than from the earlier generated plot.
 xleft=34.063;xright=329.514;xscale=(xright-xleft)/10
 fy0=99.775;fscale=(139.525-99.775)/10
 ry0=306.4;rscale=(306.4-226.9)/200
 x=xleft+(raw['time_s'][sel]-lo)*xscale
 f_y=fy0-forward[sel]*fscale;r_y=ry0-yaw[sel]*rscale
 def path(xx,yy):return 'M'+' L'.join(f'{a:.6f},{b:.6f}' for a,b in zip(xx,yy))
 found=[];shading=[]
 for e in root.iter('{'+NS+'}path'):
  style=e.get('style','')
  if 'stroke:#000000' in style and 'stroke-width:1.1' in style:
   e.set('d',path(x,f_y));e.set('id','acr_forward_velocity');found.append('forward')
  elif 'stroke:#004DD9' in style and 'stroke-width:1.1' in style:
   e.set('d',path(x,r_y));e.set('id','acr_rotation_velocity');found.append('rotation')
  elif 'fill:#F0E1FA' in style:shading.append(e)
 assert sorted(found)==['forward','rotation'] and len(shading)==2
 light_x=xleft+(m['light_start_s']-lo)*xscale
 width=(m['light_end_s']-m['light_start_s'])*xscale
 for e,y in zip(shading,[.4,226.9]):e.set('d',f'M{light_x:.6f},{y:.6f}h{width:.6f}v159h{-width:.6f}Z')
 event=(raw['time_s']>=m['start_s'])&(raw['time_s']<=m['end_s'])
 center_time=(m['start_s']+m['end_s'])/2;center_yaw=(min(yaw[event])+max(yaw[event]))/2
 radius_t=(m['end_s']-m['start_s'])/2+.18;radius_yaw=(max(yaw[event])-min(yaw[event]))/2+35
 ET.SubElement(root,'{'+NS+'}ellipse',dict(id='original_detected_event',cx=f'{xleft+(center_time-lo)*xscale:.6f}',cy=f'{ry0-center_yaw*rscale:.6f}',rx=f'{radius_t*xscale:.6f}',ry=f'{radius_yaw*rscale:.6f}',fill='none',stroke='#D9331A',**{'stroke-width':'1.1'}))
 # The reference has no clipping regions; retaining complete curves preserves
 # true peaks beyond the labeled rotation axis range at the same gain.
 assert min(f_y)>=.4 and max(f_y)<=159.4
 assert min(r_y)>159.4 and max(r_y)<420.403
 metadata=dict(reference=str(REFERENCE),recording_window_s=[lo,hi],relative_window_s=[lo-m['light_start_s'],hi-m['light_start_s']],window_shift_s=-1,gaussian_sigma_s=sigma,
  viewBox=root.get('viewBox'),plot_width_svg_units=xright-xleft,time_svg_units_per_s=xscale,forward_svg_units_per_mm_s=fscale,rotation_svg_units_per_deg_s=rscale,
  light_on_x=light_x,light_on_width=width,clipping='None; full traces preserved at reference gain')
 ET.SubElement(root,'{'+NS+'}metadata',id='analysis_provenance').text=json.dumps(metadata)
 dest=OUT/'ACR_Fly22C_circled_shift1s_reference_scale.svg'
 tree.write(dest,encoding='utf-8',xml_declaration=True)
 (OUT/'ACR_Fly22C_circled_shift1s_reference_scale.json').write_text(json.dumps(metadata,indent=2))
 print(json.dumps(metadata,indent=2))

if __name__=='__main__':main()
