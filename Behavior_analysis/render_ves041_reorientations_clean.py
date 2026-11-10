"""Presentation SVGs with editable text and native trajectory coordinates."""
from pathlib import Path
from html import escape
import csv
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'outputs/ves041_example_reorientation_timeseries_20261104'
OUT=ROOT/'outputs/ves041_reorientations_clean_20261104'
PURPLE='#963DCE'; SHADE='#F0E1FA'; ORANGE='#C55216'; INK='#202124'; GRAY='#686B70'

class SVG:
 def __init__(self,w,h):
  self.w=w;self.h=h;self.parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">', '<rect width="100%" height="100%" fill="white"/>', '<g font-family="Arial, Helvetica, sans-serif" fill="#202124">'];self.nclip=0
 def text(self,x,y,s,size=17,color=INK,anchor='start',weight='normal'):
  self.parts.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" fill="{color}" text-anchor="{anchor}" font-weight="{weight}">{escape(str(s))}</text>')
 def line(self,x1,y1,x2,y2,color=INK,width=1,dash=None):
  self.parts.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{width}"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>')
 def rect(self,x,y,w,h,color):self.parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{color}"/>')
 def poly(self,x,y,color=INK,width=1.2):
  self.parts.append('<polyline points="'+' '.join(f'{a:.2f},{b:.2f}' for a,b in zip(x,y))+f'" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>')
 def circle(self,x,y,r,color,stroke='white'):
  self.parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" fill="{color}" stroke="{stroke}" stroke-width="1.2"/>')
 def clip(self,x,y,w,h):
  self.nclip+=1;name=f'clip{self.nclip}'
  self.parts.append(f'<defs><clipPath id="{name}"><rect x="{x}" y="{y}" width="{w}" height="{h}"/></clipPath></defs><g clip-path="url(#{name})">')
 def endclip(self):self.parts.append('</g>')
 def save(self,p):p.parent.mkdir(parents=True,exist_ok=True);p.write_text('\n'.join(self.parts+['</g></svg>']))

def rows(p):
 with p.open() as f:return [{k:(v if k=='phase' else float(v)) for k,v in r.items()} for r in csv.DictReader(f)]

def load(id):
 return dict(id=id,label='VES041_1 × '+('ACR' if id==8 else 'GFP'),d=np.genfromtxt(SRC/f'fly{id}_timeseries.csv',delimiter=',',names=True),e=rows(SRC/f'fly{id}_events.csv'),light=rows(SRC/f'fly{id}_light_intervals.csv'))

def legend(s,y=91):
 s.rect(80,y-13,22,16,SHADE);s.text(112,y,'Light on',17)
 s.line(230,y-6,262,y-6,ORANGE,3);s.text(273,y,'Reorientation',17)
 s.line(465,y-6,497,y-6,INK,1.5);s.text(508,y,'Recorded behavior',17)

def title(s,r,sub):
 s.text(80,43,f'Fly {r["id"]}  |  {r["label"]}',28,weight='bold')
 s.text(80,70,sub,17,color=GRAY)

def bounds_labels(events,lo,hi,x,w):
 positions=[];lanes=[-1e9]*3
 for e in events:
  if not lo<=e['marker_s']<hi:continue
  px=x+(e['marker_s']-lo)/(hi-lo)*w
  lane=next((i for i,v in enumerate(lanes) if px-v>34),int(np.argmin(lanes)))
  lanes[lane]=px;positions.append((e,px,lane))
 return positions

def trace(s,r,lo,hi,col,box,yl,ticks,label,xticks=False,event_labels=False):
 x,y,w,h=box;d=r['d'];take=(d['time_s']>=lo)&(d['time_s']<=hi)
 tx=lambda t:x+(np.asarray(t)-lo)/(hi-lo)*w
 ty=lambda v:y+h-(np.asarray(v)-yl[0])/(yl[1]-yl[0])*h
 s.clip(x,y,w,h)
 for a in r['light']:
  start=max(lo,a['start_s']);end=min(hi,a['end_s'])
  if end>start:s.rect(float(tx(start)),y,float(tx(end)-tx(start)),h,SHADE)
 if yl[0]<0<yl[1]:s.line(x,float(ty(0)),x+w,float(ty(0)),'#DADCE0',.8)
 s.poly(tx(d['time_s'][take]),ty(d[col][take]),INK,1.15)
 for e in r['e']:
  if e['end_s']<lo or e['start_s']>hi:continue
  sel=take&(d['time_s']>=e['start_s'])&(d['time_s']<=e['end_s'])
  s.poly(tx(d['time_s'][sel]),ty(d[col][sel]),ORANGE,2.4)
  s.line(float(tx(e['marker_s'])),y,float(tx(e['marker_s'])),y+h,ORANGE,.75,'3 5')
 s.endclip()
 s.line(x,y,x,y+h,'#888B90',.85)
 for v in ticks:
  py=float(ty(v));s.line(x-5,py,x,py,'#888B90',.8);s.text(x-12,py+5,f'{v:g}',15,GRAY,'end')
 # Horizontal labels are clearer than long rotated axis titles.
 s.text(x,y-13,label,17,weight='bold')
 if xticks:
  s.line(x,y+h,x+w,y+h,'#888B90',.85)
  step=50 if hi-lo>100 else 10 if hi-lo>30 else 5 if hi-lo>10 else 1
  for t in np.arange(np.ceil(lo/step)*step,hi+1e-6,step):
   px=float(tx(t));s.line(px,y+h,px,y+h+5,'#888B90',.8);s.text(px,y+h+25,f'{t:g}',15,GRAY,'middle')
 if event_labels:
  for e,px,lane in bounds_labels(r['e'],lo,hi,x,w):
   s.text(px,y-13-20*lane,f'#{int(e["event_id"])}',15,ORANGE,'middle','bold')

def sheet(r,start):
 s=SVG(1280,1200);title(s,r,f'{start}–{start+600} s  ·  Full recording, continued in three 200-second sections');legend(s,103)
 for i in range(3):
  lo=start+200*i;hi=lo+200;y=166+i*320
  s.text(80,y,f'{lo}–{hi} s',20,weight='bold')
  trace(s,r,lo,hi,'forward_mm_s',(110,y+37,1100,88),(-12,25),[-10,0,10,20],'Forward velocity (mm/s)',event_labels=False)
  trace(s,r,lo,hi,'yaw_deg_s',(110,y+172,1100,91),(-1000,1000),[-1000,0,1000],'Yaw velocity (deg/s)',xticks=True,event_labels=True)
 s.text(1210,1160,'Time from recording start (s)',17,GRAY,'end')
 s.text(80,1160,'Event IDs match the flat path and close-ups.',17,GRAY)
 s.save(OUT/f'Fly{r["id"]:02d}_timeseries_{start:04d}_{start+600:04d}s.svg')

def path(s,r,box,lo=0,hi=1200,labels=True,callouts=False):
 x,y,w,h=box;d=r['d'];sel=(d['time_s']>=lo)&(d['time_s']<=hi)
 xx=d['x_mm'][sel];yy=d['y_mm'][sel]
 spanx=np.ptp(xx);spany=np.ptp(yy);scale=min(w/(spanx+max(2,spanx*.10)),h/(spany+max(2,spany*.10)))
 cx=(max(xx)+min(xx))/2;cy=(max(yy)+min(yy))/2
 px=lambda a:x+w/2+(np.asarray(a)-cx)*scale
 py=lambda a:y+h/2-(np.asarray(a)-cy)*scale
 s.poly(px(xx),py(yy),INK,1.3)
 light=sel&(d['light_on']==1)
 inds=np.flatnonzero(light)
 for run in np.split(inds,np.flatnonzero(np.diff(inds)>1)+1):
  if len(run)>1:s.poly(px(d['x_mm'][run]),py(d['y_mm'][run]),PURPLE,1.8)
 events=[e for e in r['e'] if lo<=e['marker_s']<hi]
 for e in events:
  use=(d['time_s']>=e['start_s'])&(d['time_s']<=e['end_s'])
  s.poly(px(d['x_mm'][use]),py(d['y_mm'][use]),ORANGE,3.5)
  s.circle(float(px(e['x_mm'])),float(py(e['y_mm'])),4,ORANGE)
 s.circle(float(px(xx[0])),float(py(yy[0])),5,'#31A354')
 # Scale bar retains equal x/y geometry without a large empty axis box.
 bar=100 if max(spanx,spany)>300 else 50 if max(spanx,spany)>100 else 1
 s.line(x+15,y+h-10,x+15+bar*scale,y+h-10,INK,2.5)
 s.text(x+15,y+h-22,f'{bar} mm',17)
 if labels and callouts and events:
  placed=[]
  points=[(float(px(e['x_mm'])),float(py(e['y_mm']))) for e in events]
  for e,(ex,ey) in zip(events,points):
   label=f'#{int(e["event_id"])}';tw=9*len(label)+8;th=22
   best=None;best_cost=float('inf')
   for radius in [17,28,42,60,80,105]:
    for angle in np.arange(16)*np.pi/8:
     lx=ex+radius*np.cos(angle);ly=ey+radius*np.sin(angle)
     box=(lx-tw/2,ly-th/2,lx+tw/2,ly+th/2)
     if box[0]<25 or box[2]>s.w-25 or box[1]<y-10 or box[3]>y+h+15:continue
     overlap=sum(box[0]<b[2] and box[2]>b[0] and box[1]<b[3] and box[3]>b[1] for b in placed)
     cover=sum(box[0]<a<box[2] and box[1]<b<box[3] for a,b in points)
     cost=overlap*10000+cover*100+radius
     if cost<best_cost:best_cost=cost;best=(lx,ly,box)
   if best:
    lx,ly,box=best;placed.append(box)
    s.line(ex,ey,lx,ly,'#CDA78F',.8)
    s.rect(box[0],box[1],tw,th,'white')
    s.text(lx,ly+5,label,15,ORANGE,'middle','bold')

 elif labels:
  for e in events:s.text(float(px(e['x_mm']))+10,float(py(e['y_mm']))-9,f'#{int(e["event_id"])}',18,ORANGE,weight='bold')

def full_path(r):
 s=SVG(1280,860);title(s,r,'Full 20-minute flat path  ·  Event numbers match the time series');legend(s,103)
 path(s,r,(80,150,1120,610),callouts=True)
 s.text(80,815,'Measured x/y path; equal spatial scale in both directions. Green dot: recording start.',17,GRAY)
 s.save(OUT/f'Fly{r["id"]:02d}_full_flatpath.svg')

def closeup(r,e):
 lo=max(0,e['start_s']-3);hi=min(1200,e['end_s']+3)
 s=SVG(1280,800);title(s,r,f'Event #{int(e["event_id"])}  ·  {e["marker_s"]:.2f} s  ·  '+('Pre-light' if e['phase']=='pre' else 'Light on'));legend(s,103)
 s.text(80,150,'Flat path',20,weight='bold');path(s,r,(65,175,380,480),lo,hi)
 s.text(540,150,'Time series',20,weight='bold')
 trace(s,r,lo,hi,'forward_mm_s',(570,198,650,108),(-12,25),[-10,0,10,20],'Forward velocity (mm/s)')
 trace(s,r,lo,hi,'yaw_deg_s',(570,374,650,108),(-1000,1000),[-1000,0,1000],'Yaw velocity (deg/s)',event_labels=True)
 use=(r['d']['time_s']>=lo)&(r['d']['time_s']<=hi)
 hh=r['d']['heading_change_deg'][use];lower=np.floor(min(hh)/30)*30;upper=np.ceil(max(hh)/30)*30
 if upper<=lower:upper=lower+30
 trace(s,r,lo,hi,'heading_change_deg',(570,550,650,108),(lower-10,upper+10),[lower,(lower+upper)/2,upper],'Heading change from recording start (deg)',xticks=True)
 s.text(1220,720,'Time from recording start (s)',17,GRAY,'end')
 s.text(80,758,f'Net turn {e["rotation_deg"]:.1f}°   ·   Net displacement {e["displacement_mm"]:.2f} mm   ·   Directionality {e["directionality"]:.2f}',18)
 s.save(OUT/'events'/f'Fly{r["id"]:02d}_event{int(e["event_id"]):03d}.svg')

def main():
 OUT.mkdir(exist_ok=True)
 for id in [8,45]:
  r=load(id);full_path(r)
  for start in [0,600]:sheet(r,start)
  for e in r['e']:closeup(r,e)
  print(f'Fly {id}: full path, two time-series sheets, {len(r["e"])} close-ups',flush=True)
 (OUT/'READ_ME.txt').write_text('Purple: light on. Orange: existing reorientation detections.\nEach time-series sheet contains three consecutive 200 s sections; two sheets cover the full 20 minutes.\nEvent numbers match the full flat path and individual event close-ups.\nOnly matched 5 s pre-light and light-on windows were scored in the original analysis; the rest of the recording is unscored.\nDetection: >30 degrees net turn in approximately 1 s, <=2 mm net displacement, directionality >=0.7.\nNumbers label the midpoint of the stored representative detection window.\nNative path geometry, existing events, and stored velocity samples are unchanged. Full-path framing is fitted to each recording; scale bars show distance.\nAll labels are editable SVG text.\nData and event verification: ../ves041_example_reorientation_timeseries_20261104/verification.json\n')

if __name__=='__main__':main()
