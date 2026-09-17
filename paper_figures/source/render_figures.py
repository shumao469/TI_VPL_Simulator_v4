#!/usr/bin/env python3
"""Render six editable, data-derived TI/Allen atlas figures without OpenGL.

Run from the repository root: python paper_figures/source/render_figures.py
Every 3D mesh is projected orthographically and shaded deterministically.
"""
from __future__ import annotations
import argparse
from dataclasses import replace
import json
from pathlib import Path
import string
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
from scipy.sparse import coo_matrix
from model import Atlas, DEFAULT, TARGET, ROOT, preview_fields, crossed_montage, generate_data

COLORS = dict(ink='#243b53', gray='#84939f', vpl='#d6a12a', ch1='#287cad', ch2='#d3674b')
CMAP = plt.get_cmap('turbo')
NORM = LogNorm(.01, 1.)
LEVELS = np.geomspace(.01, 1., 42)
STATUS = 'Analytic model preview · generic atlas mapping'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
    'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'axes.linewidth':.65,'lines.linewidth':1.1,'text.color':COLORS['ink'],
    'axes.labelcolor':COLORS['ink'],'xtick.color':COLORS['ink'],'ytick.color':COLORS['ink'],
    'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none',
    'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})


def mesh(payload, key):
    v=np.array(payload[key]['vertices'],float)
    f=np.array(payload[key]['faces'],int)
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    edges=np.unique(np.sort(edges,axis=1),axis=0)
    row=np.r_[edges[:,0],edges[:,1]];col=np.r_[edges[:,1],edges[:,0]]
    adj=coo_matrix((np.ones(len(row)),(row,col)),shape=(len(v),len(v))).tocsr()
    degree=np.maximum(np.asarray(adj.sum(1)).ravel(),1)
    for _ in range(8):
        v=.85*v+.15*(adj@v)/degree[:,None]
    tri=v[f];centers=tri.mean(1)
    normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normals/=np.maximum(np.linalg.norm(normals,axis=1)[:,None],1e-15)
    return dict(v=v,f=f,tri=tri,centers=centers,normals=normals)


def basis(view):
    w=np.array(view,float);w/=np.linalg.norm(w)
    up=np.array([0,0,1.]) if abs(w[2])<.97 else np.array([1.,0,0])
    u=np.cross(up,w);u/=np.linalg.norm(u);v=np.cross(w,u)
    return np.stack([u,v,w],1)


def surface(ax, data, B, select=None, alpha=.20, color=None, zorder=1):
    select=np.ones(len(data['f']),bool) if select is None else select
    triangles=data['tri'][select]@B
    normal=data['normals'][select]
    light=np.array([-.35,-.4,.85]);light/=np.linalg.norm(light)
    shade=.64+.36*np.abs(normal@light)
    base=np.array(matplotlib.colors.to_rgb(color or '#c4cbd1'))
    rgba=np.column_stack([base[None,:]*shade[:,None],np.full(len(normal),alpha)])
    idx=np.argsort(triangles[:,:,2].mean(1))
    coll=PolyCollection(triangles[idx,:,:2],facecolors=rgba[idx],edgecolors='none',
                        antialiased=False,zorder=zorder)
    ax.add_collection(coll)


def silhouette(ax,data,B,color='#73808a',alpha=.45,lw=.35,zorder=5):
    f=data['f'];front=data['normals']@B[:,2]>0
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    edges=np.sort(edges,axis=1)
    unique,inv,counts=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
    nfront=np.bincount(inv,weights=np.tile(front.astype(int),3),minlength=len(unique))
    selected=(counts==1)|((nfront>0)&(nfront<counts))
    lines=(data['v'][unique[selected]]@B)[...,:2]
    ax.add_collection(LineCollection(lines,colors=color,linewidths=lw,alpha=alpha,zorder=zorder))


def mark(ax,point,B,color,shape='o',size=28,zorder=10):
    p=np.asarray(point)@B
    ax.scatter([p[0]],[p[1]],s=size,marker=shape,facecolor=color,edgecolor='white',linewidth=.7,zorder=zorder)
    return p[:2]


def panel_title(ax,letter,title):
    ax.text(-.08,1.035,letter,transform=ax.transAxes,fontsize=11,weight='bold',va='bottom')
    ax.set_title(title,loc='left',pad=10)


def foot(fig,extra=''):
    fig.text(.025,.018,STATUS + ('  |  '+extra if extra else ''),fontsize=6.3,color='#627180')


def colorbar(fig,rect,norm=NORM,label='mTI / baseline brain maximum  (log scale)',ticks=None):
    cax=fig.add_axes(rect)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=CMAP),cax=cax,
                   orientation='horizontal',extend='both',ticks=ticks or [.01,.03,.1,.3,1])
    cb.ax.set_xticklabels([f'{x:g}' for x in (ticks or [.01,.03,.1,.3,1])])
    cb.set_label(label,labelpad=4,fontsize=7.2);cb.outline.set_linewidth(.5)
    if cb.solids is not None: cb.solids.set_rasterized(False)
    return cb


def brain_scene(ax, atlas, meshes, kind='sagittal', params=DEFAULT,
                montage=False, labels=False, zoom=False, field=True):
    view={'sagittal':[.18,-1,.16], 'coronal':[1,-.25,.12], 'axial':[.22,-.24,1]}[kind]
    B=basis(view)
    brain,vpl=meshes
    data=atlas.plane(kind,params=params,nx=220)
    if field:
        k=data['axis'];c=data['coordinate']
        front=(brain['centers'][:,k]-c)*B[k,2]>0
        surface(ax,brain,B,~front,alpha=.25,zorder=1)
        projected=data['xyz']@B
        values=np.ma.masked_where(~data['brain'],data['mTI']/atlas.reference)
        filled=ax.contourf(projected[:,:,0],projected[:,:,1],values,levels=LEVELS,
                    norm=NORM,cmap=CMAP,extend='both',antialiased=False,zorder=2)
        filled.set_edgecolor('face');filled.set_linewidth(.12)
        surface(ax,brain,B,front,alpha=.095,zorder=3)
    else:
        surface(ax,brain,B,alpha=.11,zorder=1)
    silhouette(ax,brain,B,zorder=4)
    surface(ax,vpl,B,alpha=.30,color=COLORS['vpl'],zorder=6)
    silhouette(ax,vpl,B,color='#956813',lw=.75,alpha=.95,zorder=7)
    target=TARGET@B;center=atlas.vpl_centroid@B
    ax.scatter(target[0],target[1],s=25,marker='x',c='white',linewidth=2.3,zorder=9)
    ax.scatter(target[0],target[1],s=25,marker='x',c='#243b53',linewidth=.85,zorder=10)
    if labels:
        ax.annotate('Left VPL',xy=center[:2],xytext=center[:2]+np.array([-2.8,-2.2]),
                    fontsize=8,arrowprops=dict(arrowstyle='-',lw=.7,color='#735623'),zorder=20)
        ax.annotate('Requested target',xy=target[:2],xytext=target[:2]+np.array([.8,2.0]),
                    fontsize=7.5,arrowprops=dict(arrowstyle='-',lw=.7,color=COLORS['ink']),zorder=20)
    if montage:
        pairs=crossed_montage();electrodes=np.array([p.as_array() for pair in pairs.values() for p in pair])
        plane=np.array([[-4.8,-4.4,.5],[1.,-4.4,.5],[1.,.9,.5],[-4.8,.9,.5]])@B
        ax.add_patch(plt.Polygon(plane[:,:2],facecolor='#eaf1f5',alpha=.32,edgecolor='#9daab4',lw=.7,zorder=9))
        for i,p in enumerate(electrodes):
            c=COLORS['ch1'] if i<2 else COLORS['ch2']
            pos=mark(ax,p,B,c,size=48,zorder=12)
            ax.text(pos[0]+.20,pos[1]+.20,f'E{i+1}',fontsize=7.5,color=c,weight='bold',zorder=13)
        for pair,c in [(electrodes[:2],COLORS['ch1']),(electrodes[2:],COLORS['ch2'])]:
            line=pair@B;ax.plot(line[:,0],line[:,1],color=c,lw=.8,ls='--',zorder=11)
    pts=brain['v']@B
    xmin,ymin=pts[:,:2].min(0);xmax,ymax=pts[:,:2].max(0)
    if montage:
        q=electrodes@B;ymax=max(ymax,q[:,1].max()+.55)
    if zoom:
        mid=(target[:2]+center[:2])/2
        xmin,ymin=mid-np.array([1.6,1.6]);xmax,ymax=mid+np.array([1.6,1.6])
    ax.set_xlim(xmin-.3,xmax+.3);ax.set_ylim(ymin-.3,ymax+.3)
    ax.set_aspect('equal');ax.set_anchor('N');ax.axis('off')
    if not zoom:
        x0=xmin+.3;y0=ymin+.4
        ax.plot([x0,x0+2],[y0,y0],color=COLORS['ink'],lw=1.3,zorder=30)
        ax.text(x0+1,y0-.28,'2 mm',ha='center',va='top',fontsize=6.5,zorder=30)
    return B


def plane_plot(ax, atlas, kind, coordinate=None, params=DEFAULT, quantity='mTI', norm=NORM, target=False):
    d=atlas.plane(kind,coordinate,params=params)
    z=np.ma.masked_where(~d['brain'],d[quantity]/atlas.reference)
    filled=ax.contourf(d['x'],d['y'],z,levels=np.geomspace(norm.vmin,norm.vmax,42),
                cmap=CMAP,norm=norm,extend='both',antialiased=False)
    filled.set_edgecolor('face');filled.set_linewidth(.12)
    ax.contour(d['x'],d['y'],d['brain'].astype(float),levels=[.5],colors='#4f5c66',linewidths=.55)
    if d['vpl'].any():
        ax.contour(d['x'],d['y'],d['vpl'].astype(float),levels=[.5],colors='white',linewidths=2)
        ax.contour(d['x'],d['y'],d['vpl'].astype(float),levels=[.5],colors='#332916',linewidths=.8,linestyles='--')
    if target:
        xy={'sagittal':TARGET[[0,2]],'coronal':TARGET[[1,2]],'axial':TARGET[[1,0]]}[kind]
        ax.scatter(*xy,s=32,marker='x',color='white',linewidths=2.4,zorder=5)
        ax.scatter(*xy,s=32,marker='x',color=COLORS['ink'],linewidths=.9,zorder=6)
    ax.set_aspect('equal');ax.set_anchor('N');ax.set_xlabel(d['labels'][0]);ax.set_ylabel(d['labels'][1],labelpad=2)
    ax.tick_params(length=2.5,width=.6)
    if kind=='axial':
        ax.set_xlim(-5.2,5.2);ax.set_ylim(-7.5,5.)
    else:
        ax.set_ylim(-7.8,0.)
    return d


def save(fig, out, stem, dpi=600):
    for ext in ['pdf','svg','png']:
        temporary=out/f'.{stem}.{ext}.tmp'
        fig.savefig(temporary,format=ext,dpi=dpi)
        if ext=='pdf' and not temporary.read_bytes().rstrip().endswith(b'%%EOF'):
            raise RuntimeError(f'Incomplete PDF export: {temporary}')
        temporary.replace(out/f'{stem}.{ext}')
    plt.close(fig)
    print('Saved',stem,flush=True)


def figure_1(atlas,meshes,out):
    fig=plt.figure(figsize=(7.2,3.45))
    for i,(kind,label) in enumerate([('sagittal','Lateral cutaway'),('coronal','Coronal cutaway'),('axial','Dorsal cutaway')]):
        ax=fig.add_axes([.035+i*.325,.23,.29,.62])
        brain_scene(ax,atlas,meshes,kind,labels=i==0)
        panel_title(ax,string.ascii_lowercase[i],label)
    colorbar(fig,[.28,.155,.44,.025])
    fig.text(.5,.96,'Allen CCF anatomy with calculated TI modulation',ha='center',fontsize=10,weight='medium')
    foot(fig,'Gold: left VPL; cross: requested target')
    save(fig,out,'Figure_01_3D_TI_cutaways')


def figure_2(atlas,meshes,out):
    fig=plt.figure(figsize=(7.2,5.5))
    grid=fig.add_gridspec(2,3,left=.085,right=.965,top=.91,bottom=.20,wspace=.42,hspace=.47)
    for r in range(2):
        for c,kind in enumerate(['sagittal','coronal','axial']):
            k={'sagittal':1,'coronal':0,'axial':2}[kind]
            coord=TARGET[k] if r==0 else atlas.slices[k]
            ax=fig.add_subplot(grid[r,c]);d=plane_plot(ax,atlas,kind,coord,target=r==0)
            axisname=['AP','ML','DV'][k]
            panel_title(ax,string.ascii_lowercase[r*3+c],f'{kind.capitalize()} | {axisname} {coord:+.2f}')
            if not d['vpl'].any():
                ax.text(.5,.03,'VPL absent from this section',transform=ax.transAxes,
                        ha='center',fontsize=6.1,color='#4c5965',bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=1.5))
    fig.text(.04,.965,'Sections through the requested target (a–c) and atlas VPL centroid (d–f)',fontsize=9.3)
    colorbar(fig,[.30,.105,.4,.019])
    foot(fig)
    save(fig,out,'Figure_02_target_and_VPL_sections')


def figure_3(atlas,meshes,out):
    fig=plt.figure(figsize=(7.2,5.6))
    grid=fig.add_gridspec(2,2,left=.095,right=.94,top=.90,bottom=.23,wspace=.29,hspace=.38)
    ax=fig.add_subplot(grid[0,0]);brain_scene(ax,atlas,meshes,'sagittal',montage=True,field=False)
    panel_title(ax,'a','Four-electrode montage')
    norm=LogNorm(.01,3.)
    for i,(q,title) in enumerate([('E1','Channel 1 | E1–E2'),('E2','Channel 2 | E3–E4'),('mTI','Maximum-direction TI envelope')]):
        ax=fig.add_subplot(grid[(i+1)//2,(i+1)%2]);plane_plot(ax,atlas,'coronal',quantity=q,norm=norm)
        panel_title(ax,string.ascii_lowercase[i+1],title)
        if i==0:
            ax.set_xlabel('');ax.tick_params(labelbottom=False)
    colorbar(fig,[.29,.125,.42,.02],norm=norm,
             label='Field / baseline brain maximum mTI  (log scale)',ticks=[.01,.1,1,3])
    fig.text(.5,.963,'I1 = I2 = 0.50 mA  |  f1 = 2000 Hz; f2 = 2010 Hz',ha='center',fontsize=8)
    foot(fig,'Electrode positions are ideal coordinates, not skull projections')
    save(fig,out,'Figure_03_montage_and_channel_fields')


def figure_4(atlas,meshes,out,data_dir):
    fig=plt.figure(figsize=(7.2,5.7))
    for i,(i1,i2) in enumerate([(.25,.75),(.5,.5),(.75,.25)]):
        ax=fig.add_axes([.085+i*.315,.57,.255,.32])
        p=replace(DEFAULT,channel_1_current_mA=i1,channel_2_current_mA=i2)
        plane_plot(ax,atlas,'axial',params=p)
        panel_title(ax,string.ascii_lowercase[i],f'I1:I2 = {i1:.2f}:{i2:.2f} mA')
    colorbar(fig,[.29,.465,.42,.018])
    ax=fig.add_axes([.13,.115,.77,.215])
    sweep=np.genfromtxt(data_dir/'current_sweep.csv',delimiter=',',names=True)
    for key,label,color,ls in [('vpl_mean_relative','VPL mean',COLORS['vpl'],'-'),
        ('target_relative','Target',COLORS['ink'],'--'),
        ('brain_p95_relative','Brain P95',COLORS['ch1'],':')]:
        ax.plot(sweep['I1_mA'],sweep[key],label=label,color=color,ls=ls)
    ax.set_xlabel('Channel 1 current I1 (mA); I1 + I2 = 1.00 mA')
    ax.set_ylabel('Relative mTI');ax.set_ylim(bottom=0)
    ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,1.30),ncol=3,columnspacing=1.6)
    ax.text(-.105,1.06,'d',transform=ax.transAxes,fontsize=11,weight='bold')
    fig.text(.5,.965,'Current redistribution at a fixed total current',ha='center',fontsize=10)
    foot(fig,'All fields share the same baseline normalization')
    save(fig,out,'Figure_04_current_redistribution')


def figure_5(atlas,meshes,out,data_dir):
    fig=plt.figure(figsize=(7.2,5.0))
    # At the requested target the two fields are orthogonal. Their bisector
    # maximizes projected modulation for the equal-current baseline.
    f=preview_fields(TARGET[None],DEFAULT)
    direction=(f['E1'][0]/np.linalg.norm(f['E1'][0])+f['E2'][0]/np.linalg.norm(f['E2'][0]))
    direction/=np.linalg.norm(direction)
    a=float(f['E1'][0]@direction);b=float(f['E2'][0]@direction)
    peak=a+b;t=np.linspace(0,.25,20001)
    waveform_data={'time_s':t,'projection_direction_AP_ML_DV':direction}
    for i,df in enumerate([10.,40.,100.]):
        ax=fig.add_axes([.13,.72-i*.225,.82,.16])
        signal=(a*np.cos(2*np.pi*2000*t)+b*np.cos(2*np.pi*(2000+df)*t))/peak
        envelope=np.sqrt(np.maximum(a*a+b*b+2*a*b*np.cos(2*np.pi*df*t),0))/peak
        ax.plot(t,signal,color='#b4c3cc',lw=.27,rasterized=False)
        ax.plot(t,envelope,color=COLORS['ch1'],lw=1.1)
        ax.plot(t,-envelope,color=COLORS['ch1'],lw=1.1)
        ax.set_xlim(0,.25);ax.set_ylim(-1.15,1.15);ax.set_yticks([-1,0,1])
        ax.set_ylabel('Projected field\n(normalized)',fontsize=7)
        panel_title(ax,string.ascii_lowercase[i],f'Δf = {df:g} Hz   |   f2 = {2000+df:g} Hz')
        if i<2: ax.tick_params(labelbottom=False)
        else: ax.set_xlabel('Time (s)')
        waveform_data[f'df_{int(df)}_signal']=signal
        waveform_data[f'df_{int(df)}_envelope']=envelope
    fig.text(.13,.105,'The beat rate changes; the spatial mTI map is invariant for nonzero Δf in this resistive model.',fontsize=7.5)
    fig.text(.13,.070,'Equal carrier frequencies produce no difference-frequency modulation; this is not a neural response model.',fontsize=7)
    foot(fig)
    np.savez_compressed(data_dir/'temporal_waveforms.npz',**waveform_data)
    save(fig,out,'Figure_05_temporal_envelope')


def figure_6(atlas,meshes,out):
    fig=plt.figure(figsize=(7.2,3.65))
    ax=fig.add_axes([.055,.25,.39,.62])
    B=brain_scene(ax,atlas,meshes,'sagittal',zoom=True,field=False)
    p=TARGET@B;q=atlas.nearest@B
    ax.plot([p[0],q[0]],[p[1],q[1]],ls='--',color='#c15c4d',lw=1.2,zorder=20)
    ax.scatter(q[0],q[1],s=16,color='#c15c4d',zorder=22)
    ax.annotate('Requested target',xy=p[:2],xytext=p[:2]+[.1,.65],fontsize=7.5,
                arrowprops=dict(arrowstyle='-',color=COLORS['ink'],lw=.65),zorder=23)
    ax.annotate('Nearest VPL\nvoxel centre',xy=q[:2],xytext=q[:2]+[.45,-1.0],fontsize=7,
                arrowprops=dict(arrowstyle='-',color='#c15c4d',lw=.65),zorder=23)
    ax.text(.02,.05,f'3D separation = {atlas.distance:.3f} mm',transform=ax.transAxes,fontsize=8)
    panel_title(ax,'a','Target–atlas correspondence')
    ax=fig.add_axes([.59,.28,.35,.55])
    field=atlas.fields['mTI']/atlas.reference
    for values,label,color in [(field,'Whole brain',COLORS['ch1']),
        (field[atlas.vpl_selection],'Left VPL',COLORS['vpl'])]:
        ordered=np.sort(values);ix=np.unique(np.linspace(0,len(ordered)-1,min(1500,len(ordered))).astype(int))
        ax.step(ordered[ix],(ix+1)/len(ordered),where='post',color=color,label=label,lw=1.3)
    ax.set_xscale('log');ax.set_xlim(.008,1.1);ax.set_ylim(0,1.025)
    ax.set_xlabel('mTI / baseline brain maximum');ax.set_ylabel('Fraction of voxels')
    ax.legend(frameon=False,loc='lower right');panel_title(ax,'b','Anatomical field distributions')
    fig.text(.06,.15,f'Brain: {int(atlas.brain.sum()):,} voxels; VPL: {int(atlas.vpl.sum()):,} voxels (100 µm atlas).',fontsize=7)
    fig.text(.06,.102,'Voxel distributions are descriptive model outputs, not biological replicates or uncertainty intervals.',fontsize=7)
    fig.text(.06,.065,'Distance is to a voxel centre under the generic transform; it is not a measured registration error.',fontsize=7)
    foot(fig)
    save(fig,out,'Figure_06_target_QC_and_distributions')


def make_contact_sheet(figures_dir):
    from PIL import Image,ImageDraw
    paths=sorted(figures_dir.glob('Figure_*.png'))
    canvas=Image.new('RGB',(1800,1780),'#edf1f4')
    for i,path in enumerate(paths):
        im=Image.open(path).convert('RGB');im.thumbnail((860,550))
        x=25+(i%2)*895;y=25+(i//2)*590
        canvas.paste(im,(x+(860-im.width)//2,y+(550-im.height)//2))
        ImageDraw.Draw(canvas).text((x+10,y+556),path.stem,fill='#243b53')
    canvas.save(figures_dir/'Contact_sheet.png',dpi=(150,150))


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'paper_figures')
    p.add_argument('--only',type=int,default=0);args=p.parse_args()
    out=args.output/'figures';out.mkdir(parents=True,exist_ok=True)
    data=args.output/'data';atlas=Atlas()
    meshes=[mesh(atlas.payload,k) for k in ['brain_mesh','left_vpl_mesh']]
    if not args.only or not (data/'current_sweep.csv').exists():
        summary=generate_data(data,atlas);print(json.dumps(summary,indent=2),flush=True)
    funcs=[figure_1,figure_2,figure_3,figure_4,figure_5,figure_6]
    for i,func in enumerate(funcs,1):
        if args.only and args.only!=i:continue
        if i in [4,5]:func(atlas,meshes,out,data)
        else:func(atlas,meshes,out)
    if not args.only:
        artwork=args.output/'artwork';artwork.mkdir(exist_ok=True)
        for kind in ['sagittal','coronal','axial']:
            fig=plt.figure(figsize=(5.3,4.1));ax=fig.add_axes([0,0,1,1])
            brain_scene(ax,atlas,meshes,kind)
            fig.savefig(artwork/f'{kind}_TI_artwork.png',dpi=600,transparent=True)
            plt.close(fig)
        make_contact_sheet(out)


if __name__=='__main__':main()
