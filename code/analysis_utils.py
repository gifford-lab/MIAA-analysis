import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import subprocess
from scipy.stats import ttest_rel,wilcoxon,ranksums,ttest_ind
import scipy
import subprocess

def get_kmer(d):
    if len(d.split('/')) > 2:
        if len(d.split('/')[2].split(',')) > 2:
            return d.split('/')[2].split(',')[2]
    else:
        return 'None'
    
def get_type(d):
    if 'S' in d:
        return 'scrambled'
    else:
        return 'motif'
    
def get_bg(d):
    if len(d.split('/')) > 2:
        if len(d.split('/')[2].split(',')) > 2:
            return '/'.join(d.split('/')[2].split(',')[:2])
    else:
        return 'None'

cat2desc = {
    "CAT1": "universally opening phrases with GC content between 60 - 70 %",
"CAT2": "universally opening phrases with GC content between 30 - 50 %",
"CAT3": "universally closing phrases with GC content between 60 - 70 %",
"CAT4": "universally closing phrases with GC content between 30 - 50 %",
"CAT5": "opening ES cells using one or more occurrences of one k-mer per phrase",
"CAT6": "opening ES cells using combinations of k-mers per phrase",
"CAT7": "phrases opening ED cells using one or more occurrences of one k-mer per phrase",
"CAT8": "phrases opening ED cells using combinations of k-mers per phrase",
"CAT9": "ES-Salient-TF",
"CAT10": "ES-Salient-Top",
"CAT11": "ES-Native",
"CAT12": "ED-Salient-TF",
"CAT13": "ED-Salient-Top",
"CAT14": "ED-Native",
"CAT15": "SLOT-CNN",
"CAT16": "background"
}

def ptostar(pval):
    p1 = 'n.s.'
    if pval < 0.05:
        p1 = '*'
    if pval < 0.01:
        p1 = '**'
    if pval < 0.001:
        p1 = '***'
    return p1

def plot_scatter(x,y,xlabel="",ylabel="",kind='dot',
                 title="",transform=None,
                 correlation=pearsonr):
    if correlation != None:
        R,p = correlation(x,y)
        title += ' R='+str(round(R,4))+' p='+str.format('{0:.3g}', p)
        
    if transform != None:
        x=transform(x)
        y=transform(y)
        
    if kind == 'hex':
        h = sns.jointplot(x,y,kind='hex')
        plt.suptitle(title)
        h.ax_joint.set_xlabel(xlabel)
        h.ax_joint.set_ylabel(ylabel)

    else:
        plt.scatter(x,y)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()

def plot_dotplot(data,id_vars,data_pivot,title="",plot_boxplot=False,plot_lines=True,plot_significance=True):

    def rename(d):
        nn = ""
        cell,control = d.split('/')
        if 'ED' in cell:
            nn+='DE'
        else:
            nn += 'ESC'
        if control == 'scrambled':
            nn += ' - CTRL'
        return nn
    
    data_melt = pd.melt(data,id_vars=id_vars,value_vars=['ES Dpn ratio','ED Dpn ratio'])
    data_melt['id'] = data_melt['variable'] + '/' + data_melt['control']
   
    data_melt['id_clean'] = [rename(data_id) for data_id in data_melt['id']]
    
    cpal = {'ESC - CTRL':'lightblue','ESC':'blue','DE':'red','DE - CTRL':'salmon'}
    
    if plot_boxplot:
        sns.boxplot(x='id_clean',y='value',data=data_melt,color='lightgrey')
        
    ax = sns.swarmplot(x='id_clean',y='value',data=data_melt,palette=cpal,order=['ESC - CTRL',
                                                                            'ESC',
                                                                            'DE',
                                                                            'DE - CTRL'])
    if plot_significance:
        dist1=0.05
        dist2=0.15
        es_scram = ttest_rel(data_pivot['ES Dpn ratio'].motif,
                             data_pivot['ES Dpn ratio'].scrambled)
        ed_scram = ttest_rel(data_pivot['ED Dpn ratio'].motif,
                             data_pivot['ED Dpn ratio'].scrambled)
        es_ed = ttest_rel(data_pivot['ES Dpn ratio'].motif,
                          data_pivot['ED Dpn ratio'].motif)
        p1 = ptostar(es_ed[1])
        p2 = ptostar(es_scram[1])
        p3 = ptostar(ed_scram[1])
        x1,x2 = 2,3
        y, h, col = 1.0 + dist1, dist1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col,clip_on=False)
        plt.text((x1+x2)*0.5, y+h, p3, ha='center', va='bottom', color=col)
        x1,x2 = 0,1
        y, h, col = 1.0 + dist1, dist1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col,clip_on=False)
        plt.text((x1+x2)*0.5, y+h, p2, ha='center', va='bottom', color=col)
        x1,x2 = 1,2
        y, h, col = 1.0 + dist2, dist1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col,clip_on=False)
        plt.text((x1+x2)*0.5, y+h, p1, ha='center', va='bottom', color=col)
        plt.xlabel('')
    if plot_lines:
        for i in range(len(data_pivot)):
            plt.plot([0,1],[data_pivot['ES Dpn ratio'].scrambled.values[i],
                            data_pivot['ES Dpn ratio'].motif.values[i]],
                     color='grey')
            plt.plot([2,3],[data_pivot['ED Dpn ratio'].motif.values[i],
                             data_pivot['ED Dpn ratio'].scrambled.values[i]],
                    color='grey')
            plt.plot([1,2],[data_pivot['ES Dpn ratio'].motif.values[i],
                            data_pivot['ED Dpn ratio'].motif.values[i]],
                     color='grey')
    plt.xticks(rotation=45,ha='right')
    plt.xlabel('')
    plt.ylabel('Openness - Dpn Ratio')
    plt.text((3)*0.5,y+3*h,title,ha='center', va='bottom')
    plt.axis([-0.5,3.5,0,1.0])
    plt.tight_layout()
    plt.subplots_adjust(top=1.10)


def mean_confidence_interval_bar(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a,axis=1), scipy.stats.sem(a,axis=1)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, h

def mean_confidence_interval_dnase(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a,axis=0), scipy.stats.sem(a,axis=0)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, h

def plot_bar_scatter(data,title=""):
    data_motif = data[data['control'] == 'motif']
    data_ctrl = data[data['control'] == 'scrambled']

    #hacky use of poorly named columns
    _,x_conf = mean_confidence_interval_bar(data[[c for c in data.columns if 'Dpn Ratio' in c]])
    _,y_conf = mean_confidence_interval_bar(data[[c for c in data.columns if 'Dpn Ratio' in c]])

    x=data['ES Dpn ratio'].values
    y=data['ED Dpn ratio'].values
    
    
    for i in range(len(data)):
        plt.plot([x[i]-x_conf[i],x[i]+x_conf[i]],[y[i],y[i]],color='grey')
        plt.plot([x[i],x[i]],[y[i]-y_conf[i],y[i]+y_conf[i]],color='grey')
        
    plt.scatter(data_motif['ES Dpn ratio'],data_motif['ED Dpn ratio'],color='cyan',label='motif')
    plt.scatter(data_ctrl['ES Dpn ratio'],data_ctrl['ED Dpn ratio'],color='black',label='random control')
    #plt.legend()
    plt.plot([0,1],[0,1],color='grey',linestyle='--')
    
    plt.title(title)
    plt.xlabel('ESC Openness - Dpn Ratio')
    plt.ylabel('DE Openness - Dpn Ratio')
    plt.tight_layout()

def plot_dnase(dnase_values,tf_order,title=""):
    dnase_order = dnase_values[dnase_values['ordering'] == tf_order]
    ES_values = dnase_order[range(0,100)].values
    ED_values = dnase_order[range(100,200)].values
    plt.plot(range(-500,500,10),np.mean(ES_values,axis=0),label='ESC',color='blue')
    mean,ci =  mean_confidence_interval_dnase(ES_values, confidence=0.95)
    ci_lower = mean-ci
    ci_upper = mean+ci
    plt.fill_between(range(-500,500,10), 
                    ci_lower, 
                    ci_upper,alpha=0.2,color='blue')
    
    plt.plot(range(-500,500,10),np.mean(ED_values,axis=0),label='DE',color='red')
    mean,ci =  mean_confidence_interval_dnase(ED_values, confidence=0.95)
    ci_lower = mean-ci
    ci_upper = mean+ci
    plt.fill_between(range(-500,500,10), 
                    ci_lower, 
                    ci_upper,alpha=0.2,color='red')
   
    plt.xlabel('distance from center')
    plt.ylabel('average DNAse signal')
    plt.title( title+'N='+str(ES_values.shape[0]))
    plt.axis([-502,502,0,100])
    plt.tight_layout()
    
def average_ac(file,tf_names):
    ES_bigwig = '/cluster/krismer/projects/scm/data/esc-diff/D0/ESC-D0/signal/macs2/pooled_rep/D0_50-100_130801.bwa.mm10_pooled.pf.pval.signal.bigwig'
    ED_bigwig = '/cluster/krismer/projects/scm/data/esc-diff/D6/ESC-D6/signal/macs2/pooled_rep/D5_endo_50-100_121229.bwa.mm10_pooled.pf.pval.signal.bigwig'
    posid_tfs = {}
    header=True
    tf_order_counts = {}
    for i,line in enumerate(open(file)):
        if header:
            header=False
            continue
        data = line.strip().split()
        for tf_name in tf_names:
            if tf_name in data[3]:
                try:
                    curr = posid_tfs[data[0]]
                except KeyError:
                    posid_tfs[data[0]] = []
                posid_tfs[data[0]].append((data[3],data[1]))
    tfs_ordered=[]
    with open('tmp.txt','w') as f:
        for pos,vlist in [(k,posid_tfs[k]) for k in posid_tfs.keys() if len(posid_tfs[k]) == len(tf_names)]:
            tfs = [vp[0] for vp in vlist]
            contains_tf = [False for _ in range(len(tf_names))]
            for tf in tfs:
                for i,tf_name in enumerate(tf_names):
                    if tf_name in tf:
                        contains_tf[i] = True
            if np.all(np.array(contains_tf)):
                tfs_order = [t[0].split('(')[0] for t in sorted(vlist,key=lambda x:int(x[1]))]
                tfs_ordered.append(','.join(tfs_order))
                f.write(pos.split(':')[0] + '\t' + pos.split(':')[1].split('-')[0]  +'\t'+ pos.split(':')[1].split('-')[1].strip() + '\n')
    subprocess.call(['/data/cgs/jhammelm/env/miniconda3/bin/computeMatrix reference-point -S '+ES_bigwig + ' '+ED_bigwig+' -o '+' results.npz --outFileNameMatrix tmp.counts -R tmp.txt -a 500 -b 500 --binSize 10']
                ,shell=True)
    table = pd.read_table('tmp.counts',skiprows=3,header=None)
    table['ordering'] = tfs_ordered
    return table
