#!/users/mscherer/software/anaconda3/envs/selection/bin/python3

import argparse as ap
import os.path
import sys
import pandas as pd
import numpy as np
import yaml
import re
import subprocess
import io
from datetime import datetime, date

pars = ap.ArgumentParser()
pars.add_argument('-s','--sample_annotation',type=str,help='Path to the sample annotation sheet',
    default=None)
pars.add_argument('-p','--path',type=str,help="Path to the data directory",
    default=None)
pars.add_argument('-r','--rnbSet',type=str,help="Path to an existing rnbSet",
    default=None)
pars.add_argument('-d','--differential',type=str,help="Path to an existing differential methylation RnBeads run (folder differential_methylation_data)",
    default=None)
pars.add_argument('-c','--config',type=str,help="Path to the configuration YAML file")
pars.add_argument('-o','--output',type=str,help="Path to the output directory")

args = pars.parse_args()

if not os.path.isfile(args.config):
    sys.exit('Config file at ' + args.config + ' does not exist')        

with open(args.config) as conf: 
    config = yaml.safe_load(conf)

if args.sample_annotation is not None:
    if not os.path.isfile(args.sample_annotation):
        sys.exit('Sample annotation sheet at ' + args.sample_annotation + ' does not exist')

    if not {'sampleID','bw_file','cov_file','cell_type'}.issubset(config['sample_annotation'].keys()):
        sys.exit('Neccessary key ' + str({'sampleID','bw_file','cov_file','cell_type'}.difference(config.keys())) + ' not present in configuration file')

    s_anno = pd.read_csv(args.sample_annotation)     

if args.path is not None:    
    if not os.path.isdir(args.path):
        sys.exit('Data directory at ' + args.path + ' does not exist')        

# convert the bigWig files to bedGraph using bigWigToBedGraph
def convert_bigWigs(s_anno, config, path, out):
    bw_files = [str(n) for n in s_anno[config['sample_annotation']['bw_file']].array]
    cov_files = [str(n) for n in s_anno[config['sample_annotation']['cov_file']].array]
    for f in bw_files + cov_files:
        print(str(datetime.now()) + ': Start converting bigWigs for ' + f)
        if not os.path.isfile(path + f):
            sys.exit('Data file ' + f + ' does not exist')

        bedGraph = re.sub('.bw','.bedGraph',f)
        cmd = 'bigWigToBedGraph ' + path + f + ' ' + out + bedGraph
        proc = subprocess.run(cmd,shell=True)

    all_beds = []
    for i in range(len(bw_files)):
        print(str(datetime.now()) + ': Start generating BEDs for ' + bw_files[1])
        bw_bed = re.sub('.bw','.bedGraph',bw_files[i])
        cov_bed = re.sub('.bw','.bedGraph',cov_files[i])
        out_bed = out + 'combined_bed_' + s_anno[config['sample_annotation']['sampleID']][i] + '.bed'
        cmd = 'bedtools intersect -a ' + out + bw_bed + ' -b ' + out + cov_bed + ' -wa -wb -f 1.0 > ' + out_bed
        proc = subprocess.run(cmd,shell=True)
        out_fr = pd.read_table(out_bed,header=None)
        out_fr[8] = [int(round(n)) for n in out_fr[3].array*out_fr[7].array]
        out_fr[9] = out_fr[7].array - out_fr[8].array
        out_fr[3] = [int(round(n*100)) for n in out_fr[3].array]
        out_fr[2] = out_fr[2] - 1
        out_fr[1] = out_fr[1] + 1
        out_fr[[0,1,2,3,8,9]].to_csv(out_bed,sep='\t',header=False,index=False)
        os.remove(out + bw_bed)
        os.remove(out + cov_bed)
        all_beds.append('combined_bed_' + s_anno[config['sample_annotation']['sampleID']][i] + '.bed')

    s_anno['bismarkCov'] = all_beds
    s_anno.to_csv(out + 'sample_annotation.csv',index=False)

def run_RnBeads(config, out, path=None):
    print(str(datetime.now()) + ': Start generating RnBeads command')
    if path is None:
        path = out

    s_anno = path + 'sample_annotation.csv'
    out_name = out + "report_" + str(date.today())
    script = ["library(RnBeads)",
        "s.anno <- '" + s_anno + "'", 
        "dir.report <- '" + out_name + "'",
        "dat.dir <- '" + path + "'",
        "rnb.options(identifiers.column='" + config['sample_annotation']['sampleID'] + "',",
        "import.bed.style='bismarkCov',",
        "min.group.size=1,",
        "max.group.count=99,",
        "exploratory.region.profiles='',",
        "assembly='" + config["rnbeads"]["assembly"] + "',",
        "filtering.coverage.threshold=" + str(config["rnbeads"]["coverage"]) + ",",
        "filtering.low.coverage.masking=TRUE,",
        "filtering.missing.value.quantile=0,",
        "filtering.high.coverage.outliers=TRUE,",
        "filtering.sex.chromosomes.removal=TRUE,",
        "qc=FALSE,",
        "export.to.bed=FALSE,",
        "export.to.trackhub=NULL,",
        "exploratory=FALSE,",
        "differential.comparison.columns='" + config['sample_annotation']['cell_type'] + "')",
        "rnb.run.analysis(dir.reports=dir.report,data.source=c(dat.dir,s.anno,'bismarkCov'),data.type='bs.bed.dir')" 
    ]
    script = "\n".join(script)
    script_file = open(out + "rnbeads_script.R", "w")
    script_file.write(script)
    script_file.close()
    cmd = "Rscript " + out + "rnbeads_script.R"
    print(str(datetime.now()) + ': Start running RnBeads')
    proc = subprocess.run(cmd,shell=True)
    print(str(datetime.now()) + ': Completed running RnBeads')
    return out_name

def run_RnBeads_rnbSet(rnbSet, config, out):
    print(str(datetime.now()) + ': Start generating RnBeads command')
    tmp_dir = out + 'tmp'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    out_name = out + "report_" + str(date.today())
    script = ["library(RnBeads)",
        "options(fftempdir='" + tmp_dir + "')",
        "rnb.set <- load.rnb.set('" + rnbSet + "')", 
        "rnb.load.annotation.from.db(setdiff(summarized.regions(rnb.set),rnb.region.types(assembly='" + config['rnbeads']['assembly'] + "')),assembly='" + config['rnbeads']['assembly'] + "')",
        "dir.report <- '" + out_name + "'",
        "rnb.options(identifiers.column='" + config['rnbeads']['sampleID'] + "',",
        "min.group.size=1,",
        "max.group.count=99,",
        "exploratory.region.profiles='',",
        "filtering.coverage.threshold=" + str(config["rnbeads"]["coverage"]) + ",",
        "filtering.low.coverage.masking=TRUE,",
        "filtering.missing.value.quantile=0,",
        "filtering.high.coverage.outliers=TRUE,",
        "filtering.sex.chromosomes.removal=TRUE,",
        "qc=FALSE,",
        "export.to.bed=FALSE,",
        "export.to.trackhub=NULL,",
        "exploratory=FALSE,",
        "differential.comparison.columns='" + config['rnbeads']['cell_type'] + "')",
        "rnb.run.analysis(dir.reports=dir.report,data.source=rnb.set,data.type='rnb.set')" 
    ]
    script = "\n".join(script)
    script_file = open(out + "rnbeads_script.R", "w")
    script_file.write(script)
    script_file.close()
    cmd = "Rscript " + out + "rnbeads_script.R"
    print(str(datetime.now()) + ': Start running RnBeads')
    proc = subprocess.run(cmd,shell=True)
    print(str(datetime.now()) + ': Completed running RnBeads')
    return out_name

def return_reliable_DMCs(rnbOut, n, existing_dmcs, out, config):
    print(str(datetime.now()) + ': Start selecting reliable DMCs')
    diff_table = pd.read_csv(rnbOut)
    row_names = np.array(diff_table['Chromosome'] + ':' + [str(n) for n in diff_table['Start']])
    diff_table = diff_table.loc[[n not in existing_dmcs for n in row_names],]
    diff_table = diff_table.sort_values('combinedRank')
    group_names = np.array(diff_table.columns[['sd.' in n for n in diff_table.columns]].array)
    group_names = [re.sub('sd.','',n) for n in group_names]
    print('Processing diffTable for ' + group_names[0] + ' and ' + group_names[1])
    sds_first = diff_table['sd.'+group_names[0]]
    sds_second = diff_table['sd.'+group_names[1]]
    se_sd_first = np.std(sds_first)/np.sqrt(len(sds_first))
    se_sd_second = np.std(sds_second)/np.sqrt(len(sds_second))
    low_sds = np.logical_and(sds_first<(np.mean(sds_first)+2*se_sd_first),sds_second<(np.mean(sds_second)+2*se_sd_second))
    diff_table = diff_table.loc[low_sds,]
    diff_positive = diff_table.loc[diff_table['mean.diff']>0,]
    diff_positive['type'] = 'high_' + group_names[0]
    diff_positive = check_cut_sites(diff_positive, out, config, n//2)
    diff_negative = diff_table.loc[diff_table['mean.diff']<0,]
    diff_negative['type'] = 'low_' + group_names[0]
    diff_negative = check_cut_sites(diff_negative, out, config, n//2)
    print(str(datetime.now()) + ': Completed selecting reliable DMCs')
    return(pd.concat([diff_positive,diff_negative]))

def check_cut_sites(frame, out, config, n):
    print(str(datetime.now()) + ': Checking for enzyme cutsite')
    f = out + '/temp.csv'
    frame.to_csv(f)
    file_path = os.path.realpath(__file__)
    file_path = re.sub('_main_.py','check_cut_site.R',file_path)
    cmd = "Rscript " + file_path + " -i " + f + " -c " + config + " -n " + str(n)
    proc = subprocess.run(cmd,shell=True)
    res = pd.read_csv(f)
    return res

def return_constant_methylated(report, n, blacklist_dmcs, out, config):
    print(str(datetime.now()) + ': Start selecting always methylated and unmethylated CpGs')
    with open(config) as conf: 
        config_pd = yaml.safe_load(conf)

    diff_report = report + '/differential_methylation_data/'
    diff_results = np.array(os.listdir(diff_report))
    diff_results = np.array(diff_results[['diffMethTable_site' in n for n in diff_results]])
    group_names_all = []
    diff_all = pd.DataFrame()
    for result in diff_results:
        diff_table = pd.read_csv(diff_report + result)
        row_names = np.array(diff_table['Chromosome'] + ':' + [str(n) for n in diff_table['Start']])
        diff_table = diff_table.loc[[n not in blacklist_dmcs for n in row_names],]
        group_names = np.array(diff_table.columns[['sd.' in n for n in diff_table.columns]].array)
        group_names = [re.sub('sd.','',n) for n in group_names]
        group_names_all.extend(group_names)
        diff_table = diff_table.loc[diff_table['diffmeth.p.val']>0.05,]
        diff_all = pd.concat([diff_all,diff_table],axis=1)

    fully_methylated = diff_all[["mean." + n for n in group_names_all]] > config_pd['dmcs']['methylated']
    fully_methylated = fully_methylated.apply(all,axis=1)
    fully_methylated = diff_all.loc[fully_methylated]
    fully_methylated['type'] = 'methylated'
    mean_coverage_meth = fully_methylated[["mean.covg." + n for n in group_names_all]].apply(np.mean,axis=1)
    fully_methylated['combined_coverage'] = mean_coverage_meth
    fully_unmethylated = diff_all[["mean." + n for n in group_names_all]] < config_pd['dmcs']['unmethylated']
    fully_unmethylated = fully_unmethylated.apply(all,axis=1)
    fully_unmethylated = diff_all.loc[fully_unmethylated]
    fully_unmethylated['type'] = 'unmethylated'
    mean_coverage_unmeth = fully_unmethylated[["mean.covg." + n for n in group_names_all]].apply(np.mean,axis=1)
    fully_unmethylated['combined_coverage'] = mean_coverage_unmeth
    #fully_methylated = fully_methylated.sort_values('combined_coverage',ascending=False)
    #fully_unmethylated = fully_unmethylated.sort_values('combined_coverage',ascending=False)
    fully_methylated = fully_methylated.loc[fully_methylated['combined_coverage']>config_pd['dmcs']['coverage_threshold']]
    fully_methylated = fully_methylated.sample(fully_methylated.shape[0], axis=0)     
    fully_methylated = check_cut_sites(fully_methylated, out, config, n)
    fully_unmethylated = fully_unmethylated.loc[fully_unmethylated['combined_coverage']>config_pd['dmcs']['coverage_threshold']]    
    fully_unmethylated = fully_unmethylated.sample(fully_unmethylated.shape[0], axis=0)     
    fully_unmethylated = check_cut_sites(fully_unmethylated, out, config, n)
    return(pd.concat([fully_methylated,fully_unmethylated]))

def write_non_cut(out,config):
    print(str(datetime.now()) + ': Start constructing amplicons that are not cut')
    out = out + "/non_cut_amplicons.csv"
    file_path = os.path.realpath(__file__)
    file_path = re.sub('_main_.py','construct_non_cut.R',file_path)
    cmd = "Rscript " + file_path + " -c " + config + " -o " + out
    proc = subprocess.run(cmd,shell=True)
    res = pd.read_csv(out)
    return res

if args.differential is not None:
    report = args.differential
    re.sub("/differential_methylation_data/", "", report)

elif args.rnbSet is None:
    if args.path is None:
        sys.exit('Either --path or --rnbSet have to be specified') 

    if not {'type'}.issubset(config['input'].keys()):
        sys.exit('Neccessary key type in input category not present in configuration file')

    if config['input']['type'] == 'bigWig':
        convert_bigWigs(s_anno, config, args.path,args.output)
        report = run_RnBeads(config, args.output)

    elif config['input']['type'] == 'bed':
        report = run_RnBeads(config, args.output, args.path)

    else:
        sys.exit('Currently only bigWig and BED files possible inputs')
else:
    report = run_RnBeads_rnbSet(args.rnbSet, config, args.output)

diff_report = report + '/differential_methylation_data/'
diff_results = np.array(os.listdir(diff_report))
diff_results = np.array(diff_results[['diffMethTable_site' in n for n in diff_results]])
num_dmcs = config['dmcs']['number']//len(diff_results)
res_all = pd.DataFrame()
if 'amplicon_blacklist' in config['dmcs'].keys():
    blacklist = pd.read_csv(config['dmcs']['amplicon_blacklist'], sep='\t')
    ids = blacklist['Chromosome'] + ':' + [str(n) for n in blacklist['Start']]

else:
    ids = np.array("")

for diff_res in diff_results:
    if ids!=np.array(""):
        np.append(ids,res_all['Chromosome'] + ':' + [str(n) for n in res_all['Start']])

    res = return_reliable_DMCs(diff_report + diff_res, num_dmcs, np.array(ids), args.output, args.config)
    res_all = pd.concat([res_all,res])

res_all.to_csv(args.output + 'reliable_DMRs.csv')
np.append(ids,res_all['Chromosome'] + ':' + [str(n) for n in res_all['Start']])
meth_unmeth = return_constant_methylated(report, config['dmcs']['number_meth_unmeth'], ids, args.output, args.config)
meth_unmeth.to_csv(args.output + 'meth_unmeth.csv')
res_all = pd.read_csv(args.output + 'reliable_DMRs.csv')
meth_unmeth = pd.read_csv(args.output + 'meth_unmeth.csv')
non_cut = write_non_cut(args.output,args.config)
non_cut.columns = ['ID','Location','Chromosome','regionStart','regionEnd','CpGCount','AciSite','Type','GCContent']
non_cut['type'] = 'non-cut'
non_cut['cutsiteInRegion'] = 'none'
non_cut['Start'] = non_cut['regionStart']
final_frame = pd.concat([res_all[['Chromosome','Start','type','regionStart','regionEnd','cutsiteInRegion']],
    meth_unmeth[['Chromosome','Start','type','regionStart','regionEnd','cutsiteInRegion']],
    non_cut[['Chromosome','Start','type','regionStart','regionEnd','cutsiteInRegion']]])
final_frame.to_csv(args.output + 'final_output.csv')

