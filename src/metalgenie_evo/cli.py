"""
MetalGenie-Evo  —  HMM-based annotation of iron cycling and metal resistance genes
===================================================================================
Built on FeGenie (Garber et al. 2020). See README for full documentation.

New in this version (metagenome support):
  7.  Integrated Prodigal  --fna_dir triggers internal ORF prediction
  8.  Contig length filter  --min_contig_len
  9.  Relaxed operon thresholds  --relaxed_operons
  10. TPM coverage normalisation  --norm_coverage
  11. Contig column in all outputs
  12. Long-format tidy output  MetalGenie-Evo-results-long.tsv
"""

import argparse, csv, json, os, re, sys, fnmatch, shutil, subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ── Default operon rules ──────────────────────────────────────────────────────
_DEFAULT_OPERON_RULES = [
    {"name":"FLEET","categories":["iron_oxidation"],
     "genes":["EetA","EetB","Ndh2","FmnB","FmnA","DmkA","DmkB","PplA"],
     "rule":"require_n_of","min_genes":5,"on_fail":"passthrough_non_members"},
    {"name":"MAM","categories":["magnetosome_formation"],
     "genes":["MamA","MamB","MamE","MamK","MamP","MamM","MamQ","MamI","MamL","MamO"],
     "rule":"require_n_of","min_genes":5,"on_fail":"passthrough_non_members"},
    {"name":"FOXABC","categories":["iron_oxidation"],
     "genes":["FoxA","FoxB","FoxC"],
     "rule":"require_n_of","min_genes":2,"on_fail":"passthrough_non_members"},
    {"name":"FOXEYZ","categories":["iron_oxidation"],
     "genes":["FoxE","FoxY","FoxZ"],
     "rule":"require_anchor","anchor":"FoxE","on_fail":"passthrough_non_members"},
    {"name":"DFE1","categories":["iron_reduction","probable_iron_reduction"],
     "genes":["DFE_0448","DFE_0449","DFE_0450","DFE_0451"],
     "rule":"require_n_of","min_genes":3,"on_fail":"passthrough_non_members"},
    {"name":"DFE2","categories":["iron_reduction","probable_iron_reduction"],
     "genes":["DFE_0461","DFE_0462","DFE_0463","DFE_0464","DFE_0465"],
     "rule":"require_n_of","min_genes":3,"on_fail":"passthrough_non_members"},
    {"name":"MtrMto","categories":["iron_oxidation","iron_reduction",
     "possible_iron_oxidation_and_possible_iron_reduction"],
     "genes":["MtrA","MtrB_TIGR03509","MtrC_TIGR03507","MtoA","CymA"],
     "rule":"mtr_disambiguation","on_fail":"keep_all"},
    {"name":"SIDERO_TRANSPORT","categories":["iron_aquisition-siderophore_transport_potential",
     "iron_aquisition-heme_transport","iron_aquisition-siderophore_transport"],
     "genes":[],"rule":"require_n_cat_or_lone_trusted","min_genes":2,
     "trusted_lone":["FutA1-iron_ABC_transporter_iron-binding-rep",
                     "FutA2-iron_ABC_transporter_iron-binding-rep",
                     "FutC-iron_ABC_transporter_ATPase-rep",
                     "LbtU-LvtA-PiuA-PirA-RhtA","LbtU-LbtB-legiobactin_receptor",
                     "LbtU_LbtB-legiobactin_receptor_2","IroC-salmochelin_transport-rep"],
     "on_fail":"drop"},
    {"name":"SIDERO_SYNTH","categories":["iron_aquisition-siderophore_synthesis"],
     "genes":[],"rule":"require_n_cat","min_genes":3,"on_fail":"drop"},
    {"name":"IRON_TRANSPORT","categories":["iron_aquisition-iron_transport",
     "iron_aquisition-heme_oxygenase"],
     "genes":[],"rule":"require_n_cat","min_genes":2,"on_fail":"drop"},
]
_REPORT_ALL_PATTERNS = ["metal_resistance-*","iron_storage"]

# ── IO ────────────────────────────────────────────────────────────────────────
def read_fasta(path):
    seqs,header,parts={},None,[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith(">"):
                if header is not None: seqs[header]="".join(parts)
                header,parts=line[1:].split()[0],[]
            else: parts.append(line)
    if header is not None: seqs[header]="".join(parts)
    return seqs

def read_fasta_lengths(path):
    lengths,header,length={},None,0
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith(">"):
                if header is not None: lengths[header]=length
                header=line[1:].split()[0]; length=0
            else: length+=len(line)
    if header is not None: lengths[header]=length
    return lengths

def read_cutoffs(path):
    c={}
    if not os.path.isfile(path): return c
    with open(path) as fh:
        for line in fh:
            ls=line.rstrip().split("\t")
            if len(ls)>=2:
                try: c[ls[0]]=float(ls[1])
                except ValueError: pass
    return c

def read_map(path):
    m={}
    if not os.path.isfile(path): return m
    with open(path) as fh:
        for line in fh:
            ls=line.rstrip().split("\t")
            if len(ls)>=2: m[ls[0]]=ls[1]
    return m

def load_operon_rules(hmm_dir):
    """Return (rules, report_all_pats, json_present).
    json_present=True → use JSON rule engine (model organisms).
    json_present=False → use FeGenie exact port (default).
    """
    p=Path(hmm_dir)/"operon_rules.json"
    if p.exists():
        with open(p) as fh: data=json.load(fh)
        rules=data.get("rules",[])
        pats=data.get("report_all_categories",_REPORT_ALL_PATTERNS)
        print(f"[INFO] Loaded {len(rules)} operon rules from {p}")
        print(f"       Using JSON rule engine (model-organism mode)")
        return rules,pats,True
    print("[INFO] Using FeGenie exact operon logic (no operon_rules.json found)")
    return [],_REPORT_ALL_PATTERNS,False

# ── Prodigal ──────────────────────────────────────────────────────────────────
def _prodigal_job(args_tuple):
    """Top-level function (picklable) for ProcessPoolExecutor."""
    fna_path, faa_out, gff_out, meta_mode = args_tuple
    stem = Path(fna_path).stem
    if Path(faa_out).exists() and Path(gff_out).exists():
        return stem, True, ""
    cmd = ["prodigal", "-i", str(fna_path),
           "-a", faa_out, "-f", "gff", "-o", gff_out, "-q"]
    if meta_mode:
        cmd += ["-p", "meta"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return stem, r.returncode == 0, r.stderr if r.returncode else ""


def run_prodigal(fna_files, out_dir, meta_mode=False, threads=1):
    prodigal_dir = out_dir / "_prodigal"
    faa_dir = prodigal_dir / "faa"
    gff_dir = prodigal_dir / "gff"
    faa_dir.mkdir(parents=True, exist_ok=True)
    gff_dir.mkdir(exist_ok=True)

    # Build job tuples (all picklable — no closures)
    jobs = [
        (str(fna), str(faa_dir / f"{fna.stem}.faa"),
         str(gff_dir / f"{fna.stem}.gff"), meta_mode)
        for fna in fna_files
    ]

    print(f"[INFO] Running Prodigal on {len(fna_files)} assemblies "
          f"({'meta' if meta_mode else 'single'} mode)…")

    n_workers = min(threads, len(fna_files))
    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_prodigal_job, j): j for j in jobs}
            done = errors = 0
            for fut in as_completed(futures):
                done += 1
                stem, ok, err = fut.result()
                if not ok:
                    errors += 1
                    print(f"\n  [WARN] Prodigal failed {stem}: {err[:80]}",
                          file=sys.stderr)
                sys.stdout.write(
                    f"\r  {done}/{len(fna_files)}  ({errors} errors)  ")
                sys.stdout.flush()
        print()
    else:
        for i, job in enumerate(jobs, 1):
            stem, ok, err = _prodigal_job(job)
            sys.stdout.write(f"\r  {i}/{len(fna_files)}  ")
            sys.stdout.flush()
            if not ok:
                print(f"\n  [WARN] Prodigal failed {stem}: {err[:80]}",
                      file=sys.stderr)
        print()

    return faa_dir, gff_dir

# ── Contig lengths ─────────────────────────────────────────────────────────────
def get_contig_lengths_from_gff(gff_path):
    lengths=defaultdict(int)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                m=re.search(r"##sequence-region\s+(\S+)\s+\d+\s+(\d+)",line)
                if m: lengths[m.group(1)]=int(m.group(2))
                continue
            parts=line.rstrip().split("\t")
            if len(parts)>=5 and parts[2]=="CDS":
                try:
                    end=int(parts[4])
                    if end>lengths[parts[0]]: lengths[parts[0]]=end
                except ValueError: pass
    return dict(lengths)

def build_contig_length_dict(faa_files, gff_dir=None, fna_dir=None, fna_ext="fna"):
    contig_lengths={}
    for faa in faa_files:
        stem=faa.stem; genome=faa.name; found=None
        if fna_dir:
            for ext in [fna_ext,"fna","fasta","fa"]:
                c=Path(fna_dir)/f"{stem}.{ext}"
                if c.exists(): found=("fna",str(c)); break
        if not found and gff_dir:
            for ext in [".gff",".gff3",".prodigal.gff"]:
                c=Path(gff_dir)/(stem+ext)
                if c.exists(): found=("gff",str(c)); break
        if found:
            kind,path=found
            contig_lengths[genome]=(read_fasta_lengths(path) if kind=="fna"
                                     else get_contig_lengths_from_gff(path))
    return contig_lengths

# ── GFF parsing ───────────────────────────────────────────────────────────────
def load_prodigal_gff(gff_path):
    coords={}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts=line.rstrip().split("\t")
            if len(parts)<9 or parts[2]!="CDS": continue
            m=re.search(r"ID=([^;]+)",parts[8])
            if m: coords[m.group(1).strip()]={"contig":parts[0],"start":int(parts[3]),
                                               "end":int(parts[4]),"strand":parts[6]}
    return coords

def load_gff_dir(gff_dir,faa_files):
    gff_dir=Path(gff_dir); gc={}
    for faa in faa_files:
        stem=faa.stem; found=None
        for ext in (".gff",".gff3",".prodigal.gff"):
            c=gff_dir/(stem+ext)
            if c.exists(): found=c; break
        if found: gc[faa.name]=load_prodigal_gff(str(found))
        else: print(f"  [WARN] No GFF for {faa.name}, using index clustering",file=sys.stderr)
    return gc

# ── Clustering ────────────────────────────────────────────────────────────────
def _index_from_name(orf):
    parts=orf.rsplit("_",1)
    if len(parts)==2:
        try: return parts[0],int(parts[1])
        except ValueError: pass
    return orf,0

def _orf_to_contig(orf):
    c,_=_index_from_name(orf); return c

def cluster_by_index(orf_set,max_gap=5):
    by_c=defaultdict(list)
    for orf in orf_set:
        c,i=_index_from_name(orf); by_c[c].append((i,orf))
    clusters=[]
    for c,entries in by_c.items():
        entries.sort(key=lambda x:x[0]); group=[entries[0][1]]
        for i in range(1,len(entries)):
            if entries[i][0]-entries[i-1][0]<=max_gap: group.append(entries[i][1])
            else: clusters.append(group); group=[entries[i][1]]
        clusters.append(group)
    return clusters

def cluster_by_coordinates(orf_set,orf_coords,max_bp_gap=5000,strand_aware=False):
    by_c=defaultdict(list)
    for orf in orf_set:
        c=orf_coords.get(orf)
        if c is None:
            contig,idx=_index_from_name(orf); by_c[contig].append((idx*300,idx*300,"+",orf))
        else: by_c[c["contig"]].append((c["start"],c["end"],c["strand"],orf))
    clusters=[]
    for contig,entries in by_c.items():
        entries.sort(key=lambda x:x[0])
        if strand_aware:
            sg_dict=defaultdict(list)
            for e in entries: sg_dict[e[2]].append(e)
            strand_groups=list(sg_dict.values())
        else: strand_groups=[entries]
        for sg in strand_groups:
            if not sg: continue
            group=[sg[0][3]]
            for i in range(1,len(sg)):
                if sg[i][0]-sg[i-1][1]<=max_bp_gap: group.append(sg[i][3])
                else: clusters.append(group); group=[sg[i][3]]
            clusters.append(group)
    return clusters

def build_clusters(genome,orf_hits,orf_coords,max_gap,max_bp_gap,strand_aware):
    if orf_coords:
        return cluster_by_coordinates(orf_hits.keys(),orf_coords,max_bp_gap=max_bp_gap,strand_aware=strand_aware)
    return cluster_by_index(orf_hits.keys(),max_gap=max_gap)

# ── HMMER ─────────────────────────────────────────────────────────────────────
def _hmmsearch_job(args_tuple):
    hmm_file,faa_file,tblout_path,bitscore_cutoff,threads=args_tuple
    if Path(tblout_path).exists(): return tblout_path,True,""
    cmd=["hmmsearch","--cpu",str(threads),"-T",str(max(bitscore_cutoff,0)),
         "--tblout",tblout_path,"--noali","-o","/dev/null",hmm_file,faa_file]
    r=subprocess.run(cmd,capture_output=True,text=True)
    return tblout_path,r.returncode==0,r.stderr if r.returncode else ""

def parse_tblout(path):
    hits=[]
    if not os.path.isfile(path): return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts=line.split()
            if len(parts)<6: continue
            try: e=float(parts[4]); b=float(parts[5])
            except ValueError: continue
            if e<0.1: hits.append((parts[0],e,b))
    return hits

def run_all_hmmsearches(faa_files,cat_hmms,cutoffs,out_tmp,threads_total,hmm_threads=1):
    jobs=[(str(hmm_path),str(faa),str(out_tmp/f"{faa.name}__{stem}.tblout"),
           cutoffs.get(stem,0),hmm_threads)
          for faa in faa_files for cat,hmm_list in cat_hmms.items()
          for stem,hmm_path in hmm_list]
    n_workers=max(1,threads_total//hmm_threads); total=len(jobs); done=errors=0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures={pool.submit(_hmmsearch_job,j):j for j in jobs}
        for fut in as_completed(futures):
            done+=1; _,ok,err=fut.result()
            if not ok: errors+=1; print(f"\n  [WARN] hmmsearch: {err[:80]}",file=sys.stderr)
            sys.stdout.write(f"\r[INFO] hmmsearch: {done}/{total}  ({errors} errors)  "); sys.stdout.flush()
    print()

def collect_best_hits(faa_files,cat_hmms,out_tmp):
    bh=defaultdict(dict)
    for faa in faa_files:
        genome=faa.name
        for cat,hmm_list in cat_hmms.items():
            for stem,_ in hmm_list:
                for orf,ev,bs in parse_tblout(str(out_tmp/f"{genome}__{stem}.tblout")):
                    prev=bh[genome].get(orf)
                    if prev is None or bs>prev["bitscore"]:
                        bh[genome][orf]={"hmm_stem":stem,"cat":cat,"evalue":ev,"bitscore":bs}
    return bh

# ── Operon filtering — FeGenie exact port ────────────────────────────────────
def _cm(cat,pats): return any(fnmatch.fnmatch(cat,p) for p in pats)

def _mtr(rows):
    stems={r["hmm_stem"] for r in rows}; updated=[dict(r) for r in rows]
    if "MtoA" in stems and "MtrB_TIGR03509" in stems and "MtrC_TIGR03507" not in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrB_TIGR03509","MtoA","CymA"}: r["cat"]="iron_oxidation"
    elif "MtrA" in stems and "MtrB_TIGR03509" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA","MtrB_TIGR03509"}: r["cat"]="iron_reduction"
    elif "MtrC_TIGR03507" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA","MtrB_TIGR03509"}: r["cat"]="iron_reduction"
    return updated

# FeGenie gene sets
_FLEET_GENES   = {"EetA","EetB","Ndh2","FmnB","FmnA","DmkA","DmkB","PplA"}
_MAM_GENES     = {"MamA","MamB","MamE","MamK","MamP","MamM","MamQ","MamI","MamL","MamO"}
_FOXABC_GENES  = {"FoxA","FoxB","FoxC"}
_FOXEYZ_GENES  = {"FoxE","FoxY","FoxZ"}
_DFE1_GENES    = {"DFE_0448","DFE_0449","DFE_0450","DFE_0451"}
_DFE2_GENES    = {"DFE_0461","DFE_0462","DFE_0463","DFE_0464","DFE_0465"}
_MTR_GENES     = {"MtrA","MtrB_TIGR03509","MtrC_TIGR03507","MtoA","CymA"}
_TRUSTED_LONE  = {"FutA1-iron_ABC_transporter_iron-binding-rep",
                  "FutA2-iron_ABC_transporter_iron-binding-rep",
                  "FutC-iron_ABC_transporter_ATPase-rep",
                  "LbtU-LvtA-PiuA-PirA-RhtA","LbtU-LbtB-legiobactin_receptor",
                  "LbtU_LbtB-legiobactin_receptor_2","IroC-salmochelin_transport-rep"}

_SIDERO_TRANS_CATS = {"iron_aquisition-siderophore_transport_potential",
                      "iron_aquisition-siderophore_transport",
                      "iron_aquisition-heme_transport"}
_SIDERO_SYNTH_CATS = {"iron_aquisition-siderophore_synthesis"}
_IRON_TRANS_CATS   = {"iron_aquisition-iron_transport","iron_aquisition-heme_oxygenase"}
_IRON_ACQ_ALL      = _SIDERO_TRANS_CATS | _SIDERO_SYNTH_CATS | _IRON_TRANS_CATS


def filter_cluster_fegenie(cluster_rows, report_all_pats, all_results=False):
    """
    Exact port of FeGenie's operon context filtering logic.

    FeGenie dispatches clusters to one rule handler based on which special
    genes are present.  The key behavioral difference from a simple threshold:

      - `break`-style rules  → drop the ENTIRE cluster when threshold not met
      - `pass`-style rules   → drop only the FAILING ORFs, keep the rest

    iron_transport/siderophore rules use the pass pattern, meaning that a gene
    in one of those categories is silently skipped if its cluster doesn't have
    enough co-occurring genes of the same category — but OTHER genes in the
    cluster (e.g. iron_reduction genes) are kept.
    """
    if all_results:
        return cluster_rows

    cats  = {r["cat"]      for r in cluster_rows}
    stems = {r["hmm_stem"] for r in cluster_rows}

    # report_all bypass (metal_resistance-*, iron_storage …)
    if all(_cm(c, report_all_pats) for c in cats):
        return cluster_rows

    # ── Helper counts ─────────────────────────────────────────────────────────
    def n_unique_total():
        return len(stems)

    def n_unique_in_cats(cat_set):
        return len({r["hmm_stem"] for r in cluster_rows if r["cat"] in cat_set})

    def n_in_gene_set(gene_set):
        return len({r["hmm_stem"] for r in cluster_rows if r["hmm_stem"] in gene_set})

    # ── Dispatch: FLEET ───────────────────────────────────────────────────────
    if stems & _FLEET_GENES:
        n = n_in_gene_set(_FLEET_GENES)
        non_fleet = [r for r in cluster_rows if r["hmm_stem"] not in _FLEET_GENES]
        if n >= 5:
            return cluster_rows
        return non_fleet if non_fleet else []

    # ── MAM ───────────────────────────────────────────────────────────────────
    if stems & _MAM_GENES:
        n = n_in_gene_set(_MAM_GENES)
        non_mam = [r for r in cluster_rows if r["hmm_stem"] not in _MAM_GENES]
        if n >= 5:
            return cluster_rows
        return non_mam if non_mam else []

    # ── FoxABC ────────────────────────────────────────────────────────────────
    if stems & _FOXABC_GENES:
        n = n_in_gene_set(_FOXABC_GENES)
        non_fox = [r for r in cluster_rows if r["hmm_stem"] not in _FOXABC_GENES]
        if n >= 2:
            return cluster_rows
        return non_fox if non_fox else []

    # ── FoxEYZ ────────────────────────────────────────────────────────────────
    if stems & _FOXEYZ_GENES:
        non_fox = [r for r in cluster_rows if r["hmm_stem"] not in _FOXEYZ_GENES]
        if "FoxE" in stems:
            return cluster_rows
        return non_fox if non_fox else []

    # ── DFE1 ─────────────────────────────────────────────────────────────────
    if stems & _DFE1_GENES:
        n = n_in_gene_set(_DFE1_GENES)
        non_dfe = [r for r in cluster_rows if r["hmm_stem"] not in _DFE1_GENES]
        if n >= 3:
            return cluster_rows
        return non_dfe if non_dfe else []

    # ── DFE2 ─────────────────────────────────────────────────────────────────
    if stems & _DFE2_GENES:
        n = n_in_gene_set(_DFE2_GENES)
        non_dfe = [r for r in cluster_rows if r["hmm_stem"] not in _DFE2_GENES]
        if n >= 3:
            return cluster_rows
        return non_dfe if non_dfe else []

    # ── Cyc1 ─────────────────────────────────────────────────────────────────
    # (handled in second_pass, not here)

    # ── Mtr / Mto disambiguation ──────────────────────────────────────────────
    if "CymA" in stems:
        if stems & {"MtrA","MtoA","MtrB_TIGR03509","MtrC_TIGR03507"}:
            return _mtr(cluster_rows)
        return [r for r in cluster_rows if r["hmm_stem"] != "CymA"] or []

    if stems & _MTR_GENES:
        return _mtr(cluster_rows)

    # ── Iron acquisition categories (per-orf pass/break logic) ────────────────
    if cats & _IRON_ACQ_ALL:
        n_total = n_unique_total()
        kept = []
        skip_cluster = False

        for r in cluster_rows:
            cat = r["cat"]

            if cat in _SIDERO_TRANS_CATS:
                if n_total < 2:
                    skip_cluster = True; break          # break → drop cluster
                if n_unique_in_cats(_SIDERO_TRANS_CATS) < 2:
                    pass                                 # pass  → skip this orf
                else:
                    kept.append(r)

            elif cat in _SIDERO_SYNTH_CATS:
                if n_total < 3:
                    skip_cluster = True; break
                if n_unique_in_cats(_SIDERO_SYNTH_CATS) < 3:
                    pass
                else:
                    kept.append(r)

            elif cat in _IRON_TRANS_CATS:
                if n_total < 2:
                    skip_cluster = True; break
                if n_unique_in_cats(_IRON_TRANS_CATS) < 2:
                    pass
                else:
                    kept.append(r)

            else:
                # Non-iron-aquisition gene co-occurring in the cluster — keep it
                kept.append(r)

        if skip_cluster:
            return []

        # FeGenie's "else" branch for iron_aquisition clusters:
        # if no kept ORFs from the category checks but cluster has >1 unique HMM,
        # check for trusted lone receptors
        if not kept:
            if n_total > 1 or (stems & _TRUSTED_LONE):
                return cluster_rows
            return []

        return kept

    # ── All other categories (iron_reduction, iron_oxidation, etc.) ───────────
    # FeGenie reports these as long as the cluster has ≥1 orf — no additional rule
    return cluster_rows


def filter_cluster_json(cluster_rows, operon_rules, report_all_pats,
                        all_results=False, contig_len=None, relaxed_threshold=10000):
    """
    JSON-rule engine (for operon_rules.json users — model organisms).
    Used when operon_rules.json is present in --hmm_dir.
    """
    if all_results: return cluster_rows
    cats={r["cat"] for r in cluster_rows}
    if all(_cm(c,report_all_pats) for c in cats): return cluster_rows
    relaxed=(contig_len is not None and contig_len>0 and contig_len<relaxed_threshold)

    def _ugi(rows,gs): return len({r["hmm_stem"] for r in rows if r["hmm_stem"] in gs})
    def _ric(rows,cs): return [r for r in rows if r["cat"] in cs]

    def _apply(rows,rd):
        gs=set(rd.get("genes",[])); cs=set(rd.get("categories",[]))
        on_fail=rd.get("on_fail","keep_all"); rt=rd["rule"]
        sp={r["hmm_stem"] for r in rows}; cp={r["cat"] for r in rows}
        if not ((gs and sp&gs) or (cs and cp&cs)): return rows
        raw_min=rd.get("min_genes",1)
        min_n=max(1,raw_min//2) if relaxed else raw_min
        non_mbrs=[r for r in rows if r["hmm_stem"] not in gs]
        if rt=="require_n_of":
            if _ugi(rows,gs)>=min_n: return rows
            if on_fail=="passthrough_non_members" and non_mbrs: return non_mbrs
            return rows if on_fail=="keep_all" else []
        if rt=="require_anchor":
            anchor=rd.get("anchor","")
            if anchor in {r["hmm_stem"] for r in rows if r["hmm_stem"] in gs}: return rows
            if on_fail=="passthrough_non_members" and non_mbrs: return non_mbrs
            return rows if on_fail=="keep_all" else []
        if rt=="require_n_cat":
            if len({r["hmm_stem"] for r in _ric(rows,cs)})>=min_n: return rows
            return rows if on_fail=="keep_all" else []
        if rt=="require_n_cat_or_lone_trusted":
            tl=set(rd.get("trusted_lone",[])); cat_m=_ric(rows,cs)
            uq={r["hmm_stem"] for r in cat_m}
            if len(uq)>1 or tl&sp or len(uq)>=min_n: return rows
            return rows if on_fail=="keep_all" else []
        if rt=="mtr_disambiguation": return _mtr(rows)
        return rows

    rows=cluster_rows
    for rd in operon_rules:
        rows=_apply(rows,rd)
        if not rows: return []
    return rows


FE_REDOX={"iron_reduction","iron_oxidation"}

def count_heme(seq):
    if not seq: return 0
    return sum(len(re.findall(p,seq)) for p in
               [r"C(..)CH",r"C(...)CH",r"C(....)CH",r"C(.{14})CH",r"C(.{15})CH"])

def second_pass(cluster_rows,g2c,seq_dict,all_results=False):
    if all_results: return cluster_rows
    stems=[r["hmm_stem"] for r in cluster_rows]; kept=[]
    for r in cluster_rows:
        stem=r["hmm_stem"]; cat=r["cat"]
        if stem=="Cyc1":
            if len({h for h in set(stems) if g2c.get(h,"") in FE_REDOX})>=2: kept.append(r)
            continue
        if re.match(r"Cyc2",stem):
            seq=seq_dict.get(r["genome"],{}).get(r["orf"],"")
            if len(seq)>=365 and count_heme(seq)>0: kept.append(r)
            continue
        if cat=="iron_gene_regulation":
            if any("regulation" in g2c.get(h,"") for h in stems): kept.append(r)
            continue
        kept.append(r)
    return kept


# ── UniOP integration ─────────────────────────────────────────────────────────
# ── UniOP integration ─────────────────────────────────────────────────────────
def _is_prodigal_faa(faa_path):
    """
    Return True if the FAA file has Prodigal-format headers with embedded
    coordinates: >orf_name # start # end # strand # ID=...
    UniOP can parse these directly with -a.
    Bakta/NCBI headers have no coordinates → need FNA input (-i).
    """
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                # Prodigal header: >ID # start # end # strand # ...
                parts = line.split("#")
                if len(parts) >= 4:
                    try:
                        int(parts[1].strip())  # start
                        int(parts[2].strip())  # end
                        return True
                    except ValueError:
                        pass
                return False
    return False


def run_uniop(faa_files, fna_dir, out_dir, uniop_path, fna_ext="fna"):
    """
    Run UniOP (https://github.com/hongsua/UniOP) on each genome.

    UniOP interface:
      python UniOP -a input.faa    # when FAA has Prodigal-format coordinate headers
      python UniOP -i input.fna    # when starting from nucleotide sequence
      python UniOP -t output_dir   # specify output directory

    Output files (fixed names, one per run):
      uniop.operon   — one operon per line: orf1,orf2,...
      uniop.pred     — pair predictions:    orf_i <TAB> orf_j <TAB> probability

    Because UniOP always writes to fixed filenames, we run each genome in its
    own temporary subdirectory to avoid collisions.

    Returns dict: genome_faa_name → {orf_id → operon_id}
    """
    uniop_dir = out_dir / "_uniop"
    uniop_dir.mkdir(exist_ok=True)

    genome_operon_map = {}

    for faa in faa_files:
        stem     = faa.stem
        work_dir = uniop_dir / stem
        work_dir.mkdir(exist_ok=True)

        operon_file = work_dir / "uniop.operon"
        pred_file   = work_dir / "uniop.pred"

        # Skip if already computed (caching)
        if operon_file.exists():
            print(f"  [INFO] UniOP cache hit for {stem}")
        else:
            # Choose input mode: FAA (Prodigal headers) or FNA
            use_faa = _is_prodigal_faa(str(faa))
            if use_faa:
                cmd = ["python", str(uniop_path), "-a", str(faa), "-t", str(work_dir)]
            else:
                # Need FNA file
                fna_path = None
                if fna_dir:
                    for ext in [fna_ext, "fna", "fasta", "fa"]:
                        candidate = Path(fna_dir) / f"{stem}.{ext}"
                        if candidate.exists():
                            fna_path = candidate
                            break
                if fna_path is None:
                    print(f"  [WARN] UniOP: no Prodigal headers in {faa.name} and no "
                          f"FNA found for {stem}. Skipping.", file=sys.stderr)
                    continue
                cmd = ["python", str(uniop_path), "-i", str(fna_path), "-t", str(work_dir)]

            r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(work_dir))
            if r.returncode != 0:
                print(f"  [WARN] UniOP failed for {stem}: {r.stderr[:200]}",
                      file=sys.stderr)
                continue

        # Parse uniop.operon (authoritative: full operons)
        # Format: one operon per line, comma-separated orf names
        # e.g.   orf_1_1,orf_1_2,orf_1_3
        orf_to_op = {}
        if operon_file.exists():
            with open(operon_file) as fh:
                for op_idx, line in enumerate(fh):
                    line = line.rstrip()
                    if not line or line.startswith("#"):
                        continue
                    op_id = f"{stem}_OP{op_idx + 1:04d}"
                    for orf in line.split(","):
                        orf = orf.strip()
                        if orf:
                            orf_to_op[orf] = op_id

        # Fall back to uniop.pred (pair probabilities) if operon file is empty
        if not orf_to_op and pred_file.exists():
            # Format: orf_i <TAB> orf_j <TAB> probability
            # Assign same operon ID to pairs with prob > 0.5 using union-find
            parent = {}
            def _find(x):
                while parent.get(x, x) != x:
                    parent[x] = parent.get(parent.get(x, x), x)
                    x = parent.get(x, x)
                return x
            def _union(a, b):
                ra, rb = _find(a), _find(b)
                if ra != rb:
                    parent[rb] = ra

            with open(pred_file) as fh:
                for line in fh:
                    line = line.rstrip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        try:
                            prob = float(parts[2])
                            if prob > 0.5:
                                _union(parts[0].strip(), parts[1].strip())
                                parent.setdefault(parts[0].strip(), parts[0].strip())
                                parent.setdefault(parts[1].strip(), parts[1].strip())
                        except (ValueError, IndexError):
                            pass

            # Assign operon IDs from connected components
            comp_ids = {}
            op_counter = 1
            for orf in list(parent.keys()):
                root = _find(orf)
                if root not in comp_ids:
                    comp_ids[root] = f"{stem}_OP{op_counter:04d}"
                    op_counter += 1
                orf_to_op[orf] = comp_ids[root]

        n_operons = len(set(orf_to_op.values()))
        n_genes   = len(orf_to_op)
        print(f"  [INFO] {stem}: {n_operons} operons, {n_genes} genes assigned")

        genome_operon_map[faa.name] = orf_to_op
        genome_operon_map[stem]     = orf_to_op

    return genome_operon_map


# ── Bakta ↔ Prodigal coordinate mapping ──────────────────────────────────────
def parse_gff_coords(gff_path, source_hint=""):
    """
    Parse a GFF/GFF3 file and return:
      dict: contig → [(start, end, strand, gene_id)]

    Works with both Prodigal GFF (ID=contig_1_5) and Bakta GFF3
    (ID=AMXMAG_00053). source_hint is only used for logging.
    """
    coords = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            try:
                start  = int(parts[3])
                end    = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            attrs  = parts[8]
            # Extract ID from attributes
            gene_id = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("ID="):
                    gene_id = attr[3:].strip()
                    break
                if attr.startswith("locus_tag=") and gene_id is None:
                    gene_id = attr[10:].strip()
            if gene_id:
                coords[contig].append((start, end, strand, gene_id))
    return dict(coords)


def build_prodigal_bakta_map(bakta_gff_dir, prodigal_gff_dir, faa_files,
                              coord_tol=3):
    """
    Build a mapping: genome → {prodigal_orf_id → bakta_gene_id}
    by matching CDS entries on (contig, start, end, strand).

    coord_tol: tolerance in bp for coordinate matching (Bakta and Prodigal
    sometimes differ by 1-3 bp at gene boundaries due to post-processing).

    For each FAA file, looks for:
      - Bakta GFF3:    bakta_gff_dir / stem.gff3
      - Prodigal GFF:  prodigal_gff_dir / stem.gff  (or .gff3, .prodigal.gff)
    """
    prodigal_to_bakta = {}   # genome_name → {prodigal_id → bakta_id}
    stats = {"matched": 0, "unmatched_prodigal": 0, "unmatched_bakta": 0}

    for faa in faa_files:
        stem   = faa.stem
        genome = faa.name

        # Find Bakta GFF3
        bakta_gff = None
        for ext in (".gff3", ".gff"):
            c = Path(bakta_gff_dir) / (stem + ext)
            if c.exists():
                bakta_gff = c
                break
        if bakta_gff is None:
            continue

        # Find Prodigal GFF
        prodigal_gff = None
        for ext in (".gff", ".gff3", ".prodigal.gff"):
            c = Path(prodigal_gff_dir) / (stem + ext)
            if c.exists():
                prodigal_gff = c
                break
        if prodigal_gff is None:
            continue

        bakta_coords    = parse_gff_coords(str(bakta_gff),    "bakta")
        prodigal_coords = parse_gff_coords(str(prodigal_gff), "prodigal")

        # Build index: contig → {(start, end, strand) → bakta_id}
        bakta_index = defaultdict(dict)
        for contig, entries in bakta_coords.items():
            for start, end, strand, gene_id in entries:
                bakta_index[contig][(start, end, strand)] = gene_id

        # Match Prodigal ORFs to Bakta genes
        orf_map = {}
        for contig, entries in prodigal_coords.items():
            b_idx = bakta_index.get(contig, {})
            for start, end, strand, prod_id in entries:
                # Exact match first
                key = (start, end, strand)
                if key in b_idx:
                    orf_map[prod_id] = b_idx[key]
                    stats["matched"] += 1
                    continue
                # Fuzzy match within coord_tol bp
                found = None
                for (bs, be, bst), bid in b_idx.items():
                    if (bst == strand and
                            abs(bs - start) <= coord_tol and
                            abs(be - end)   <= coord_tol):
                        found = bid
                        break
                if found:
                    orf_map[prod_id] = found
                    stats["matched"] += 1
                else:
                    stats["unmatched_prodigal"] += 1

        prodigal_to_bakta[genome] = orf_map

    total = stats["matched"] + stats["unmatched_prodigal"]
    if total > 0:
        pct = stats["matched"] / total * 100
        print(f"[INFO] Prodigal↔Bakta mapping: "
              f"{stats['matched']}/{total} ORFs matched ({pct:.1f}%)")
        if stats["unmatched_prodigal"] > 0:
            print(f"       {stats['unmatched_prodigal']} Prodigal ORFs had no "
                  f"matching Bakta gene (will show ORF name in Bakta column)")

    return prodigal_to_bakta


def write_operon_structure(path, final_rows, genome_operon_map,
                           prodigal_to_bakta=None):
    """
    OperonStructure_geneCall.tsv — HMM hits linked to UniOP operon predictions.
    If prodigal_to_bakta provided, adds bakta_gene_id column.
    """
    op_to_orfs = defaultdict(list)
    for r in final_rows:
        op_map = genome_operon_map.get(r["genome"], {})
        op_id  = op_map.get(r["orf"], f"singleton_{r['genome']}_{r['orf']}")
        op_to_orfs[op_id].append(r["orf"])

    use_bakta = bool(prodigal_to_bakta)
    fields = ["operon_id","genome","contig","orf","gene","category",
              "hmm_stem","bitscore","e_value","unioperon_members"]
    if use_bakta:
        fields.insert(fields.index("orf") + 1, "bakta_gene_id")

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in final_rows:
            genome = r["genome"]; orf = r["orf"]
            op_map = genome_operon_map.get(genome, {})
            op_id  = op_map.get(orf, f"singleton_{genome}_{orf}")
            members = [o for o in op_to_orfs[op_id] if o != orf]
            row = {
                "operon_id":         op_id,
                "genome":            genome,
                "contig":            r["contig"],
                "orf":               orf,
                "gene":              r["gene_name"],
                "category":          r["cat"],
                "hmm_stem":          r["hmm_stem"],
                "bitscore":          f"{r['bitscore']:.1f}",
                "e_value":           f"{r['evalue']:.2e}",
                "unioperon_members": ",".join(members) if members else "",
            }
            if use_bakta:
                row["bakta_gene_id"] = prodigal_to_bakta.get(genome, {}).get(orf, orf)
            w.writerow(row)


# ── Anvi'o output ─────────────────────────────────────────────────────────────
def write_anvio_functions(path, final_rows, prodigal_to_bakta=None):
    """
    Functions table for anvi-import-functions (tab-delimited):
      gene_callers_id    source    accession    function    e_value

    If prodigal_to_bakta is provided (Bakta external gene calls workflow),
    gene_callers_id contains Bakta gene IDs (e.g. AMXMAG_00053) which map
    directly to Anvi'o gene_callers_id when the db was built with
    anvi-gen-contigs-database --external-gene-calls from Bakta output.

    Without Bakta mapping, contains Prodigal ORF names — map to integer IDs
    via anvi-export-gene-calls before importing.
    """
    use_bakta = bool(prodigal_to_bakta)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_callers_id", "source", "accession", "function", "e_value"])
        seen = set()
        for r in final_rows:
            orf = r["orf"]; genome = r["genome"]
            if orf in seen: continue
            seen.add(orf)
            if use_bakta:
                caller = prodigal_to_bakta.get(genome, {}).get(orf, orf)
            else:
                caller = orf
            w.writerow([caller, "MetalGenie-Evo", r["hmm_stem"],
                        f"{r['gene_name']} [{r['cat']}]",
                        f"{r['evalue']:.2e}"])



# ── Coverage ──────────────────────────────────────────────────────────────────
def compute_coverage_from_bam(bam_path,out_dir,label):
    result={}; df=out_dir/f"{label}.coverage"
    r=subprocess.run(["samtools","coverage","-H","-o",str(df),bam_path],capture_output=True,text=True)
    if r.returncode==0 and df.exists():
        with open(df) as fh:
            for line in fh:
                if line.startswith("#"): continue
                p=line.rstrip().split("\t")
                if len(p)>=7:
                    try: result[p[0]]={"depth":float(p[6]),"reads":int(p[3]),"length":int(p[2])-int(p[1])+1}
                    except (ValueError,IndexError): pass
        return result
    print(f"  [WARN] samtools coverage failed for {label}, trying depth",file=sys.stderr)
    r=subprocess.run(["samtools","depth","-a",bam_path],capture_output=True,text=True)
    if r.returncode!=0:
        print(f"  [ERROR] samtools depth failed: {r.stderr[:80]}",file=sys.stderr); return {}
    sums=defaultdict(float); counts=defaultdict(int)
    for line in r.stdout.splitlines():
        p=line.split("\t")
        if len(p)>=3:
            try: sums[p[0]]+=float(p[2]); counts[p[0]]+=1
            except ValueError: pass
    for c in sums:
        n=counts[c]; result[c]={"depth":sums[c]/n if n else 0.0,"reads":0,"length":n}
    return result

def load_depth_file(path):
    result={}
    if not os.path.isfile(path):
        print(f"  [WARN] Depth file not found: {path}",file=sys.stderr); return result
    with open(path) as fh: lines=[l.rstrip() for l in fh if l.strip()]
    if not lines: return result
    cols=lines[0].split("\t")
    def _row(d=0.0,r=0,l=0): return {"depth":float(d),"reads":int(r),"length":int(l)}
    if cols[0].lower() in ("contigname","#contigname"):
        for line in lines[1:]:
            p=line.split("\t")
            if len(p)>=3 and p[0]!="contigName":
                try: result[p[0]]=_row(p[2],0,int(p[1]))
                except (ValueError,IndexError): pass
        return result
    if cols[0].lower()=="rname" or (len(cols)>=7 and cols[6].lower()=="meandepth"):
        for line in lines[1:]:
            p=line.split("\t")
            if len(p)>=7:
                try: result[p[0]]=_row(p[6],p[3],int(p[2])-int(p[1])+1)
                except (ValueError,IndexError): pass
        return result
    if cols[0].startswith("#"):
        for line in lines[1:]:
            if line.startswith("#"): continue
            p=line.split("\t")
            if len(p)>=3:
                try: result[p[0]]=_row(p[1],0,int(p[2]))
                except (ValueError,IndexError): pass
        return result
    for line in lines:
        p=line.split("\t")
        if len(p)>=2:
            try: result[p[0]]=_row(p[1])
            except ValueError: pass
    return result

def load_bams_tsv(path):
    d={}
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if not line or line.startswith("#"): continue
            p=line.split("\t")
            if len(p)>=2: d[p[0].strip()]=p[1].strip()
    return d

def build_contig_coverage(faa_files,bam_map,depth_map,out_dir):
    cov_dir=out_dir/"_coverage"; cov_dir.mkdir(exist_ok=True); gc={}
    for faa in faa_files:
        label=faa.name; stem=faa.stem
        dp=depth_map.get(label) or depth_map.get(stem)
        if dp: print(f"  [INFO] Loading depth for {stem}…"); gc[label]=load_depth_file(dp); continue
        bp=bam_map.get(label) or bam_map.get(stem)
        if bp: print(f"  [INFO] Computing coverage for {stem}…"); gc[label]=compute_coverage_from_bam(bp,cov_dir,stem)
    return gc

def compute_tpm(genome_cov,contig_lens):
    rates={}
    for c,info in genome_cov.items():
        length=contig_lens.get(c,info.get("length",1)) or 1
        rates[c]=info["depth"]/length
    total=sum(rates.values()) or 1.0
    return {c:(r/total)*1e6 for c,r in rates.items()}

# ── Writers ───────────────────────────────────────────────────────────────────
def write_summary(path,rows):
    with open(path,"w") as fh:
        fh.write("category,genome/assembly,contig,orf,gene,bitscore,"
                 "bitscore_cutoff,cluster_id,heme_c_motifs,protein_sequence\n")
        prev=None
        for r in rows:
            if prev is not None and r["cluster_id"]!=prev: fh.write("#,#,#,#,#,#,#,#,#\n")
            fh.write(f"{r['cat']},{r['genome']},{r['contig']},{r['orf']},"
                     f"{r['gene_name']},{r['bitscore']:.1f},{r['cutoff']},"
                     f"{r['cluster_id']},{r['heme_motifs']},{r['sequence']}\n")
            prev=r["cluster_id"]

def write_gene_summary(path,rows):
    with open(path,"w") as fh:
        fh.write("process,assembly,contig,orf,gene,bitscore,cluster_id\n")
        prev=None
        for r in rows:
            if prev is not None and r["cluster_id"]!=prev: fh.write("#,#,#,#,#,#,#\n")
            fh.write(f"{r['cat']},{r['genome']},{r['contig']},{r['orf']},"
                     f"{r['gene_name']},{r['bitscore']:.1f},{r['cluster_id']}\n")
            prev=r["cluster_id"]

def write_long_format(path,rows):
    fields=["category","genome","contig","orf","gene","bitscore",
            "bitscore_cutoff","cluster_id","heme_c_motifs","contig_len"]
    with open(path,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({"category":r["cat"],"genome":r["genome"],"contig":r["contig"],
                        "orf":r["orf"],"gene":r["gene_name"],
                        "bitscore":f"{r['bitscore']:.1f}","bitscore_cutoff":r["cutoff"],
                        "cluster_id":r["cluster_id"],"heme_c_motifs":r["heme_motifs"],
                        "contig_len":r.get("contig_len","")})

def write_heatmap(path,rows,all_genomes,norm_dict=None):
    all_cats=sorted({r["cat"] for r in rows})
    cm=defaultdict(lambda:defaultdict(set))
    for r in rows: cm[r["cat"]][r["genome"]].add(r["cluster_id"])
    with open(path,"w") as fh:
        fh.write("X,"+",".join(all_genomes)+"\n")
        for cat in all_cats:
            vals=[]
            for g in all_genomes:
                raw=len(cm[cat].get(g,set()))
                vals.append(f"{raw/norm_dict[g]*1000:.4f}" if norm_dict and norm_dict.get(g,0)>0 else str(raw))
            fh.write(cat+","+",".join(vals)+"\n")

def write_coverage_heatmap(path,rows,all_genomes,genome_cov,norm_coverage=False,contig_lengths=None):
    all_cats=sorted({r["cat"] for r in rows})
    cm=defaultdict(lambda:defaultdict(float))
    for r in rows:
        info=genome_cov.get(r["genome"],{}).get(r["contig"],{})
        if not info: continue
        if norm_coverage and contig_lengths:
            tpm=compute_tpm(genome_cov.get(r["genome"],{}),contig_lengths.get(r["genome"],{}))
            val=tpm.get(r["contig"],0.0)
        else: val=info.get("depth",0.0)
        cm[r["cat"]][r["genome"]]+=val
    label="TPM" if norm_coverage else "mean_depth_sum"
    with open(path,"w") as fh:
        fh.write(f"# coverage metric: {label}\n")
        fh.write("X,"+",".join(all_genomes)+"\n")
        for cat in all_cats:
            fh.write(cat+","+",".join(f"{cm[cat].get(g,0.0):.4f}" for g in all_genomes)+"\n")

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    p=argparse.ArgumentParser(prog="MetalGenie-Evo",
        description="HMM-based annotation of iron cycling and metal resistance genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ig=p.add_mutually_exclusive_group(required=True)
    ig.add_argument("--faa_dir",help="Directory of ORF .faa files (Prodigal output)")
    ig.add_argument("--fna_dir",help="Directory of nucleotide assemblies — Prodigal run internally")
    p.add_argument("--faa_ext",default="faa")
    p.add_argument("--fna_ext",default="fna")
    p.add_argument("--meta",action="store_true",help="Prodigal metagenomic mode (-p meta)")
    p.add_argument("--gff_dir",help="Prodigal GFF files for bp clustering")
    p.add_argument("--hmm_dir",required=True)
    p.add_argument("--out",default="metalgenie_evo_out")
    p.add_argument("--min_contig_len",type=int,default=0,
                   help="Skip ORFs on contigs shorter than this (bp). 0=no filter")
    p.add_argument("--max_gap",type=int,default=5,help="Max ORF-index gap (index mode)")
    p.add_argument("--max_bp_gap",type=int,default=5000,help="Max bp gap (GFF mode)")
    p.add_argument("--strand_aware",action="store_true")
    p.add_argument("--relaxed_operons",action="store_true",
                   help="Halve operon min_genes for contigs < --relaxed_threshold")
    p.add_argument("--relaxed_threshold",type=int,default=10000,
                   help="Contig length (bp) threshold for relaxed operon rules")
    p.add_argument("--all_results",action="store_true")
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--hmm_threads",type=int,default=1)
    p.add_argument("--norm",action="store_true",help="Normalise gene-count heatmap")
    p.add_argument("--bam",help="Single BAM file (requires samtools >=1.10)")
    p.add_argument("--bams",help="TSV: genome<TAB>bam_path")
    p.add_argument("--depth",help="Pre-computed depth file (jgi/BBMap/samtools/plain)")
    p.add_argument("--depths",help="TSV: genome<TAB>depth_path")
    p.add_argument("--norm_coverage",action="store_true",help="TPM-normalise coverage heatmap")
    p.add_argument("--keep_tblout",action="store_true")
    # ── Operon prediction (UniOP) ─────────────────────────────────────────────
    p.add_argument("--operon_prediction",action="store_true",
                   help="Run UniOP operon prediction and produce "
                        "OperonStructure_geneCall.tsv. Requires --fna_dir "
                        "(or nucleotide files) and UniOP installed.")
    p.add_argument("--uniop_path",default="uniop",
                   help="Path to UniOP executable (default: 'uniop', assumed in PATH)")
    # ── Anvi'o output ─────────────────────────────────────────────────────────
    p.add_argument("--anvio",action="store_true",
                   help="Write MetalGenie-Evo-anvio-functions.tsv, compatible with "
                        "anvi-import-functions. gene_callers_id column contains ORF "
                        "names — map to integer IDs from your Anvi'o contigs database "
                        "before importing (see README).")
    p.add_argument("--bakta_gff_dir",
                   help="Directory of Bakta GFF3 files. When provided together with "
                        "--fna_dir, MetalGenie-Evo runs Prodigal internally for HMM "
                        "search and UniOP, then maps Prodigal ORF names back to Bakta "
                        "gene IDs via coordinate matching. The --anvio output will use "
                        "Bakta IDs directly, compatible with Anvi'o databases built "
                        "from Bakta external gene calls.")
    args=p.parse_args()

    out_dir=Path(args.out); out_dir.mkdir(parents=True,exist_ok=True)
    tblout_dir=out_dir/"_tblout_cache"; tblout_dir.mkdir(exist_ok=True)
    gff_dir_path=None

    if args.fna_dir:
        fna_dir=Path(args.fna_dir)
        fna_files=sorted(fna_dir.glob(f"*.{args.fna_ext}"))
        if not fna_files: sys.exit(f"[ERROR] No .{args.fna_ext} in {fna_dir}")
        faa_dir,gff_dir_path=run_prodigal(fna_files,out_dir,meta_mode=args.meta,threads=args.threads)
        faa_ext="faa"
    else:
        faa_dir=Path(args.faa_dir); faa_ext=args.faa_ext
        if args.gff_dir: gff_dir_path=Path(args.gff_dir)

    faa_files=sorted(faa_dir.glob(f"*.{faa_ext}"))
    if not faa_files: sys.exit(f"[ERROR] No .{faa_ext} in {faa_dir}")
    print(f"[INFO] {len(faa_files)} genome/bin FAA files")

    gene_map=read_map(str(Path(args.hmm_dir)/"MetalGenie-map.txt"))
    if not gene_map: gene_map=read_map(str(Path(args.hmm_dir)/"FeGenie-map.txt"))
    cutoffs=read_cutoffs(str(Path(args.hmm_dir)/"HMM-bitcutoffs.txt"))
    cat_hmms=defaultdict(list); h2c={}
    for entry in sorted(Path(args.hmm_dir).iterdir()):
        if entry.is_dir() and not entry.name.startswith("."):
            for hf in sorted(entry.glob("*.hmm")):
                cat_hmms[entry.name].append((hf.stem,hf)); h2c[hf.stem]=entry.name
    if not cat_hmms: sys.exit(f"[ERROR] No HMM dirs in {args.hmm_dir}")
    total_hmms=sum(len(v) for v in cat_hmms.values())
    print(f"[INFO] {total_hmms} HMMs across {len(cat_hmms)} categories")
    operon_rules,report_all_pats,json_mode=load_operon_rules(args.hmm_dir)

    contig_lengths={}
    if args.min_contig_len>0 or args.relaxed_operons or args.norm_coverage:
        print("[INFO] Loading contig lengths…")
        contig_lengths=build_contig_length_dict(
            faa_files,gff_dir=gff_dir_path,
            fna_dir=Path(args.fna_dir) if args.fna_dir else None,
            fna_ext=args.fna_ext if args.fna_dir else "fna")
        n_cl=sum(1 for f in faa_files if f.name in contig_lengths)
        print(f"       {n_cl}/{len(faa_files)} genomes have contig length data")

    genome_coords={}
    if gff_dir_path:
        print(f"[INFO] Loading GFF from {gff_dir_path}…")
        genome_coords=load_gff_dir(gff_dir_path,faa_files)
        print(f"       {sum(1 for f in faa_files if f.name in genome_coords)}/{len(faa_files)} with GFF")
    else:
        print("[INFO] No GFF: using ORF-index clustering (FeGenie-compatible)")

    print("[INFO] Loading protein sequences…")
    seq_dict={faa.name:read_fasta(str(faa)) for faa in faa_files}
    print(f"[INFO] Launching hmmsearch ({args.threads} threads, {args.hmm_threads} per job)…")
    run_all_hmmsearches(faa_files,cat_hmms,cutoffs,tblout_dir,args.threads,args.hmm_threads)
    print("[INFO] Collecting best HMM hits…")
    best_hit=collect_best_hits(faa_files,cat_hmms,tblout_dir)

    print("[INFO] Clustering and filtering…")
    if args.min_contig_len>0: print(f"       Skipping contigs < {args.min_contig_len} bp")
    if args.relaxed_operons: print(f"       Relaxed thresholds for contigs < {args.relaxed_threshold} bp")
    cluster_id=0; final_rows=[]; n_filt=0

    for faa in faa_files:
        genome=faa.name; orf_hits=best_hit.get(genome,{})
        if not orf_hits: continue
        orf_coords=genome_coords.get(genome,{}); clen=contig_lengths.get(genome,{})
        if args.min_contig_len>0 and clen:
            before=len(orf_hits)
            orf_hits={o:h for o,h in orf_hits.items()
                      if clen.get(_orf_to_contig(o),args.min_contig_len+1)>=args.min_contig_len}
            n_filt+=before-len(orf_hits)
        raw_clusters=build_clusters(genome,orf_hits,orf_coords,
                                    args.max_gap,args.max_bp_gap,args.strand_aware)
        for orf_group in raw_clusters:
            cluster_rows=[]
            for orf in orf_group:
                hit=orf_hits.get(orf)
                if hit is None: continue
                contig=_orf_to_contig(orf)
                cluster_rows.append({"cat":hit["cat"],"genome":genome,"contig":contig,
                    "orf":orf,"hmm_stem":hit["hmm_stem"],"bitscore":hit["bitscore"],
                    "cutoff":cutoffs.get(hit["hmm_stem"],0),"evalue":hit["evalue"],
                    "cluster_id":cluster_id,"contig_len":clen.get(contig,0)})
            if not cluster_rows: cluster_id+=1; continue
            min_clen=(min(r["contig_len"] for r in cluster_rows if r["contig_len"]>0)
                      if any(r["contig_len"]>0 for r in cluster_rows) else None)
            rel_thr=args.relaxed_threshold if args.relaxed_operons else 0
            if json_mode:
                filtered=filter_cluster_json(cluster_rows,operon_rules,report_all_pats,
                                              args.all_results,contig_len=min_clen,
                                              relaxed_threshold=rel_thr)
            else:
                filtered=filter_cluster_fegenie(cluster_rows,report_all_pats,
                                                 args.all_results)
            filtered=second_pass(filtered,h2c,seq_dict,args.all_results)
            for r in filtered:
                r["gene_name"]=gene_map.get(r["hmm_stem"],r["hmm_stem"])
                r["sequence"]=seq_dict.get(genome,{}).get(r["orf"],"")
                r["heme_motifs"]=count_heme(r["sequence"])
                final_rows.append(r)
            cluster_id+=1

    if n_filt: print(f"[INFO] {n_filt} ORFs removed (short contig)")
    final_rows.sort(key=lambda r:(r["cluster_id"],r["orf"]))
    norm_dict={faa.name:len(seq_dict.get(faa.name,{})) for faa in faa_files} if args.norm else None
    all_genomes=sorted(f.name for f in faa_files)

    for path,fn in [(out_dir/"MetalGenie-Evo-summary.csv",write_summary),
                    (out_dir/"MetalGenie-Evo-geneSummary-clusters.csv",write_gene_summary),
                    (out_dir/"MetalGenie-Evo-results-long.tsv",write_long_format)]:
        print(f"[INFO] Writing {path.name}…"); fn(str(path),final_rows)
    print(f"[INFO] Writing MetalGenie-Evo-heatmap-data.csv…")
    write_heatmap(str(out_dir/"MetalGenie-Evo-heatmap-data.csv"),final_rows,all_genomes,norm_dict)

    bam_map={}; depth_map={}
    if args.bams: bam_map=load_bams_tsv(args.bams)
    elif args.bam:
        for faa in faa_files: bam_map[faa.name]=args.bam; bam_map[faa.stem]=args.bam
    if args.depths: depth_map=load_bams_tsv(args.depths)
    elif args.depth:
        for faa in faa_files: depth_map[faa.name]=args.depth; depth_map[faa.stem]=args.depth
    if bam_map or depth_map:
        print("[INFO] Computing coverage…")
        gc=build_contig_coverage(faa_files,bam_map,depth_map,out_dir)
        if gc:
            print("[INFO] Writing MetalGenie-Evo-coverage-heatmap.csv…")
            write_coverage_heatmap(str(out_dir/"MetalGenie-Evo-coverage-heatmap.csv"),
                                   final_rows,all_genomes,gc,
                                   norm_coverage=args.norm_coverage,
                                   contig_lengths=contig_lengths if args.norm_coverage else None)
        else: print("[WARN] No coverage data loaded.",file=sys.stderr)

    if not args.keep_tblout: shutil.rmtree(tblout_dir,ignore_errors=True)
    else: print(f"[INFO] tblout cache at {tblout_dir}/")

    # ── UniOP operon prediction (optional) ────────────────────────────────────
    genome_operon_map  = {}
    prodigal_to_bakta  = {}

    if args.operon_prediction:
        fna_dir_for_uniop = Path(args.fna_dir) if args.fna_dir else None
        print(f"[INFO] Running UniOP on {len(faa_files)} genomes…")
        print(f"       UniOP path: {args.uniop_path}")
        if fna_dir_for_uniop is None:
            print("       No --fna_dir: will use FAA files directly "
                  "(requires Prodigal-format headers)")
        genome_operon_map = run_uniop(
            faa_files,
            fna_dir    = fna_dir_for_uniop,
            out_dir    = out_dir,
            uniop_path = args.uniop_path,
            fna_ext    = args.fna_ext if args.fna_dir else "fna",
        )
        if genome_operon_map and args.bakta_gff_dir and gff_dir_path:
            print("[INFO] Building Prodigal↔Bakta coordinate map…")
            prodigal_to_bakta = build_prodigal_bakta_map(
                args.bakta_gff_dir,
                str(gff_dir_path),
                faa_files)
            op_path = out_dir / "MetalGenie-Evo-OperonStructure.tsv"
            print(f"[INFO] Writing {op_path.name}…")
            write_operon_structure(str(op_path), final_rows, genome_operon_map,
                                   prodigal_to_bakta=prodigal_to_bakta)
        else:
            print("[WARN] UniOP produced no predictions.", file=sys.stderr)

    # ── Bakta ↔ Prodigal mapping (without UniOP) ──────────────────────────────
    # Build the map even if --operon_prediction was not requested,
    # so --anvio can use Bakta IDs directly.
    if args.bakta_gff_dir and not args.operon_prediction:
        if gff_dir_path and args.fna_dir:
            print("[INFO] Building Prodigal↔Bakta coordinate map…")
            prodigal_to_bakta = build_prodigal_bakta_map(
                args.bakta_gff_dir,
                str(gff_dir_path),
                faa_files)

    # ── Anvi'o functions output (optional) ───────────────────────────────────
    if args.anvio:
        anvio_path = out_dir / "MetalGenie-Evo-anvio-functions.tsv"
        print(f"[INFO] Writing {anvio_path.name}…")
        write_anvio_functions(str(anvio_path), final_rows,
                              prodigal_to_bakta=prodigal_to_bakta)
        id_note = ("Bakta gene IDs" if prodigal_to_bakta
                   else "Prodigal ORF names — map to int IDs before import")
        print(f"       gene_callers_id: {id_note}")
        print(f"       Import with: anvi-import-functions -c CONTIGS.db "
              f"-i {anvio_path.name} -p MetalGenie-Evo")

    n_hit=len({r["genome"] for r in final_rows}); cc=defaultdict(int)
    for r in final_rows: cc[r["cat"]]+=1
    print(f"\n{'─'*60}\n  MetalGenie-Evo  —  run complete")
    print(f"  {len(final_rows)} ORFs in {n_hit}/{len(faa_files)} genomes")
    print(f"\n  Hits per category:")
    for cat,n in sorted(cc.items()): print(f"    {n:5d}  {cat}")
    print(f"\n  Outputs  →  {out_dir}/\n{'─'*60}")

if __name__=="__main__": main()