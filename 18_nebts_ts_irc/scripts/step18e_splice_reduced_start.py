#!/usr/bin/env python3
# step18e splice: build a reduced-model START geometry = common frozen bulk (ts_start.pdb)
# with ONLY the 102 active-atom coordinates overwritten by a target state (reactant/product),
# so the subsequent reduced-region Opt relaxes to that state within the identical model as
# the step18c OptTS. Verbatim column copy preserves PDB fixed-width fields; a per-atom
# name/resname check aborts if the two files' atom ordering ever disagrees.
# Pure-Python (math only), pre-3.8 compatible. No numpy.
import sys, math
SUB0 = 6207
COORD = slice(30, 54); NAME = slice(12, 16); RESN = slice(17, 20)
def is_atom(l): return l.startswith("ATOM  ") or l.startswith("HETATM")
def load_active(p):
    with open(p) as fh: return set(int(t) for t in fh.read().split())
def collect(p, aset):
    coords={}; names={}; ai=-1; n=0
    with open(p) as fh:
        for l in fh:
            if is_atom(l):
                ai+=1; n+=1
                if ai in aset: coords[ai]=l[COORD]; names[ai]=(l[NAME],l[RESN])
    return coords, names, n
def main():
    active_path,bulk_path,target_path,out_path=sys.argv[1:5]
    aset=load_active(active_path); print("active atoms: %d"%len(aset))
    tc,tn,tnat=collect(target_path,aset); print("target atoms: %d"%tnat)
    if aset-set(tc): raise SystemExit("FAIL target missing active indices")
    out_lines=[]; ai=-1; bnat=0; changed=0; mism=[]
    with open(bulk_path) as fh:
        for l in fh:
            if is_atom(l):
                ai+=1; bnat+=1
                if ai in aset:
                    if (l[NAME],l[RESN])!=tn[ai]:
                        mism.append((ai,l[NAME].strip(),tn[ai][0].strip()))
                    l=l[:30]+tc[ai]+l[54:]; changed+=1
            out_lines.append(l)
    print("bulk atoms:   %d"%bnat)
    if bnat!=tnat: raise SystemExit("FAIL atom-count mismatch")
    if mism: raise SystemExit("FAIL %d active-atom identity mismatches (ordering differs)"%len(mism))
    if changed!=len(aset): raise SystemExit("FAIL changed != active count")
    with open(out_path,"w") as fh: fh.writelines(out_lines)
    print("spliced atoms: %d\nwrote: %s"%(changed,out_path))
    def xyz(idx):
        a=-1
        for l in out_lines:
            if is_atom(l):
                a+=1
                if a==idx: return (float(l[30:38]),float(l[38:46]),float(l[46:54]))
    if bnat>SUB0+12:
        def d(a,b):
            pa,pb=xyz(a),xyz(b); return math.sqrt(sum((pa[k]-pb[k])**2 for k in range(3)))
        print("substrate START: break(O3-C4)=%.3f form(C1-C6)=%.3f"%(d(SUB0+7,SUB0+8),d(SUB0,SUB0+12)))
if __name__=="__main__": main()
