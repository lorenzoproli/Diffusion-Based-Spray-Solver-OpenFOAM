#!/usr/bin/env bash

fnum() {
    awk -v x="$1" 'BEGIN {
        if (x == "" || x ~ /nan|inf/i) { printf "%s", x; exit }
        if (x != 0 && (x < 1e-3 || x > 1e3)) printf "%.3e", x;
        else printf "%.6f", x;
    }'
}

print_sep() {
    printf '%*s\n' "${COLUMNS:-120}" '' | tr ' ' '-'
}

print_last_data_row() {
    local file="$1"
    [[ -f "$file" ]] || return 1
    awk 'NF && $1 !~ /^#/ { last=$0 } END { if (last!="") print last; else exit 1 }' "$file"
}

echo "Monitor breakup/debug - $(date)"
print_sep

# ============================================================
# DRAG SUMMARY
# ============================================================
echo "[DRAG SUMMARY]"
if row=$(print_last_data_row dragDebug.log 2>/dev/null); then
    tc=$(echo "$row" | awk -F, '{print $2}')
    tcb=$(echo "$row" | awk -F, '{print $5}')
    prog=$(awk -v a="$tc" -v b="$tcb" 'BEGIN{ if (b>0) printf "%.2f", 100*a/b; else print "nan"}')
    printf "%-12s %-16s %-12s\n" "tc" "tColumnBreakup" "progress%"
    printf "%-12s %-16s %-12s\n" "$(fnum "$tc")" "$(fnum "$tcb")" "$prog"
else
    echo "dragDebug.log non disponibile"
fi
print_sep

# ============================================================
# FILES
# ============================================================
echo "[FILES]"
ls -lh breakupDebug_*.log dragDebug.log 2>/dev/null | awk '{printf "%-28s %8s %s %s %s\n",$9,$5,$6,$7,$8}'
print_sep

# ============================================================
# WAVE - last 3
# ============================================================
echo "[WAVE - last 3 events]"
if [[ -f breakupDebug_wave.log ]]; then
    awk -F, '
    BEGIN {
        printf "%-10s %-10s %-10s %-12s %-12s %-10s %-10s %-12s %-12s %-10s %-10s\n",
               "tc","d","dOld","dChildTgt","mStrip","ms","WeGas","lambdaKH","tauKH","Uchild","Uparent"
    }
    NR>1 {rows[n%3]=$0; n++}
    END {
        start=(n>3?n-3:0);
        for(i=start;i<n;i++){
            split(rows[i%3],a,",");
            printf "%-10.4e %-10.4e %-10.4e %-12.4e %-12.4e %-10.4e %-10.4f %-12.4e %-12.4e %-10.4f %-10.4f\n",
                   a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11];
        }
    }' breakupDebug_wave.log
else
    echo "breakupDebug_wave.log non disponibile"
fi
print_sep

# ============================================================
# STDPE - last 5
# ============================================================
echo "[STDPE - last 5 events]"
if [[ -f breakupDebug_stdpe.log ]]; then
    awk -F, '
    BEGIN {
        printf "%-10s %-10s %-10s %-10s %-10s %-10s %-10s %-8s\n",
               "event","type","tc","d","Urmag","We","WeCrit","allow"
    }
    NR>1 {rows[n%5]=$0; n++}
    END {
        start=(n>5?n-5:0);
        for(i=start;i<n;i++){
            split(rows[i%5],a,",");
            printf "%-10s %-10s %-10.3e %-10.3e %-10.3f %-10.3f %-10.3f %-8s\n",
                   a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8];
        }
    }' breakupDebug_stdpe.log
else
    echo "breakupDebug_stdpe.log non disponibile"
fi
print_sep

# ============================================================
# STAGE2 - last 5
# ============================================================
echo "[STAGE2 - last 5 events]"
if [[ -f breakupDebug_stage2.log ]]; then
    awk -F, '
    BEGIN {
        printf "%-10s %-10s %-10s %-10s %-12s %-12s %-10s\n",
               "tc","tCb","rhoc","UgRef","Dparent0","UrelPE0","Urmag"
    }
    NR>1 && NF>=7 {rows[n%5]=$0; n++}
    END {
        start=(n>5?n-5:0);
        for(i=start;i<n;i++){
            split(rows[i%5],a,",");
            printf "%-10.3e %-10.3e %-10.4f %-10.3f %-12.4e %-12.4f %-10.4f\n",
                   a[1],a[2],a[3],a[4],a[5],a[6],a[7];
        }
    }' breakupDebug_stage2.log
else
    echo "breakupDebug_stage2.log non disponibile"
fi
print_sep

# ============================================================
# BREAKUP - last 5
# ============================================================
echo "[BREAKUP - last 5 events]"
if [[ -f breakupDebug_breakup.log ]]; then
    awk -F, '
    BEGIN {
        printf "%-8s %-10s %-12s %-10s %-10s %-10s %-10s %-10s %-8s %-10s\n",
               "mode","tc","Dparent0","UrelPE0","tElapsed","tDef","tb","dSMD","FLig","dFrag0"
    }
    NR>1 {rows[n%5]=$0; n++}
    END {
        start=(n>5?n-5:0);
        for(i=start;i<n;i++){
            split(rows[i%5],a,",");
            printf "%-8s %-10.3e %-12.4e %-10.4f %-10.3e %-10.3e %-10.3e %-10.4e %-8.3f %-10.4e\n",
                   a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[10],a[11],a[12];
        }
    }' breakupDebug_breakup.log
else
    echo "breakupDebug_breakup.log non disponibile"
fi
print_sep
