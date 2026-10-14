#!/bin/bash
# parallel, resumable ranged download: pdl.sh URL OUT [N]
url=$1; out=$2; n=${3:-8}
size=$(curl -sIL "$url" | grep -i '^content-length' | tail -1 | tr -dc 0-9)
echo "size=$size"
chunk=$(( (size + n - 1) / n ))
for i in $(seq 0 $((n-1))); do
  s=$((i*chunk)); e=$(( s+chunk-1 )); [ $e -ge $size ] && e=$((size-1))
  f="$out.p$i"; need=$((e-s+1))
  ( while true; do
      have=$(stat -c %s "$f" 2>/dev/null || echo 0)
      [ "$have" -ge "$need" ] && break
      curl -sSfL -r $((s+have))-$e "$url" >> "$f" || sleep 3
    done ) &
done
wait
cat $(for i in $(seq 0 $((n-1))); do echo "$out.p$i"; done) > "$out" && rm -f $out.p*
echo "done $out $(stat -c %s $out)"
