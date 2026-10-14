# NOTE ON INPUTS: this script reads its inputs from ../../02_derived_data/, which is
# NOT distributed with this repository (see README section 1). Rebuild it by running
# the upstream step, or obtain the tables from the manuscript's supplementary files.
"""GO/KEGG over-representation for the 169-gene intersection via Enrichr (multipart)."""
import urllib.request, json, time, sys, uuid
import os
_HERE = os.path.dirname(os.path.abspath(__file__))

sys.stdout.reconfigure(encoding="utf-8")

GENEFILE = os.path.join(_HERE, "..", "..", "02_derived_data", "169_genes.csv")
genes = [l.strip() for l in open(GENEFILE, encoding="utf-8-sig") if l.strip()]
print("input genes:", len(genes))

BASE = "https://maayanlab.cloud/Enrichr"


def post_multipart(url, fields):
    b = "----wb" + uuid.uuid4().hex
    body = b""
    for k, v in fields.items():
        body += f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode()
    body += f"--{b}--\r\n".encode()
    req = urllib.request.Request(
        url, data=body,
        headers={"User-Agent": "Mozilla/5.0",
                 "Content-Type": f"multipart/form-data; boundary={b}"})
    return urllib.request.urlopen(req, timeout=90)


res = json.load(post_multipart(BASE + "/addList",
                               {"list": "\n".join(genes), "description": "intersection169"}))
uid = res["userListId"]
print("userListId:", uid)

LIBS = ["GO_Biological_Process_2023", "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023", "KEGG_2021_Human"]
out = {}
for lib in LIBS:
    url = f"{BASE}/enrich?userListId={uid}&backgroundType={lib}"
    rows = []
    for _ in range(4):
        try:
            rows = json.load(urllib.request.urlopen(url, timeout=90)).get(lib, [])
            break
        except Exception as e:
            print("retry", lib, e); time.sleep(4)
    out[lib] = rows
    sig = [r for r in rows if r[6] is not None and r[6] < 0.05]
    print(f"\n########## {lib}  total={len(rows)}  adjP<0.05: {len(sig)} ##########")
    for r in sig[:20]:
        g = r[5] or []
        print(f"  adjP={r[6]:.2e}  n={len(g):>3}  {r[1]}")

json.dump(out, open(os.path.join(_HERE, "_enrich169_raw.json"),
                    "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved _enrich169_raw.json")
