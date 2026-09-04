#!/usr/bin/env bash
## One-shot local run for the sleuthing pass. From the biolai-peers folder:
##   bash run_local.sh
## Needs R with jsonlite and curl (installed below if missing) and open web
## access. Takes roughly 15 to 30 minutes, almost all of it the OpenAlex pull.
set -euo pipefail
cd "$(dirname "$0")"
export BIOLAI_KEYS="${BIOLAI_KEYS:-tamu,purdue,ncstate,lsu,vt,auburn,kstate,texastech}"
export BIOLAI_MAILTO="${BIOLAI_MAILTO:-coleoguy@gmail.com}"
Rscript -e 'for (p in c("jsonlite","curl")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")'
echo "== 08 directory counts (faculty and grad student pages)"; Rscript R/08_directory_counts.R
echo "== 01 IPEDS completions and Pell";                          Rscript R/01_ipeds.R
echo "== 03 OpenAlex research-active faculty";                    Rscript R/03_openalex_faculty.R
echo "== 02 ProQuest strings + OpenAlex dissertations";           Rscript R/02_proquest.R
echo
echo "Done. Commit and push the outputs so the remote session can finish the report:"
echo "  git add data/raw data/derived && git commit -m 'Local pulls: directories, IPEDS, OpenAlex' && git push"
