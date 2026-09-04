# BiolAI peer departments: context

## Goal
Pick 8 to 10 structural peers of TAMU Biology (Arts and Sciences), lock a 2019 to 2025 baseline, and run a department-year difference-in-differences (plus synthetic control) on PhDs, grant dollars, publications, trainee output, and undergrad PhD placement, re-pulled annually with the same scripts.

## Status (2026-09-04)
- Pipeline written and offline-tested end to end on synthetic fixtures (scripts 01, 02, 05, 06, 07 executed; 03 and 04 parse-checked only, their API calls need an unblocked network). Base R plus jsonlite and curl.
- Seed pool roster (`data/seed_pool.csv`, 25 candidates plus TAMU) filled with colleges, unit boundaries, IPEDS unitids, ROR ids, regexes for affiliation and ProQuest matching, contamination verdicts, coordinator contacts, and a `verification` column.
- No quantitative data pulled yet: this session's egress policy blocked nces.ed.gov, api.openalex.org, api.nsf.gov, ncses.nsf.gov, and every .edu host. Web search budget ran out at 200 calls, so 9 of 25 contamination checks are UNKNOWN.
- `data/peers.csv` does not exist yet; `peers-v1` is not tagged. Do not tag until 05 has run on real data and the UNKNOWN verdicts are resolved.
- Memo (`memo/peer_justification_memo.docx`) is a structural draft: filters and unit boundaries only, distance table to be inserted after 05 runs.

## Bottlenecks
1. Network: run 01, 03, 04 from a machine with open egress (laptop or HPRC). Expect 03 to take ~1 to 2 h for 26 institutions at OpenAlex polite-pool rates.
2. ProQuest has no API: 02 prints one search string per institution; export each result list as CSV into `data/raw/proquest/`.
3. HERD is interactive tables only: download biological and biomedical sciences R&D by institution for FY2019 to FY2024 into `data/raw/herd/herd_bio_<FY>.csv` (columns unitid or institution, amount_k).
4. Faculty hand-count for calibration (tamu, purdue, florida, colostate, kstate) not done: directory URLs are in `data/internal/faculty_handcount_template.csv`. Ten minutes of clicking.
5. Contamination check still UNKNOWN for: okstate, missouri, lsu, auburn, msstate, oregonstate, cornell, minnesota, wisconsin. Six searches each (NRT, NRT AI, AI biology fellowship, data science biology scholars, HHMI, computational biology T32).

## Next actions
1. Run `Rscript R/01_ipeds.R` and `Rscript R/03_openalex_faculty.R` on a laptop; commit `data/raw/` and `data/derived/`.
2. Finish the 9 UNKNOWN contamination checks; update `ai_bio_verdict` and `ai_bio_note` in `data/seed_pool.csv`.
3. Hand-count faculty at the five calibration departments; fill `data/internal/faculty_handcount.csv`; re-run 03.
4. Pull HERD tables; run ProQuest searches; run 02, 04, 05. Inspect `data/derived/mahalanobis_table.csv` and `grants_pi_match_diagnostic.csv` (PI match rate at TAMU should exceed 50% of NSF BIO awards, otherwise the OpenAlex attribution is too weak and 04 should fall back to directorate plus NIH dept_type).
5. Lock: `git add data/peers.csv && git commit && git tag -a peers-v1 -m "peer set locked <date>"`.
6. Send the coordinator email (`email/grad_coordinator_email.md`; Gmail draft created 2026-09-04) to the locked peers; enter replies in `data/internal/coordinator_replies.csv`.
7. Insert the distance table and covariate means into the memo; send to the Dean's office.
8. Annual: re-run 03, 04, 06, 07 each September; figures land in `figures/`.

## Decisions
- Project lives in `biolai-peers/` inside the evobiR repo on the designated branch, excluded from the package build via `.Rbuildignore`. Putting analysis scripts in the package `R/` would break R CMD check.
- Unit of comparison is defined by content (general biology: cell/molecular, EEB, physiology, micro, neuro), not by college label. Entomology and biochemistry units are excluded everywhere even when they sit inside the A&S college (Illinois SIB entomology, Illinois MCB biochemistry, CSU and K-State biochemistry, Penn State BMB, Minnesota BMBB, Iowa State BBMB, Georgia BCMB). Tennessee BCMB and LSU Biological Sciences cannot be split from biochemistry; both flagged.
- IPEDS covariate uses CIP 26.0101 + 26.0199 (first major, bachelor's). 26.0102 Biomedical Sciences is inside 26.01 but at TAMU it is the vet college BIMS major (~1,000 degrees/yr), which the design excludes; it is reported as its own column and kept out of the covariate. Split-biology peers (UC Davis, Minnesota, Illinois, Ohio State, Georgia) confer general-biology degrees under 26.05, 26.07, 26.13, so the narrow count undercounts them; 05 reports an all-26 sensitivity distance.
- Penn State reports all campuses under one unitid since 2019-20, so its completions include Commonwealth Campuses; it also has no vet school. Both noted in the roster.
- Florida's Biology BS is co-administered with CALS and Cornell's Biological Sciences major spans A&S and CALS; IPEDS overcounts the A&S unit at both.
- Contamination verdicts: IN-DEPARTMENT for Michigan State (NRT-IMPACTS in Plant Biology) and Georgia (Genetics T32 with quantitative emphasis, Odum IDEAS NRT). Florida is UNIVERSITY-LEVEL ONLY but has the highest institution-wide AI exposure in the pool (AI Scholars, AI across the curriculum, 100+ AI hires); recommend a sensitivity run without UF. Tennessee flagged for NIMBioS and Bredesen data-science exposure.
- Hard filter 1 as written ("land-grant or state flagship") drops Texas Tech (neither) and Cornell (private). Penn State fails the vet filter. All three still appear in the distance table.
- Stipend study frame (public AAU plus SEC) overlaps the seed pool at 15 institutions; SEC members without vet schools or ag colleges (Alabama, Arkansas, Kentucky, Ole Miss, Oklahoma, South Carolina) were not added. For Missouri the stipend study's primary row was the CVM Biomedical Sciences PhD; the correct analog here is the Division of Biological Sciences (its secondary row).
- Peer count rule in 05: take 10 when the 9th and 10th nearest eligible are within 10% of the 8th's distance, else 8.
- Faculty proxy: OpenAlex authors with 3+ works 2021 to 2025 whose raw affiliation string matches the department regex; rescaled by the median hand-count ratio across five departments. Cannot separate faculty from productive postdocs; the calibration spread is the honest error bar.
- Grants: NSF and NIH awards attributed to the department by PI-to-OpenAlex-author name match; NIH dept_type and NSF BIO directorate kept as alternative attributions. Outcome is the 3-year rolling sum.
- Inference in 07: CR1 cluster-robust SEs by institution plus permutation (in-space placebo) p-values, because one treated unit makes the cluster count small. Synthetic control weights by softmax-parameterised optim (base R), placebo RMSPE-ratio test.
- Treatment year set to 2026 (`treat_year` in 00_config.R); pre-period 2019 to 2025; event study omits 2025.

## People
- Grad coordinators for every candidate: `data/seed_pool.csv` (coordinator_name, coordinator_email). Missing: Oregon State, Cornell, Minnesota, Wisconsin.
- Dean's office: memo recipient.

## Deadlines
- None fixed. Suggested: lock peers-v1 before the fall 2026 semester data close so the 2026 outcome year is the first post-treatment observation.
