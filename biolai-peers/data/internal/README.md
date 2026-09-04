# data/internal (not committed)

Files the pipeline reads from here, all optional:

| file | columns | produced by |
|---|---|---|
| `faculty_handcount.csv` | inst_key, n_tt_faculty, directory_url, date_counted | hand count of tenure-track faculty on five department directory pages (tamu, purdue, florida, colostate, kstate); calibrates 03 |
| `coordinator_replies.csv` | inst_key, year, phd_granted, phd_enrolled, reply_date, source_email | grad coordinator email replies; overrides ProQuest in 05 and 06 |
| `placement.csv` | year, n_grads, n_phd | TAMU Biology internal rosters: BS graduates and how many entered PhD programs |

Templates with the right headers are in this directory as `*_template.csv`.
