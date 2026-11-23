# Herb-Drug Interaction Checker (Streamlit)

Mechanism-transparent pharmacokinetic interaction screening for South Asian
medicinal plants and prescription drugs. Research screening tool, not clinical advice.

## Files

    app.py            Streamlit entry point; serves index.html with engine.js and data.js inlined
    index.html        interface
    engine.js         prediction functions (no I/O)
    data.js           curated dataset with per-parameter provenance
    requirements.txt  pinned Streamlit version
    LICENSE           MIT, code only

## Run locally

    pip install -r requirements.txt
    streamlit run app.py

## Deploy on Streamlit Community Cloud

1. Push this folder to a GitHub repository (the folder root is the repo root).
2. Go to share.streamlit.io, choose "Create app", and select the repo, branch `main`,
   and main file path `app.py`.
3. Deploy. Community Cloud installs `requirements.txt` automatically.

## Licence and data

Code is MIT. Seed data is a separate question and is not covered by the MIT grant.
This deployment contains curated kinetic values with cited PubMed/DOI sources and
does not contain any IMPPAT or NPASS-derived data. IMPPAT is CC BY-NC-ND 4.0 and its
derived data must not be published without permission. Do not add the ingest output
(`imppat_index.json`) to a public deployment.

Parameters marked `illustrative` are placeholders and cap every prediction that uses
them at LOW confidence.
