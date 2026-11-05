#!/bin/bash
# ISOLDE's ParmEd conversions of ff14SB and of the Bryce-group set (Meagher 2003 GTP) are not redistributed here
# (custom ISOLDE licence).  Fetch them at the commit used for every calculation, then build gtp_meagher.xml.
set -e
cd "$(dirname "$0")"
SHA=2b9372ad2cc5c05292e6d02b928a5f76a9ad482b
for f in amberff14SB.xml bryce_set.xml; do
  curl -sSfL -o $f https://raw.githubusercontent.com/tristanic/isolde/$SHA/isolde/src/openmm/amberff/$f
done
curl -sSfL -o openmm_amber14_tip3p.xml https://raw.githubusercontent.com/openmm/openmm/8.1.1/wrappers/python/openmm/app/data/amber14/tip3p.xml
python3 ../scripts/a00_make_ff.py
md5sum amberff14SB.xml bryce_set.xml gtp_meagher.xml 2>/dev/null || md5 amberff14SB.xml bryce_set.xml gtp_meagher.xml
