# Requirements

* Download the E. coli K-12 MG1655 version 3 genbank reference genome (https://www.ncbi.nlm.nih.gov/nuccore/556503834) and place it in ./data with the name "NC_000913_3.gb".  
* Run the script aledb_parse_annot_assoc.sh to generate necessary intermediate files used by other analysis.  
* Remaining notebooks process data and generate results for ALE mutation analysis.

The most recent published figures used from repo used the following software versions:  
* Ubuntu 20.04
* R 3.6.3
    - ComplexHeatmap 2.9.1
    - trackViewer 1.22.1
* Python 3.7.7
    - pandas 1.1.3
    - numpy 1.19.2
    - matplotlib 3.3.2
    - scipy 1.5.2
    - seaborn 0.11.0
* javascript
    - NGL 2.0.0-dev
