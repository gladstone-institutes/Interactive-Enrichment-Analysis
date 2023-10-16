# Fetches GMT files for building database collections (GH Actions)

## Datestamp for downloads
DATE=$(date --utc +%Y%m%d) #change argto --utc (linux), -u (macOS)

## WikiPathways
WP_REL=$(curl --silent https://data.wikipathways.org/current/gmt/ | grep -o -m 1 '>wikipathways-.*-gmt' | sed 's/>wikipathways-//;s/-gmt//' | head -n 1)  
wget -O "shiny_run/databases/gmts/wp_hs_"${DATE}".gmt" "https://data.wikipathways.org/current/gmt/wikipathways-"${WP_REL}"-gmt-Homo_sapiens.gmt"
wget -O "shiny_run/databases/gmts/wp_mm_"${DATE}".gmt" "https://data.wikipathways.org/current/gmt/wikipathways-"${WP_REL}"-gmt-Mus_musculus.gmt"

##PFOCR
PFOCR_DATE=$(curl --silent https://data.wikipathways.org/pfocr/ | grep -oP '<a class=file-link href="\K\d+' | tail -n 1)
wget -O "shiny_run/databases/gmts/pfocr_hs_"${DATE}".gmt" "https://data.wikipathways.org/pfocr/"${PFOCR_DATE}"/pfocr-"${PFOCR_DATE}"-gmt-Homo_sapiens.gmt"
wget -O "shiny_run/databases/gmts/pfocr_mm_"${DATE}".gmt" "https://data.wikipathways.org/pfocr/"${PFOCR_DATE}"/pfocr-"${PFOCR_DATE}"-gmt-Mus_musculus.gmt"

#Gene Ontology
wget -O "shiny_run/databases/gmts/go_hs_"${DATE}".gmt" "http://download.baderlab.org/EM_Genesets/current_release/Human/entrezgene/GO/Human_GO_bp_no_GO_iea_entrezgene.gmt"
wget -O "shiny_run/databases/gmts/go_mm_"${DATE}".gmt" "http://download.baderlab.org/EM_Genesets/current_release/Mouse/entrezgene/GO/MOUSE_GO_bp_no_GO_iea_entrezgene.gmt"

