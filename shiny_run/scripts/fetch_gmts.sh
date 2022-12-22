# Fetches GMT files for building database collections (GH Actions)

## Datestamp for downloads
DATE=$(date --utc +%Y%m%d) #change argto --utc (linux), -u (macOS)

## WikiPathways
WP_REL=$(curl --silent https://wikipathways-data.wmcloud.org/current/gmt/ | grep -o -m 1 '>wikipathways-.*-gmt' | sed 's/>wikipathways-//;s/-gmt//' | head -n 1)  
wget -O "shiny_run/databases/gmts/wp_hs_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-"${WP_REL}"-gmt-Homo_sapiens.gmt"
wget -O "shiny_run/databases/gmts/wp_mm_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-"${WP_REL}"-gmt-Mus_musculus.gmt"

##PFOCR
PFOCR_REL=$(curl --silent https://wikipathways-data.wmcloud.org/pfocr/ | grep -o -m 1 '>pfocr-.*-gmt' | sed 's/>pfocr-//;s/-gmt//' | tail -n 1)
wget -O "shiny_run/databases/gmts/pfocr_hs_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/pfocr/pfocr-"${PFOCR_REL}"-gmt-Homo_sapiens.gmt"
wget -O "shiny_run/databases/gmts/pfocr_mm_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/pfocr/pfocr-"${PFOCR_REL}"-gmt-Mus_musculus.gmt"

#Gene Ontology
wget -O "shiny_run/databases/gmts/go_hs_"${DATE}".gmt" "http://download.baderlab.org/EM_Genesets/current_release/Human/entrezgene/GO/Human_GO_bp_no_GO_iea_entrezgene.gmt"
wget -O "shiny_run/databases/gmts/go_mm_"${DATE}".gmt" "http://download.baderlab.org/EM_Genesets/current_release/Mouse/entrezgene/GO/MOUSE_GO_bp_no_GO_iea_entrezgene.gmt"

