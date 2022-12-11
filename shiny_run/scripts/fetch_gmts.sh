# Fetches GMT files for building database collections (GH Actions)

## Datestamp for downloads
DATE=$(date --utc +%Y%m%d) #change argto --utc (linux), -u (macOS)

## WikiPathways
WP_REL=$(curl --silent https://wikipathways-data.wmcloud.org/current/gmt/ | grep -o -m 1 'wikipathways-\d*' | sed 's/wikipathways-//' | head -n 1)  
wget -O "shiny_run/databases/gmts/wp_hs_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-"${WP_REL}"-gmt-Homo_sapiens.gmt"
wget -O "shiny_run/databases/gmts/wp_mm_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-"${WP_REL}"-gmt-Mus_musculus.gmt"

##PFOCR
PFOCR_REL=$(curl --silent https://wikipathways-data.wmcloud.org/pfocr/ | grep -o -m 1 '>pfocr-.*-gmt' | sed 's/>pfocr-//;s/-gmt//' | tail -n 1)
wget -O "shiny_run/databases/gmts/pfocr_hs_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/pfocr/pfocr-"${PFOCR_REL}"-gmt-Homo_sapiens.gmt"
wget -O "shiny_run/databases/gmts/pfocr_mm_"${DATE}".gmt" "https://wikipathways-data.wmcloud.org/pfocr/pfocr-"${PFOCR_REL}"-gmt-Mus_musculus.gmt"

#Gene Ontology
GO_HS_REL=$(curl --silent https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp | grep -o 'c5\.go\.bp\.v.*\.entrez\.gmt' | sed 's/c5\.go\.bp\.v//;s/\.entrez\.gmt//')
GO_MM_REL=$(curl --silent https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp | grep -o 'm5\.go\.bp\.v.*\.entrez\.gmt' | sed 's/m5\.go\.bp\.v//;s/\.entrez\.gmt//')
wget -O "shiny_run/databases/gmts/go_hs_"${DATE}".gmt" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"${GO_HS_REL}"/c5.go.bp.v"${GO_HS_REL}".entrez.gmt"
wget -O "shiny_run/databases/gmts/go_mm_"${DATE}".gmt" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"${GO_MM_REL}"/m5.go.bp.v"${GO_MM_REL}".entrez.gmt"


