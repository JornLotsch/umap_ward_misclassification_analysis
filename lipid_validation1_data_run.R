xx <- read.csv("/home/joern/Schreibtisch/ST000048.txt", skip = 186, sep = "\t")

# Define the path to your data file
file_path <- "/home/joern/Schreibtisch/ST000048.txt"



# 1. Read the header line (line 185) to get column names
headers <- read.table(file_path, 
                      header = FALSE, 
                      sep = "\t", 
                      skip = 3,     # Skip lines before header
                      nrows = 1,      # Read only the header line
                      stringsAsFactors = FALSE)

# 2. Read the actual data (lines 187–358)
metabaolmics_study_ST004430_data <- read.table(file_path, 
                   header = FALSE, 
                   sep = "\t", 
                   skip = 186,        # Skip up to line before data
                   nrows = 358 - 186, # Read lines 187–358
                   stringsAsFactors = FALSE,
                   fill = TRUE)       # Fills missing entries

# 3. Apply headers to the data
colnames(metabaolmics_study_ST004430_data) <- unlist(headers[1, ])


# Inspect your imported dataset
str(metabaolmics_study_ST004430_data)
head(metabaolmics_study_ST004430_data)

library(tidyverse)

# 1. You already have your metabolomics data loaded
# metabaolmics_study_ST004430_data <- read.table(...)

# 2. Transpose the file so that each row = sample, each column = metabolite
metabolomics_wide <- metabaolmics_study_ST004430_data %>%
  pivot_longer(
    cols = -Samples,
    names_to = "local_sample_id",
    values_to = "intensity"
  ) %>%
  pivot_wider(
    names_from = Samples,
    values_from = intensity
  )

# Make sure IDs are numeric (to match metadata)
metabolomics_wide$local_sample_id <- as.numeric(metabolomics_wide$local_sample_id)

# 3. Read in metadata CSV
metadata <- read.csv("/home/joern/Schreibtisch/metabaolmics_study_ST004430_metadata.csv", stringsAsFactors = FALSE)

# 4. Align to ensure matching order by local_sample_id
metadata <- metadata[order(metadata$local_sample_id), ]
metabolomics_wide <- metabolomics_wide[order(metabolomics_wide$local_sample_id), ]

# 5. Merge both objects (if desired)
combined_data <- left_join(metadata, metabolomics_wide, by = "local_sample_id")

# 6. Check that samples match
stopifnot(all(metadata$local_sample_id == metabolomics_wide$local_sample_id))

# Optional: quick peek
head(combined_data)

apply(combined_data,1,function(x) sum(is.na(x)))
