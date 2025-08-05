#!/usr/bin/env Rscript

# generate_azimuth_swarm.R
# Generate swarm commands for batch Azimuth annotation

library(stringr)

# Define classes and reference sets
classes <- c("other", "inn", "exn_upper", "exn_lower")
refs <- c(
  "allen_M1", "allen_MCA",
  "BRAIN-other", "BRAIN-exn", "BRAIN-inn"
)

# Create or overwrite swarm file
swarm_file <- "azimuth.swarm"
if (file.exists(swarm_file)) {
  file.remove(swarm_file)
}

# Write swarm lines
for (ref in refs) {
  for (class in classes) {
    line <- paste("Rscript azimuth.R", class, ref)
    write(line, file = swarm_file, append = TRUE)
  }
}

cat("Swarm file 'azimuth.swarm' generated.\n")
