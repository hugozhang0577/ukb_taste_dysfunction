#!/bin/bash
# Stage the PRSice-2 release into the project first; RAP workers have no
# internet, so the binary is fetched from project storage rather than a URL.
dx download "/gwas/cohort_primary/PRS/PRSice_linux.zip" -o PRSice_linux.zip && unzip -o PRSice_linux.zip && chmod +x PRSice_linux && ./PRSice_linux --version
