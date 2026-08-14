#!/bin/bash
# =============================================================================
# PRS — put PRSice-2 in the JupyterLab workspace
# =============================================================================
#
# PRSice-2 is distributed as a zip holding a precompiled Linux binary, so there
# is nothing to build and nothing to install: upload the release zip to the
# project once, then unzip it in the workspace.
#
# Two things to get right, both of which fail quietly otherwise:
#   - unzip does not reliably preserve the executable bit, so chmod is not
#     optional even though the file came out of an official release;
#   - the binary is dynamically linked, so it must run on the worker, not be
#     copied to some other machine.
#
# The version is pinned by which zip was uploaded, which is the point: a URL
# pointing at "latest" would silently change the software between runs.
#
# Run these lines in the JupyterLab terminal.
#
# Input:  PRSice_linux.zip, uploaded to $PROJECT_PRS in the project
# Output: ./PRSice_linux (executable) and ./PRSice.R in the workspace
# =============================================================================

# Bring the release into the workspace and unpack it

dx download "/gwas/cohort_primary/PRS/PRSice_linux.zip" -o PRSice_linux.zip && unzip -o PRSice_linux.zip && chmod +x PRSice_linux && ./PRSice_linux --version
