'''
Run rifdock, either locally or on SGE.

Usage:
    run_rifdock.py <folder> [Options]

Options:
    --job-distributor=STR, -d  Running on SGE or locally? Defaultse to
    local.
    --tasks=NUM, -n  How many jobs to split into if running on SGE.
'''
