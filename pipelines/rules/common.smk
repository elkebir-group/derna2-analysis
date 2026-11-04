# Provenance: binaries and deposits checked against pins/, and derna2 built at a pin.


# The deposited results of one experiment, against pins/results-{exp}.sha256.
rule check_results:
    input:
        manifest=ROOT / "pins" / "results-{exp}.sha256",
    output:
        STAMPS / "results_{exp}.ok",
    params:
        root=lambda w: RESULTS / w.exp,
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {params.root} {output}"


# The two derna2 binaries every run target uses, built from source at `derna2.pin`.
# A host that has them already names their directory as `derna2.bin_dir` instead.
rule build_derna2:
    output:
        derna=DERIVED / "bin" / "derna2-{pin}" / "derna",
        fold=DERIVED / "bin" / "derna2-{pin}" / "mfe_fold",
        sums=DERIVED / "bin" / "derna2-{pin}" / "SHA256SUMS",
        commit=DERIVED / "bin" / "derna2-{pin}" / "COMMIT",
    params:
        repo=config["derna2"]["repo"],
        out=lambda w: DERIVED / "bin" / f"derna2-{w.pin}",
    log:
        DERIVED / "logs" / "build_derna2_{pin}.log",
    shell:
        "bash {SCRIPTS}/build_derna2.sh {params.repo} {wildcards.pin} {params.out} > {log} 2>&1"


# The deposit of one experiment as a tar of exactly its manifest's files. Unpack it under
# data/results/ (the archive's paths start with the experiment's name).
rule package_deposit:
    input:
        manifest=ROOT / "pins" / "results-{exp}.sha256",
    output:
        DERIVED / "deposit" / "{exp}.tar",
    params:
        root=lambda w: RESULTS / w.exp,
    shell:
        "{PY} {SCRIPTS}/package_deposit.py {input.manifest} {params.root} {output} {wildcards.exp}"
