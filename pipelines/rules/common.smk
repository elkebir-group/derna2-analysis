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


# Frozen binaries of one experiment, against pins/{exp}.sha256.
rule check_bins:
    input:
        manifest=ROOT / "pins" / "{exp}.sha256",
    output:
        STAMPS / "bins_{exp}.ok",
    params:
        root=lambda w: Path(config[w.exp]["bin_dir"]),
    wildcard_constraints:
        exp="e1|e2|fig1",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {params.root} {output}"


# derna2's benchmark binaries rebuilt from source at a pin. Whether a rebuild reproduces
# the frozen binaries byte for byte is what check_derna2_build reports.
rule build_derna2:
    output:
        sums=DERIVED / "bin" / "derna2-{pin}" / "SHA256SUMS",
    params:
        repo=config["derna2"]["repo"],
        out=lambda w: DERIVED / "bin" / f"derna2-{w.pin}",
    log:
        DERIVED / "logs" / "build_derna2_{pin}.log",
    shell:
        "bash {SCRIPTS}/build_derna2.sh {params.repo} {wildcards.pin} {params.out} > {log} 2>&1"


rule check_derna2_build:
    input:
        sums=lambda w: DERIVED / "bin" / f"derna2-{config[w.exp]['pin']}" / "SHA256SUMS",
        manifest=ROOT / "pins" / "{exp}.sha256",
    output:
        DERIVED / "bin" / "derna2-{exp}.compare.txt",
    wildcard_constraints:
        exp="e1|e2",
    run:
        pinned = dict(reversed(line.split()) for line in open(input.manifest) if line.strip())
        lines = []
        for line in open(input.sums):
            digest, name = line.split()
            want = pinned.get(name)
            verdict = "absent from pin" if want is None else (
                "identical" if want == digest else "DIFFERS")
            lines.append(f"{name}\t{verdict}")
        Path(output[0]).write_text("\n".join(lines) + "\n")
        print("\n".join(lines))


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
