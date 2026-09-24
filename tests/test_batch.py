"""Batch mode: many families in one run, each input read once."""

from __future__ import annotations

import csv

import pytest

from teaid import cli

FASTA = (">fam1#LTR/Gypsy\n" + "ACGTTGCA" * 40 + "\n>fam2#Unknown\n" + "TTTTAAAA" * 20
         + "\n>fam3#DNA/hAT\n" + "GGCCAATT" * 30 + "\n")

RMOUT = """\
   SW   perc perc perc  query       position in query              matching          repeat                position in repeat
score   div. del. ins.  sequence    begin    end          (left)   repeat            class/family      begin   end    (left)     ID

 1000    5.0  0.0  0.0  chr1            1000     1319    (5000) + fam1              LTR/Gypsy               1    320     (0)     1
  800   12.5  0.0  0.0  chr1            2000     2159    (4000) C fam1              LTR/Gypsy             (160)   160       1     2
  700    8.0  0.0  0.0  chr2             100      339    (4000) + fam3              DNA/hAT                 1    240     (0)     3
"""

QUIET = ["--no-dotplot", "--no-orfs", "--no-homology"]


@pytest.fixture
def inputs(tmp_path):
    fasta = tmp_path / "families.fa"
    fasta.write_text(FASTA)
    annot = tmp_path / "genome.fa.out"
    annot.write_text(RMOUT)
    return fasta, annot


def manifest(path):
    with (path / "results.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TestBatch:
    def test_all_draws_every_family_with_copies_and_says_what_it_skipped(
            self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        out = tmp_path / "sheets"
        code = cli.main(["--annot", str(annot), "-c", str(fasta), "--all", "-o", str(out), *QUIET])
        assert code == cli.EXIT_OK
        rows = manifest(out)
        assert [r["family"] for r in rows] == ["fam1", "fam3"]
        assert all(r["exit"] == "0" and r["slug"] == "" for r in rows)
        assert {r["copies"] for r in rows} == {"2", "1"}
        assert (out / "fam1.teaid.html").exists() and (out / "fam3.teaid.html").exists()
        assert (out / "index.html").exists()
        assert "1 library families have no annotated copies" in capsys.readouterr().err

    def test_a_failed_family_is_recorded_and_the_batch_exits_8(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        wanted = tmp_path / "wanted.txt"
        wanted.write_text("fam1\tany\textra columns\n# a comment\nfam2\nabsent\n")
        out = tmp_path / "sheets"
        code = cli.main(["--annot", str(annot), "-c", str(fasta), "--families", str(wanted),
                         "-o", str(out), *QUIET])
        assert code == cli.EXIT_BATCH_INCOMPLETE
        assert "[batch-incomplete]" in capsys.readouterr().err
        rows = {r["family"]: r for r in manifest(out)}
        assert rows["fam1"]["exit"] == "0"
        assert rows["fam2"]["slug"] == "no-evidence"
        assert rows["absent"]["slug"] == "no-family"

    def test_stdout_keeps_one_line_per_file_and_stderr_one_per_family(
            self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        out = tmp_path / "sheets"
        cli.main(["--annot", str(annot), "-c", str(fasta), "--all", "-o", str(out), *QUIET])
        captured = capsys.readouterr()
        wrote = [ln for ln in captured.out.splitlines() if ln.startswith("wrote ")]
        assert sorted(ln.rsplit("/", 1)[1] for ln in wrote) == [
            "fam1.teaid.html", "fam3.teaid.html", "index.html", "results.tsv"]
        progress = [ln for ln in captured.err.splitlines() if "/2] " in ln]
        assert progress[0].startswith("teaid: [1/2] fam1: ok (")

    def test_family_and_all_do_not_combine(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        code = cli.main(["--annot", str(annot), "-c", str(fasta), "--all", "-f", "fam1",
                         "-o", str(tmp_path)])
        assert code == cli.EXIT_USAGE
        assert "[usage]" in capsys.readouterr().err

    def test_parallel_workers_give_the_same_manifest(self, inputs, tmp_path):
        fasta, annot = inputs
        serial, parallel = tmp_path / "one", tmp_path / "two"
        cli.main(["--annot", str(annot), "-c", str(fasta), "--all", "-o", str(serial), *QUIET])
        cli.main(["--annot", str(annot), "-c", str(fasta), "--all", "-j", "2",
                  "-o", str(parallel), *QUIET])

        def comparable(rows):
            return [{k: v for k, v in r.items() if k != "seconds"} for r in rows]

        assert comparable(manifest(serial)) == comparable(manifest(parallel))

    def test_no_index_skips_the_contact_sheet(self, inputs, tmp_path):
        fasta, annot = inputs
        cli.main(["--annot", str(annot), "-c", str(fasta), "--all", "--no-index",
                  "-o", str(tmp_path), *QUIET])
        assert (tmp_path / "results.tsv").exists() and not (tmp_path / "index.html").exists()


class TestStems:
    def test_names_that_would_share_a_file_get_suffixes(self):
        """One family per run could leave this to the caller; a batch cannot."""
        stems = cli._unique_stems(["fam/a", "fam:a", "fam_a", "other"])
        assert stems == {"fam/a": "fam_a", "fam:a": "fam_a.2", "fam_a": "fam_a.3",
                         "other": "other"}


SEEDS = "".join(
    "# STOCKHOLM 1.0\n#=GF ID    {name}\n#=GC RF    {rf}\n"
    "chr1:1-{n}_+    {rf}\nchr1:1001-{end}_-    {rf}\n//\n".format(
        name=name, rf=rf, n=len(rf), end=1000 + len(rf))
    for name, rf in (("seedA", "ACGTTGCA" * 30), ("seedB", "GGCCAATT" * 25))
)


class TestSeedBatch:
    def test_every_record_of_a_seed_file_gets_a_sheet(self, tmp_path):
        stk = tmp_path / "families.stk"
        stk.write_text(SEEDS)
        out = tmp_path / "sheets"
        assert cli.main(["--stk", str(stk), "--all", "-o", str(out), *QUIET]) == cli.EXIT_OK
        assert [r["family"] for r in manifest(out)] == ["seedA", "seedB"]
        assert (out / "seedA.teaid.html").exists() and (out / "seedB.teaid.html").exists()
