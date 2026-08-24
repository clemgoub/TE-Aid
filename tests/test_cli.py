"""End-to-end CLI tests, including the exit codes the pipeline contract relies on.

TE-Aid is invoked per family by a seed-building pipeline that must keep going
when one family fails, so a failure has to be distinguishable from a success by
exit code alone.
"""

from __future__ import annotations

import pytest

from teaid import cli

FASTA = ">fam1#LTR/Gypsy\n" + "ACGTTGCA" * 40 + "\n>fam2#Unknown\n" + "TTTTAAAA" * 20 + "\n"

RMOUT = """\
   SW   perc perc perc  query       position in query              matching          repeat                position in repeat
score   div. del. ins.  sequence    begin    end          (left)   repeat            class/family      begin   end    (left)     ID

 1000    5.0  0.0  0.0  chr1            1000     1319    (5000) + fam1              LTR/Gypsy               1    320     (0)     1
  800   12.5  0.0  0.0  chr1            2000     2159    (4000) C fam1              LTR/Gypsy             (160)   160       1     2
"""


@pytest.fixture
def inputs(tmp_path):
    fasta = tmp_path / "families.fa"
    fasta.write_text(FASTA)
    annot = tmp_path / "genome.fa.out"
    annot.write_text(RMOUT)
    return fasta, annot


def run(args) -> int:
    return cli.main([str(a) for a in args])


class TestSuccessPath:
    def test_writes_a_sheet_and_exits_zero(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        out = tmp_path / "sheets"
        code = run(["--annot", annot, "-c", fasta, "-f", "fam1", "-o", out, "--no-dotplot"])
        assert code == cli.EXIT_OK
        assert (out / "fam1.teaid.html").exists()
        assert "2 copies" in capsys.readouterr().out

    def test_html_is_self_describing(self, inputs, tmp_path):
        fasta, annot = inputs
        out = tmp_path / "sheets"
        run(["--annot", annot, "-c", fasta, "-f", "fam1", "-o", out, "--no-dotplot"])
        html = (out / "fam1.teaid.html").read_text()
        assert "fam1" in html
        assert "divergence from consensus" in html
        # The sheet renders evidence and never asserts a classification.
        assert "predicted class" not in html.casefold()

    def test_family_resolves_across_a_class_suffix(self, inputs, tmp_path):
        fasta, annot = inputs
        code = run(
            ["--annot", annot, "-c", fasta, "-f", "fam1#LTR/Gypsy", "-o", tmp_path, "--no-dotplot"]
        )
        assert code == cli.EXIT_OK


class TestFailurePaths:
    def test_missing_annotation_exits_no_input(self, inputs, tmp_path):
        fasta, _ = inputs
        code = run(["--annot", tmp_path / "absent.out", "-c", fasta, "-f", "fam1", "-o", tmp_path])
        assert code == cli.EXIT_NO_INPUT

    def test_unknown_family_exits_no_family_and_suggests(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        code = run(["--annot", annot, "-c", fasta, "-f", "fam11", "-o", tmp_path])
        assert code == cli.EXIT_NO_FAMILY
        assert "did you mean" in capsys.readouterr().err

    def test_family_with_no_annotated_copies_exits_no_evidence(self, inputs, tmp_path, capsys):
        """fam2 has a consensus but no rows in the annotation."""
        fasta, annot = inputs
        code = run(["--annot", annot, "-c", fasta, "-f", "fam2", "-o", tmp_path])
        assert code == cli.EXIT_NO_EVIDENCE
        assert "no copies" in capsys.readouterr().err

    def test_family_is_required_for_a_multi_entry_library(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        code = run(["--annot", annot, "-c", fasta, "-o", tmp_path])
        assert code == cli.EXIT_USAGE
        assert "--family is required" in capsys.readouterr().err

    def test_family_is_optional_for_a_single_entry_library(self, tmp_path):
        fasta = tmp_path / "one.fa"
        fasta.write_text(">fam1#LTR/Gypsy\n" + "ACGTTGCA" * 40 + "\n")
        annot = tmp_path / "genome.fa.out"
        annot.write_text(RMOUT)
        assert run(["--annot", annot, "-c", fasta, "-o", tmp_path, "--no-dotplot"]) == cli.EXIT_OK

    def test_out_of_range_threshold_is_rejected(self, inputs, tmp_path):
        fasta, annot = inputs
        with pytest.raises(SystemExit):
            run(["--annot", annot, "-c", fasta, "-f", "fam1", "-t", "1.5", "-o", tmp_path])


class TestLegacyAlias:
    def test_alias_warns_and_still_works(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        code = cli.main_legacy_alias(
            ["--annot", str(annot), "-c", str(fasta), "-f", "fam1", "-o", str(tmp_path), "--no-dotplot"]
        )
        assert code == cli.EXIT_OK
        assert "'TE-Aid' is now 'teaid'" in capsys.readouterr().err
