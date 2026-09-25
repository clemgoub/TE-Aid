"""End-to-end CLI tests, including the exit codes the pipeline contract relies on.

TE-Aid is invoked per family by a seed-building pipeline that must keep going
when one family fails, so a failure has to be distinguishable from a success by
exit code alone.
"""

from __future__ import annotations

import re
from pathlib import Path

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
        code = run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", out, "--no-dotplot"])
        assert code == cli.EXIT_OK
        assert (out / "fam1.teaid.html").exists()
        assert "2 copies" in capsys.readouterr().out

    def test_html_is_self_describing(self, inputs, tmp_path):
        fasta, annot = inputs
        out = tmp_path / "sheets"
        run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", out, "--no-dotplot"])
        html = (out / "fam1.teaid.html").read_text()
        assert "fam1" in html
        assert "divergence from consensus" in html
        # The sheet renders evidence and never asserts a classification.
        assert "predicted class" not in html.casefold()

    def test_family_resolves_across_a_class_suffix(self, inputs, tmp_path):
        fasta, annot = inputs
        code = run(
            ["--annot", annot, "-c", fasta, "-n", "fam1#LTR/Gypsy", "-o", tmp_path, "--no-dotplot"]
        )
        assert code == cli.EXIT_OK


class TestFailurePaths:
    def test_missing_annotation_exits_no_input(self, inputs, tmp_path):
        fasta, _ = inputs
        code = run(["--annot", tmp_path / "absent.out", "-c", fasta, "-n", "fam1", "-o", tmp_path])
        assert code == cli.EXIT_NO_INPUT

    def test_unknown_family_exits_no_family_and_suggests(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        code = run(["--annot", annot, "-c", fasta, "-n", "fam11", "-o", tmp_path])
        assert code == cli.EXIT_NO_FAMILY
        assert "did you mean" in capsys.readouterr().err

    def test_family_with_no_annotated_copies_exits_no_evidence(self, inputs, tmp_path, capsys):
        """fam2 has a consensus but no rows in the annotation."""
        fasta, annot = inputs
        code = run(["--annot", annot, "-c", fasta, "-n", "fam2", "-o", tmp_path])
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
            run(["--annot", annot, "-c", fasta, "-n", "fam1", "-f", "1.5", "-o", tmp_path])


SEED = """\
# STOCKHOLM 1.0
#=GF ID    fam1
#=GF TP    Interspersed_Repeat;Transposable_Element
#=GF SQ    3
#=GC RF    {rf}
chr1:1-{n}_+    {rf}
chr1:1001-{end}_-    {rf}
chr1:2001-{end2}_+    {rf}
//
""".format(rf="ACGTTGCA" * 40, n=320, end=1320, end2=2320)


@pytest.fixture
def seed_file(tmp_path):
    path = tmp_path / "families.stk"
    path.write_text(SEED)
    return path


class TestSeedInput:
    def test_seed_alone_produces_a_sheet(self, seed_file, tmp_path, capsys):
        """The seed carries copies, loci, consensus and alignment, so no
        separate FASTA is needed."""
        code = run(["--stk", seed_file, "-n", "fam1", "-o", tmp_path, "--no-dotplot", "--no-orfs"])
        assert code == cli.EXIT_OK
        assert (tmp_path / "fam1.teaid.html").exists()
        assert "3 copies" in capsys.readouterr().out

    def test_seed_qc_is_opt_in(self, seed_file, tmp_path):
        base = ["--stk", seed_file, "-n", "fam1", "-o", tmp_path, "--no-dotplot", "--no-orfs"]
        run(base)
        without = (tmp_path / "fam1.teaid.html").read_text()
        run(base + ["--seed-qc"])
        with_qc = (tmp_path / "fam1.teaid.html").read_text()
        assert "Seed depth" not in without
        assert "Seed depth" in with_qc
        assert "#=GF TP" in with_qc

    def test_pipeline_mode_prefers_the_seed(self, seed_file, inputs, tmp_path):
        """--pipeline reverses the standalone priority to seed-first, so the
        seed supplies the consensus and the seed-QC panels."""
        fasta, annot = inputs
        code = run([
            "--pipeline", "--stk", seed_file, "--annot", annot, "-c", fasta,
            "-n", "fam1", "-o", tmp_path, "--seed-qc", "--no-dotplot", "--no-orfs",
        ])
        assert code == cli.EXIT_OK
        assert "Seed depth" in (tmp_path / "fam1.teaid.html").read_text()

    def test_seed_and_annotation_together_are_not_alternatives(
        self, seed_file, inputs, tmp_path, capsys
    ):
        """The seed says which copies were *chosen*; the annotation says which
        *exist*. Given both, panels 1-2 show the genome and panel 6 shows the
        seed's sampling — the comparison is the whole point of seed QC."""
        fasta, annot = inputs
        code = run([
            "--stk", seed_file, "--annot", annot, "-n", "fam1",
            "-o", tmp_path, "--seed-qc", "--no-dotplot", "--no-orfs",
        ])
        assert code == cli.EXIT_OK
        # The .out fixture holds 2 rows for fam1; the seed holds 3 sequences.
        out = capsys.readouterr().out
        assert "2 copies" in out
        html = (tmp_path / "fam1.teaid.html").read_text()
        assert "2 annotated genomic copies" in html
        assert "the seed uses 3" in html

    def test_a_family_absent_from_the_annotation_still_renders(
        self, seed_file, inputs, tmp_path, capsys
    ):
        """Worth warning about, not worth refusing the sheet for."""
        fasta, annot = inputs
        seed = tmp_path / "other.stk"
        seed.write_text(SEED.replace("fam1", "orphan"))
        code = run([
            "--stk", seed, "--annot", annot, "-n", "orphan",
            "-o", tmp_path, "--no-dotplot", "--no-orfs",
        ])
        assert code == cli.EXIT_OK
        assert "fall back to the seed" in capsys.readouterr().err

    def test_consensus_is_ignored_with_a_seed(self, seed_file, inputs, tmp_path, capsys):
        """A FASTA passed alongside would silently disagree with the alignment
        the depth panel is drawn from."""
        fasta, _ = inputs
        code = run([
            "--stk", seed_file, "-c", fasta, "-n", "fam1",
            "-o", tmp_path, "--no-dotplot", "--no-orfs",
        ])
        assert code == cli.EXIT_OK
        assert "--consensus ignored" in capsys.readouterr().err

    def test_unknown_record_exits_bad_seed_and_lists_what_is_there(
        self, seed_file, tmp_path, capsys
    ):
        code = run(["--stk", seed_file, "-n", "absent", "-o", tmp_path])
        assert code == cli.EXIT_BAD_SEED
        err = capsys.readouterr().err
        assert "[bad-seed]" in err
        assert "fam1" in err

    def test_malformed_seed_exits_bad_seed(self, tmp_path, capsys):
        path = tmp_path / "broken.stk"
        path.write_text("this is not stockholm\n")
        assert run(["--stk", path, "-o", tmp_path]) == cli.EXIT_BAD_SEED
        assert "[bad-seed]" in capsys.readouterr().err

    def test_missing_seed_file_exits_no_input(self, tmp_path, capsys):
        assert run(["--stk", tmp_path / "nope.stk", "-o", tmp_path]) == cli.EXIT_NO_INPUT
        assert "[no-input]" in capsys.readouterr().err


class TestFailSoftContract:
    """The pipeline queues a failed packet for a curator, so every failure has
    to be machine-readable: a distinct exit code and a stable stderr slug."""

    def test_each_failure_carries_a_distinct_code_and_slug(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        cases = [
            (["--annot", tmp_path / "absent.out", "-c", fasta, "-n", "fam1"],
             cli.EXIT_NO_INPUT, "no-input"),
            (["--annot", annot, "-c", fasta, "-n", "nope"],
             cli.EXIT_NO_FAMILY, "no-family"),
            (["--annot", annot, "-c", fasta, "-n", "fam2"],
             cli.EXIT_NO_EVIDENCE, "no-evidence"),
            (["--annot", annot, "-c", fasta],
             cli.EXIT_USAGE, "usage"),
        ]
        seen = set()
        for args, code, slug in cases:
            assert run(args + ["-o", tmp_path]) == code
            assert f"[{slug}]" in capsys.readouterr().err
            seen.add(code)
        assert len(seen) == len(cases), "exit codes must be distinguishable"

    def test_an_argparse_error_also_carries_the_usage_slug(self, capsys):
        with pytest.raises(SystemExit) as raised:
            run(["--no-such-flag"])
        assert raised.value.code == cli.EXIT_USAGE
        assert "[usage]" in capsys.readouterr().err

    def test_the_documented_exit_table_matches_the_code(self):
        # INTEGRATION.md is what callers read, so its table is the contract.
        doc = Path(__file__).resolve().parents[1] / "docs" / "INTEGRATION.md"
        rows = re.findall(r"^\| (\d) \| ([^|]+?) \|", doc.read_text(), re.M)
        documented = {int(code): cell.strip().strip("`") for code, cell in rows}
        for code, slug in cli._SLUG.items():
            assert documented.get(code) == slug, f"exit {code}"
        assert documented[1] == "*none*"
        assert set(documented) == set(cli._SLUG) | {0, 1}

    def test_slugs_and_codes_stay_in_step(self):
        for code in (
            cli.EXIT_USAGE, cli.EXIT_NO_INPUT, cli.EXIT_NO_FAMILY,
            cli.EXIT_NO_EVIDENCE, cli.EXIT_BAD_SEED,
        ):
            assert code in cli._SLUG
        assert len(set(cli._SLUG.values())) == len(cli._SLUG)


class TestLegacyAlias:
    def test_alias_warns_and_still_works(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        code = cli.main_legacy_alias(
            ["--annot", str(annot), "-c", str(fasta), "-n", "fam1", "-o", str(tmp_path), "--no-dotplot"]
        )
        assert code == cli.EXIT_OK
        err = capsys.readouterr().err
        assert "'TE-Aid' is now 'teaid'" in err and "Every v1 flag is accepted" in err


class TestProteinProvenance:
    """Panel 5 must say which library produced it, and why tier 3 did not run."""

    @pytest.fixture
    def captured(self, monkeypatch):
        from teaid import homology, orfs, report

        seen = {"repeatpeps": None}

        def fake_repeatpeps(found, path, **kw):
            seen["repeatpeps"] = path
            return []

        monkeypatch.setattr(homology, "search", lambda *a, **k: [])
        monkeypatch.setattr(homology, "search_repeatpeps", fake_repeatpeps)
        monkeypatch.setattr(orfs, "find_orfs", lambda *a, **k: [])
        original = report.build_figure

        def keep(data, theme):
            seen["data"] = data
            return original(data, theme)

        monkeypatch.setattr(report, "build_figure", keep)
        return seen

    def test_a_prebuilt_library_still_honours_repeatpeps(self, inputs, tmp_path, captured):
        fasta, annot = inputs
        prebuilt = tmp_path / "mine.bhmm"
        prebuilt.write_text("")
        peps = tmp_path / "RepeatPeps.lib"
        peps.write_text(">x\nMK\n")
        run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "--no-dotplot",
             "--proteins", prebuilt, "--repeatpeps", peps])
        assert captured["repeatpeps"] == peps

    def test_the_sheet_names_the_protein_library_it_searched(self, inputs, tmp_path, captured):
        fasta, annot = inputs
        prebuilt = tmp_path / "mine.bhmm"
        prebuilt.write_text("")
        run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "--no-dotplot",
             "--proteins", prebuilt])
        assert any("mine.bhmm" in s for s in captured["data"].sources)

    def test_a_tier_three_that_did_not_run_says_why(self, inputs, tmp_path, captured, monkeypatch):
        from teaid import proteins

        fasta, annot = inputs
        library = proteins.Library(hmm=tmp_path / "x.bhmm", repeatpeps=None,
                                   tiers=("pfam:130", "repeatpeps-gaps"))
        monkeypatch.setattr(proteins, "build", lambda **k: library)
        base = ["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "--no-dotplot"]
        for extra, reason in (([], "RepeatPeps.lib not found"),
                              (["--no-orfs"], "--no-orfs"),
                              (["--deep"], "--deep")):
            run(base + extra)
            notes = " | ".join(captured["data"].notes)
            assert "tier 3" in notes and reason in notes, (extra, notes)


SPLIT = """\
   SW   perc perc perc  query       position in query              matching          repeat                position in repeat
score   div. del. ins.  sequence    begin    end          (left)   repeat            class/family      begin   end    (left)     ID

 1000    5.0  0.0  0.0  chr1            1000     1159    (5000) + fam1              LTR/Gypsy               1    160   (160)     4
  900    5.0  0.0  0.0  chr1            2160     2319    (4000) + fam1              LTR/Gypsy             161    320     (0)     4
"""


class TestFragmentGrouping:
    """Two rows sharing an ID are one copy interrupted by 1,000 bp of genome."""

    @pytest.fixture
    def split(self, tmp_path):
        fasta = tmp_path / "families.fa"
        fasta.write_text(FASTA)
        annot = tmp_path / "split.fa.out"
        annot.write_text(SPLIT)
        return fasta, annot

    def test_rows_sharing_an_id_are_one_copy(self, split, tmp_path, capsys):
        fasta, annot = split
        assert run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path,
                    "--no-dotplot", "--no-orfs", "--no-homology"]) == cli.EXIT_OK
        assert "fam1: 1 copy (2 annotation rows), 1 full length" in capsys.readouterr().out

    def test_no_group_fragments_counts_rows(self, split, tmp_path, capsys):
        fasta, annot = split
        run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path,
             "--no-dotplot", "--no-orfs", "--no-homology", "--no-group-fragments"])
        assert "fam1: 2 copies, 0 full length" in capsys.readouterr().out

    def test_the_sheet_says_how_rows_became_copies(self, split, tmp_path, monkeypatch):
        from teaid import report

        seen = {}
        original = report.build_figure

        def keep(data, theme):
            seen["notes"] = data.notes
            return original(data, theme)

        monkeypatch.setattr(report, "build_figure", keep)
        fasta, annot = split
        run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path,
             "--no-dotplot", "--no-orfs", "--no-homology"])
        assert any("2 annotation rows grouped into 1 copy" in n for n in seen["notes"])


class TestTables:
    """v1's -t and -T, as --tables and --all-tables."""

    def test_tables_are_written_with_headers_and_1_based_coordinates(self, tmp_path, capsys):
        fasta = tmp_path / "families.fa"
        fasta.write_text(FASTA)
        annot = tmp_path / "split.fa.out"
        annot.write_text(SPLIT)
        assert run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path,
                    "--no-dotplot", "--no-orfs", "--no-homology", "--all-tables"]) == cli.EXIT_OK
        names = {p.name for p in tmp_path.glob("fam1.*.tsv")}
        assert names == {"fam1.selfhits.tsv", "fam1.orfs.tsv", "fam1.protein_hits.tsv",
                         "fam1.copies.tsv"}
        rows = (tmp_path / "fam1.copies.tsv").read_text().splitlines()
        assert rows[0].split("\t")[:4] == ["copy", "chrom", "start", "end"]
        first, second = (r.split("\t") for r in rows[1:])
        # The .out rows are 1000-1159 and 2160-2319, 1-based, and one copy.
        assert first[:4] == ["1", "chr1", "1000", "1159"]
        assert second[:4] == ["1", "chr1", "2160", "2319"]
        assert first[-2:] == ["2", "yes"]

    def test_tables_only_when_asked(self, inputs, tmp_path):
        fasta, annot = inputs
        run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "--no-dotplot"])
        assert not list(tmp_path.glob("*.tsv"))


class TestV1ShortFlags:
    """-f and -t mean what they meant in v1 (decided 2026-09-25)."""

    def test_f_is_the_full_length_threshold_and_t_writes_tables(self, tmp_path, capsys):
        fasta = tmp_path / "families.fa"
        fasta.write_text(FASTA)
        annot = tmp_path / "split.fa.out"
        annot.write_text(SPLIT)
        # One copy covering the whole 320 bp consensus: full length at any threshold.
        assert run(["--annot", annot, "-c", fasta, "-n", "fam1", "-f", "1.0", "-t",
                    "-o", tmp_path, *("--no-dotplot", "--no-orfs", "--no-homology")]) == cli.EXIT_OK
        assert "1 full length" in capsys.readouterr().out
        assert (tmp_path / "fam1.orfs.tsv").exists()
        assert not (tmp_path / "fam1.copies.tsv").exists()

    def test_the_family_is_n(self):
        args = cli.build_parser().parse_args(["-n", "fam1", "-f", "0.8"])
        assert args.family == "fam1" and args.full_length_threshold == 0.8


class TestV1DisplayFlags:
    def test_alpha_full_length_alpha_and_auto_y_reach_panel_1(self, tmp_path, monkeypatch):
        from teaid import report

        seen = {}
        original = report.build_figure

        def keep(data, theme):
            seen["fig"] = original(data, theme)
            seen["data"] = data
            return seen["fig"]

        monkeypatch.setattr(report, "build_figure", keep)
        fasta = tmp_path / "families.fa"
        fasta.write_text(FASTA)
        annot = tmp_path / "genome.fa.out"
        annot.write_text(RMOUT)
        assert run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "-a", "0.3",
                    "-F", "1", "-y", "30", "--no-dotplot", "--no-orfs", "--no-homology"]) == cli.EXIT_OK
        opacities = {t.name.split(" ")[0]: t.opacity for t in seen["fig"].data
                     if t.name and t.name.split(" ")[0] in ("partial", "full")}
        assert opacities == {"partial": 0.3, "full": 1.0}
        assert list(seen["fig"].layout.yaxis.range) == [0, 30]

    def test_dotmatcher_is_accepted_and_ignored_with_a_note(self, inputs, tmp_path, capsys):
        fasta, annot = inputs
        assert run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "-D",
                    "--no-dotplot", "--no-orfs", "--no-homology"]) == cli.EXIT_OK
        assert "-D ignored" in capsys.readouterr().err

    def test_an_alpha_outside_0_1_is_a_usage_error(self, inputs, tmp_path):
        fasta, annot = inputs
        with pytest.raises(SystemExit):
            run(["--annot", annot, "-c", fasta, "-n", "fam1", "-o", tmp_path, "-a", "2"])
