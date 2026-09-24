"""The v1 route: copies found by blastn of the consensus against a genome."""

from __future__ import annotations

import random
import shutil

import pytest

from teaid import cli, genome
from teaid.records import Annotation, Copy, DivergenceKind

# Real blastn -outfmt 6 output for the synthetic genome built below: two whole
# copies on chr1 (one reverse), and one copy on chr2 with a 400 bp internal
# deletion, which blastn reports as two HSPs.
OUTFMT6 = (
    "famX\tchr1\t98.748\t1198\t15\t0\t3\t1200\t3003\t4200\t0.0\t2130\n"
    "famX\tchr1\t95.167\t1200\t58\t0\t1\t1200\t7400\t6201\t0.0\t1895\n"
    "famX\tchr2\t97.211\t502\t14\t0\t1\t502\t2501\t3002\t0.0\t850\n"
    "famX\tchr2\t98.000\t300\t6\t0\t901\t1200\t3001\t3300\t1.52e-150\t521\n"
)


class TestParse:
    def parsed(self):
        return genome.parse(OUTFMT6, family="famX", class_label="DNA/hAT", consensus_length=1200)

    def test_one_copy_per_hsp(self):
        assert len(self.parsed()) == 4

    def test_coordinates_become_half_open_and_strand_comes_from_order(self):
        forward, reverse = self.parsed().copies[:2]
        assert (forward.start, forward.end, forward.strand) == (3002, 4200, "+")
        assert (reverse.start, reverse.end, reverse.strand) == (6200, 7400, "-")
        assert (forward.consensus_start, forward.consensus_end) == (2, 1200)
        assert forward.consensus_left == 0

    def test_identity_is_not_called_divergence(self):
        """v1 plotted 100 - pident as divergence. It is an uncorrected identity."""
        first = self.parsed().copies[0]
        assert first.divergence == pytest.approx(100 - 98.748)
        assert first.divergence_kind is DivergenceKind.BLASTN_IDENTITY
        assert "blastn identity" in DivergenceKind.BLASTN_IDENTITY.axis_label


def hsp(start, end, qstart=0, chrom="chr1"):
    return Copy(chrom=chrom, start=start, end=end, strand="+", family="f",
                consensus_start=qstart, consensus_end=qstart + end - start)


class TestRemoveRedundant:
    """v1's reduce.cpp, ported: a hit inside a kept hit's locus goes."""

    def kept(self, *copies):
        return [(c.chrom, c.start, c.end) for c in genome.remove_redundant(Annotation(copies=list(copies)))]

    def test_a_nested_hit_is_dropped(self):
        assert self.kept(hsp(0, 1000), hsp(100, 500)) == [("chr1", 0, 1000)]

    def test_overlapping_hits_that_do_not_nest_are_both_kept(self):
        assert len(self.kept(hsp(0, 600), hsp(400, 1000))) == 2

    def test_identical_loci_follow_v1s_order_dependent_rule(self):
        """reduce.cpp drops a later identical locus only when the kept one starts
        no later in the query, so the outcome depends on blastn's output order.
        Ported as it is, quirk included."""
        first_smaller = [hsp(0, 500, qstart=10), hsp(0, 500, qstart=90)]
        kept = genome.remove_redundant(Annotation(copies=first_smaller))
        assert [c.consensus_start for c in kept] == [10]
        first_larger = [hsp(0, 500, qstart=90), hsp(0, 500, qstart=10)]
        assert len(genome.remove_redundant(Annotation(copies=first_larger))) == 2

    def test_the_same_locus_on_another_sequence_is_kept(self):
        assert len(self.kept(hsp(0, 1000), hsp(100, 500, chrom="chr2"))) == 2


def _synthetic(tmp_path):
    rng = random.Random(7)
    base = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    cons = base(1200)

    def mutate(s, rate):
        return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in s)

    def rc(s):
        return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

    chr1 = base(3000) + mutate(cons, 0.02) + base(2000) + rc(mutate(cons, 0.05)) + base(3000)
    chr2 = base(2500) + mutate(cons[:500] + cons[900:], 0.03) + base(2500)
    genome_fa = tmp_path / "genome.fa"
    genome_fa.write_text(f">chr1 synthetic\n{chr1}\n>chr2\n{chr2}\n")
    consensus = tmp_path / "cons.fa"
    consensus.write_text(f">famX#DNA/hAT\n{cons}\n")
    return consensus, genome_fa


needs_blast = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="BLAST+ not installed",
)


@needs_blast
class TestTheV1CommandLine:
    def test_a_v1_command_line_runs(self, tmp_path, capsys):
        consensus, genome_fa = _synthetic(tmp_path)
        code = cli.main(["-q", str(consensus), "-g", str(genome_fa), "-o", str(tmp_path),
                         "-m", "400", "-R", "--no-dotplot", "--no-homology"])
        out = capsys.readouterr()
        assert code == cli.EXIT_OK, out.err
        # Four HSPs, four copies: the chr2 copy is split at its deletion, which
        # is the leakiness the sheet warns about.
        assert "famX: 4 copies" in out.out

    def test_the_database_is_built_beside_the_genome_as_v1_did(self, tmp_path):
        consensus, genome_fa = _synthetic(tmp_path)
        cli.main(["-q", str(consensus), "-g", str(genome_fa), "-o", str(tmp_path),
                  "--no-dotplot", "--no-orfs", "--no-homology"])
        assert genome._is_database(genome_fa)

    def test_the_blastn_spelling_is_the_same_route(self, tmp_path, capsys):
        consensus, genome_fa = _synthetic(tmp_path)
        cli.main(["-c", str(consensus), "--blastn", str(genome_fa), "-o", str(tmp_path),
                  "--no-dotplot", "--no-orfs", "--no-homology"])
        assert "famX: 4 copies" in capsys.readouterr().out

    def test_remove_redundant_is_accepted(self, tmp_path, capsys):
        consensus, genome_fa = _synthetic(tmp_path)
        code = cli.main(["-q", str(consensus), "-g", str(genome_fa), "-r", "-o", str(tmp_path),
                         "--no-dotplot", "--no-orfs", "--no-homology"])
        assert code == cli.EXIT_OK
        assert "famX: 4 copies" in capsys.readouterr().out  # none of the four nest

    def test_the_sheet_says_how_its_copies_were_found(self, tmp_path, monkeypatch):
        from teaid import report

        seen = {}
        original = report.build_figure

        def keep(data, theme):
            seen["data"] = data
            return original(data, theme)

        monkeypatch.setattr(report, "build_figure", keep)
        consensus, genome_fa = _synthetic(tmp_path)
        cli.main(["-q", str(consensus), "-g", str(genome_fa), "-o", str(tmp_path),
                  "--no-dotplot", "--no-orfs", "--no-homology"])
        data = seen["data"]
        assert any("genome.fa" in s and "blastn" in s for s in data.sources)
        assert any("misses diverged copies" in n for n in data.notes)
        assert data.divergence_kind is DivergenceKind.BLASTN_IDENTITY


class TestGenomeRouteFailures:
    def test_a_missing_genome_is_no_input(self, tmp_path, capsys):
        consensus, _ = _synthetic(tmp_path)
        code = cli.main(["-q", str(consensus), "-g", str(tmp_path / "absent.fa"), "-o", str(tmp_path)])
        assert code == cli.EXIT_NO_INPUT
        assert "[no-input]" in capsys.readouterr().err

    def test_the_genome_route_needs_a_consensus(self, tmp_path, capsys):
        _, genome_fa = _synthetic(tmp_path)
        assert cli.main(["-g", str(genome_fa), "-o", str(tmp_path)]) == cli.EXIT_USAGE
        assert "--consensus is required with --genome" in capsys.readouterr().err
