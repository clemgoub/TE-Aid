"""The v1 route: copies found by blastn of the consensus against a genome."""

from __future__ import annotations

import random
import shutil

import pytest

from teaid import cli, genome
from teaid.records import Annotation, Copy, DivergenceKind, JunctionKind

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


# Real blastn -outfmt 6 output for a synthetic LTR element (2,800 bp with 400 bp
# LTRs) in a 20 kb chr1 carrying three insertions: a full copy at 3001-5800, a
# solo LTR at 8801-9200, and a copy interrupted by 2 kb of other sequence
# (12201-13400 and 15401-17000). Each LTR also matches the other LTR of the
# consensus, which is where the extra hits come from.
LTR_OUTFMT6 = (
    "ltrX\tchr1\t96.643\t2800\t94\t0\t1\t2800\t3001\t5800\t0.0\t4650\n"
    "ltrX\tchr1\t96.812\t1600\t51\t0\t1201\t2800\t15401\t17000\t0.0\t2673\n"
    "ltrX\tchr1\t96.583\t1200\t41\t0\t1\t1200\t12201\t13400\t0.0\t1989\n"
    "ltrX\tchr1\t97.250\t400\t11\t0\t2401\t2800\t12201\t12600\t0.0\t678\n"
    "ltrX\tchr1\t97.000\t400\t12\t0\t1\t400\t8801\t9200\t0.0\t673\n"
    "ltrX\tchr1\t97.000\t400\t12\t0\t2401\t2800\t8801\t9200\t0.0\t673\n"
    "ltrX\tchr1\t96.535\t404\t14\t0\t1\t404\t5401\t5804\t0.0\t669\n"
    "ltrX\tchr1\t96.750\t400\t13\t0\t1\t400\t16601\t17000\t0.0\t667\n"
    "ltrX\tchr1\t96.000\t400\t16\t0\t2401\t2800\t3001\t3400\t0.0\t651\n"
)


def hsp(start, end, qstart=0, chrom="chr1", strand="+", score=100.0):
    return Copy(chrom=chrom, start=start, end=end, strand=strand, family="f",
                consensus_start=qstart, consensus_end=qstart + end - start, score=score)


class TestRemoveRedundant:
    """-r: one placement per stretch of genome, and interrupted copies joined."""

    def reduce(self, *copies):
        return genome.remove_redundant(Annotation(copies=list(copies)))

    def test_three_ltr_insertions_come_out_as_three_copies(self):
        """Counting every hit gives 9 copies here, and v1's reduce.cpp gave 5."""
        hits = genome.parse(LTR_OUTFMT6, family="ltrX", class_label="LTR/Gypsy",
                            consensus_length=2800)
        reduced, done = genome.remove_redundant(hits)
        assert (done.hits, done.placed, done.copies, done.ties) == (9, 4, 3, 1)
        insertions = reduced.insertions()
        assert [(i.start, i.end, i.n_fragments) for i in insertions] == [
            (3000, 5800, 1), (8800, 9200, 1), (12200, 17000, 2)]
        # The solo LTR matched both LTRs equally; it is drawn at the lower start.
        assert insertions[1].fragments[0].consensus_start == 0
        (join,) = insertions[2].junctions
        assert join.kind is JunctionKind.INSERTION and join.size == 2000

    def test_a_weaker_hit_inside_a_stronger_one_goes(self):
        reduced, done = self.reduce(hsp(0, 1000, score=900), hsp(100, 500, qstart=600, score=300))
        assert [(c.start, c.end) for c in reduced.copies] == [(0, 1000)]
        assert done.ties == 0

    def test_the_outcome_does_not_depend_on_hit_order(self):
        a = hsp(0, 400, qstart=0, score=673)
        b = hsp(0, 400, qstart=2400, score=673)
        for order in ((a, b), (b, a)):
            reduced, done = self.reduce(*order)
            assert [c.consensus_start for c in reduced.copies] == [0]
            assert done.ties == 1

    def test_hits_on_opposite_strands_do_not_compete(self):
        _, done = self.reduce(hsp(0, 1000), hsp(0, 1000, strand="-"))
        assert done.placed == 2

    def test_the_same_locus_on_another_sequence_is_kept(self):
        _, done = self.reduce(hsp(0, 1000), hsp(100, 500, chrom="chr2"))
        assert done.placed == 2

    def test_pieces_in_consensus_order_within_5_kb_are_joined(self):
        reduced, done = self.reduce(hsp(0, 500, qstart=0), hsp(5500, 6000, qstart=500))
        assert done.copies == 1
        assert len({c.hit_id for c in reduced.copies}) == 1

    def test_pieces_further_apart_stay_separate(self):
        _, done = self.reduce(hsp(0, 500, qstart=0), hsp(5501, 6001, qstart=500))
        assert done.copies == 2

    def test_tandem_copies_are_not_joined(self):
        """Two whole copies in a row share most of the consensus, so they are two."""
        _, done = self.reduce(hsp(0, 995, qstart=0), hsp(1200, 2195, qstart=5))
        assert done.copies == 2

    def test_on_the_minus_strand_the_later_piece_sits_earlier_in_the_consensus(self):
        _, joined = self.reduce(hsp(0, 500, qstart=500, strand="-"),
                                hsp(1000, 1500, qstart=0, strand="-"))
        assert joined.copies == 1
        _, apart = self.reduce(hsp(0, 500, qstart=0, strand="-"),
                               hsp(1000, 1500, qstart=500, strand="-"))
        assert apart.copies == 2


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

    def test_every_v1_flag_is_accepted_through_the_alias(self, tmp_path, capsys):
        """The whole v1 flag surface, in one v1-style command line."""
        consensus, genome_fa = _synthetic(tmp_path)
        out = tmp_path / "out"
        code = cli.main_legacy_alias([
            "-q", str(consensus), "-g", str(genome_fa), "-o", str(out), "-t", "-T", "-r",
            "-e", "10e-8", "-f", "0.9", "-m", "300", "-R", "-a", "0.3", "-F", "1", "-y", "10",
            "-D", "--no-homology"])
        assert code == cli.EXIT_OK, capsys.readouterr().err
        assert (out / "famX.copies.tsv").exists() and (out / "famX.teaid.html").exists()

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

    def test_remove_redundant_joins_the_split_copy(self, tmp_path, capsys, monkeypatch):
        from teaid import report

        seen = {}
        original = report.build_figure

        def keep(data, theme):
            seen["data"] = data
            return original(data, theme)

        monkeypatch.setattr(report, "build_figure", keep)
        consensus, genome_fa = _synthetic(tmp_path)
        code = cli.main(["-q", str(consensus), "-g", str(genome_fa), "-r", "-o", str(tmp_path),
                         "--no-dotplot", "--no-orfs", "--no-homology"])
        assert code == cli.EXIT_OK
        # The copy with a 400 bp deletion came back as two hits; -r joins them.
        assert "famX: 3 copies (4 blastn hits kept)" in capsys.readouterr().out
        notes = seen["data"].notes
        assert any(n.startswith("-r: 4 blastn hits, 4 after") and "3 copies" in n for n in notes)
        assert not any("annotation rows" in n for n in notes)

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
        titles = [a.text for a in report.build_figure(data, __import__("teaid").theme.LIGHT)
                  .layout.annotations if a.text and a.text.startswith("1 · ")]
        assert titles and titles[0].startswith("1 · Genomic hits vs identity")


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
