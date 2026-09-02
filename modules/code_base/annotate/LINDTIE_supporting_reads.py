'''
Module      : LINDTIE_supporting_reads
Description : Counts, per variant, the case reads that carry the variant
              (supporting_read_count) and the case reads that span its breakpoint
              (spanning_reads), directly from the genome alignment.

Rationale   : oarfish num_reads is whole-contig abundance and overstates variant support.
              This version counts support from the ALIGNMENT itself -- the same alignment
              IGV displays -- so the numbers match IGV by construction. A read supports the
              variant if, at the breakpoint, its alignment reproduces the variant's
              signature (the same thing the contig does):

                RI          : aligns THROUGH the exon->intron boundary, staying aligned
                              min(intron, through_depth) bp into the retained intron
                EE          : aligns THROUGH the boundary between the reference exon and the
                              NOVEL extension of the block, min(extension, through_depth) bp
                              into the extension (the block midpoint is NOT discriminating --
                              it usually sits inside the reference exon)
                NE          : splices INTO the novel exon at the junction
                AS/NEJ/PNJ  : has a gap whose start ~pos1, end ~pos2 AND whose deleted length
                              matches the variant, and that is not better explained by an
                              annotated (canonical) junction
                DEL         : same, read from the RAW cigar so few-bp deletions register
                FUS, IGR    : is a split read (SA) bridging pos1(chrom1) and pos2(chrom2)
                UN          : soft/hard-clips at ~the breakpoint
                INS, ITD,PTD: carries ~varsize inserted bases at ~pos1 (reliability=low
                              -- ONT insertions align unreliably)
                FUS/IGR/INS/small-DEL are reported supporting_read_count_reliability=low.

              Accuracy comes from three structural requirements, not from tuning:
                1. the read's feature must be in the right PLACE (position tolerance),
                2. it must be the right SIZE (a 17 bp deletion is not a 1 bp ONT indel, a
                   1252 bp junction is not the 1249 bp canonical one),
                3. it must not be better explained by a COMPETING reference structure (an
                   annotated intron / a reference exon), which is what separates a variant
                   junction from the canonical junction it sits a few bp away from.

              Requirements 2 and 3 need the reference annotation (--tx_annotation, the same
              GTF the annotation step uses). Without it the module falls back to the
              position-only behaviour and marks the affected rows reliability=low.

              spanning_reads (the denominator) is the alignment envelope over the counted
              junction (reference_start <= pos <= reference_end), unioned with the supporting
              set -- always defined, and it is what IGV shows as coverage.

Coordinates : pos1/pos2 in the variant table are 1-BASED; pysam blocks and CIGAR offsets are
              0-BASED half-open. The conventions below were verified against the alignments
              (a read carrying the variant matches with offset 0):

                AS/NEJ/PNJ  gap = (pos1, pos2 - 1)   deleted length = pos2 - pos1 - 1
                DEL         gap = (pos1, pos2 + 1)   deleted length = pos2 - pos1 + 1
                RI/EE/NE    block = [pos1 - 1, pos2) (0-based half-open)

Copyright   : (c) Jia Wei Tan, 2026
License     : MIT
'''

import logging
import re
import sys
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import pysam

# refined_annotated_contigs_info.tsv stores pos1/pos2 as "chr:pos(strand)" e.g.
# "chr1:945043(-)". parse_locus() also accepts a bare integer position.
_LOCUS_RE = re.compile(r'^\s*([^:\s]+):(\d+)(?:\(([+-])\))?')
_TX_RE = re.compile(r'transcript_id "([^"]+)"')
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

# Variant classes -> the alignment signature a supporting read must reproduce.
THROUGH_TYPES = {'RI'}                  # retained intron: reads align contiguously across
INSIDE_TYPES = {'NE'}                  # novel exon: reads splice INTO the segment
INTERIOR_TYPES = {'EE'}                # extended exon: reads align THROUGH the novel extension
GAP_TYPES = {'AS', 'NEJ', 'PNJ', 'DEL'}  # a splice/deletion gap pos1<->pos2
SPLIT_TYPES = {'FUS', 'IGR'}            # split read bridging two loci
INSERTION_TYPES = {'INS', 'ITD', 'PTD'}  # extra sequence -- low confidence (ONT indels)
CLIP_TYPES = {'UN'}                     # soft/hard-clip at the breakpoint

CIGAR_INS, CIGAR_DEL, CIGAR_SKIP = 1, 2, 3
CIGAR_CLIPS = (4, 5)                    # soft-clip, hard-clip
CIGAR_REF_CONSUMING = (0, 2, 3, 7, 8)   # M D N = X

# Below this size a gap is a deletion, whose size the assembly and the reads often disagree
# about; at or above it the gap is a splice junction, which both place identically. sig_gap()
# matches the two regimes differently -- see its docstring.
SMALL_EVENT_BP = 200

# Smallest insertion this module will call, in bp. ONT indel noise lives below it, so an
# insertion smaller than this cannot be supported by an alignment signature at all.
MIN_INSERTION_BP = 15

OUT_COLUMNS = ['variant_id', 'contig_id', 'variant_type', 'supporting_read_count',
               'spanning_reads', 'junction_VAF', 'candidate_depth', 'counted_locus',
               'support_signature', 'supporting_read_count_reliability']


def parse_args():
    p = ArgumentParser(description='Count alignment-verified variant-supporting reads and '
                                   'breakpoint-spanning reads in the case sample.')
    p.add_argument('variant_info',
                   help='refined_annotated_contigs_info.tsv (contig_id, variant_id, '
                        'variant_type, pos1, pos2, varsize, cpos, contig_varsize)')
    p.add_argument('case_bam', help='case genome BAM (reads_all_sorted.bam), indexed')
    p.add_argument('--tx_annotation', default=None,
                   help='transcriptome GTF (the same reference the annotation step uses). '
                        'Supplies the reference exons that define which part of a block is '
                        'NOVEL (EE/RI) and the annotated junctions a read may be better '
                        'explained by (AS/NEJ/PNJ/DEL). Without it those rules fall back '
                        'to position-only matching and are reported reliability=low')
    p.add_argument('--window', type=int, default=300,
                   help='bp either side of a breakpoint to fetch candidate reads (default 300)')
    p.add_argument('--min_anchor', type=int, default=20,
                   help='min aligned bp required on the REFERENCE side of a boundary for '
                        'THROUGH / INTERIOR / INSIDE (default 20)')
    p.add_argument('--through_depth', type=int, default=100,
                   help='how far a read must stay aligned INTO the novel segment (retained '
                        'intron / exon extension / novel exon) to count as support: '
                        'min(segment length, this). Rejects reads that merely spill a few bp '
                        'past an exon edge (default 100)')
    p.add_argument('--min_gap', type=int, default=25,
                   help='alignment gaps <= this many bp are merged into the adjacent block, so '
                        'ONT indel noise does not fragment reads; larger gaps are splice/'
                        'deletion junctions (default 25)')
    p.add_argument('--gap_merge', type=int, default=30,
                   help='when matching a read to a splice/deletion gap, consecutive D/N ops '
                        'separated by <= this many aligned bp are treated as ONE event, so a '
                        'junction plus an adjacent small deletion (and an ONT-fragmented '
                        'deletion) match the single gap the contig shows (default 30)')
    p.add_argument('--tol', type=int, default=35,
                   help='max breakpoint position tolerance in bp when matching a read gap/clip/'
                        'insertion to the annotated junction. For a gap it is tightened to '
                        'min(tol, max(5, 10%% of the gap)) so a small deletion cannot absorb '
                        'unrelated indel noise (default 35)')
    p.add_argument('--bp_tol', type=int, default=10,
                   help='breakpoint tolerance for a clip position (UN); the contig gives an '
                        'exact breakpoint, so this is tighter than --tol (default 10)')
    p.add_argument('--near_canonical', type=int, default=10,
                   help='if an annotated junction sits within this many bp of a variant '
                        'junction (but is not identical), the variant is a shifted canonical '
                        'junction -- reported reliability=low (default 10)')
    p.add_argument('--min_clip', type=int, default=25,
                   help='min soft/hard-clip length at the breakpoint for UN support (default 25)')
    p.add_argument('--inside_window', type=int, default=250,
                   help='for NE/EE, how far past the counted junction to look for spliced-in '
                        'reads (default 250)')
    p.add_argument('--ins_tol', type=int, default=100,
                   help='breakpoint tolerance for INS: inserted bases within this many bp of '
                        'pos1 are pooled (ONT insertions are placed at scattered positions; '
                        'default 100)')
    p.add_argument('--split_tol', type=int, default=20000,
                   help='max bp between a split read\'s supplementary alignment (either end) '
                        'and pos2 for it to count as bridging a fusion. Fusion partners can be '
                        'far apart and tightening this only cost recall on the simulation '
                        '(default 20000)')
    p.add_argument('--support_reads_out', default=None,
                   help='optional TSV of variant_id -> supporting read name')
    p.add_argument('--out', required=True, help='output TSV')
    p.add_argument('--log', default=None, help='record program progress in LOG_FILE')
    # Accepted for backward compatibility with the old process command line; unused now.
    for dead in ('--contig_fasta', '--contig_bam', '--genome_fasta', '--flank',
                 '--max_edit_frac', '--both_ends_max_size', '--both_ends_min_depth',
                 '--max_pool', '--max_tandem', '--min_entropy', '--probes_out'):
        p.add_argument(dead, default=None, help='(deprecated, ignored)')
    return p.parse_args()


def init_logging(log_filename):
    fmt = '%(asctime)s %(levelname)s - %(message)s'
    if log_filename:
        logging.basicConfig(filename=log_filename, level=logging.INFO, filemode='w',
                            format=fmt, datefmt="%Y-%m-%dT%H:%M:%S%z")
    else:
        logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                            format=fmt, datefmt="%Y-%m-%dT%H:%M:%S%z")
    logging.info('program started: %s', ' '.join(sys.argv))


def parse_locus(field):
    '''(chrom, pos, strand) from a "chr:pos(strand)" string, else (None, None, None).'''
    if field is None:
        return (None, None, None)
    m = _LOCUS_RE.match(str(field).strip())
    if not m:
        return (None, None, None)
    return (m.group(1), int(m.group(2)), m.group(3))


#------------------------------------------------------------------------------
# Reference annotation: merged exons + annotated junctions
#------------------------------------------------------------------------------

class RefAnnotation:
    '''Merged reference exons and annotated introns per chromosome, in 0-based half-open
    coordinates -- the same exon intervals the annotation step derives (get_gene_lookup /
    block_helper.is_novel_block), so "novel" means the same thing here as it does there.

    Two questions are answered:
      * which part of a contig block is NOVEL (not covered by any reference exon) -- the only
        part where a variant read can differ from a wild-type one (EE/RI);
      * is a read's gap better explained by an ANNOTATED junction than by the variant
        (AS/NEJ/PNJ/DEL) -- i.e. is it just a wild-type read.
    '''

    def __init__(self, gtf_path):
        self.exons = {}     # chrom -> (starts, ends)      merged, sorted, 0-based half-open
        self.introns = {}   # chrom -> (donors, acceptors) sorted by donor, 0-based
        self._load(gtf_path)

    def _load(self, gtf_path):
        logging.info('Loading reference exons/junctions from %s', gtf_path)
        gtf = pd.read_csv(gtf_path, comment='#', sep='\t', header=None, low_memory=False,
                          usecols=[0, 2, 3, 4, 8],
                          names=['chrom', 'feature', 'start', 'end', 'attr'])
        ex = gtf[gtf.feature == 'exon'].copy()
        del gtf
        ex['start'] = ex.start.values - 1          # GTF is 1-based inclusive -> 0-based
        ex['tx'] = ex.attr.str.extract(_TX_RE, expand=False)

        for chrom, grp in ex.groupby('chrom', sort=False):
            grp = grp.sort_values('start')
            starts, ends = grp.start.values, grp.end.values
            reach = np.maximum.accumulate(ends)
            new = np.r_[True, starts[1:] > reach[:-1]]
            self.exons[chrom] = (starts[new],
                                 np.maximum.reduceat(ends, np.flatnonzero(new)))

        # Annotated introns: the gap between consecutive exons of the same transcript.
        ex = ex.sort_values(['tx', 'start'])
        tx, chrom = ex.tx.values, ex.chrom.values
        starts, ends = ex.start.values, ex.end.values
        is_intron = (tx[:-1] == tx[1:]) & (starts[1:] > ends[:-1])
        for c in np.unique(chrom[:-1][is_intron]):
            keep = is_intron & (chrom[:-1] == c)
            # Transcripts of a gene share introns; keep each junction once so that
            # introns_near() returns a handful of distinct rivals, not hundreds of copies.
            pairs = np.unique(np.stack((ends[:-1][keep], starts[1:][keep]), axis=1), axis=0)
            order = np.lexsort((pairs[:, 1], pairs[:, 0]))
            self.introns[c] = (pairs[order, 0], pairs[order, 1])

        logging.info('Loaded %d merged exons and %d distinct annotated junctions on %d contigs.',
                     sum(len(s) for s, _ in self.exons.values()),
                     sum(len(d) for d, _ in self.introns.values()), len(self.exons))

    def exons_overlapping(self, chrom, start, end):
        '''Merged reference exons overlapping [start, end), in order.'''
        if chrom not in self.exons or end <= start:
            return []
        starts, ends = self.exons[chrom]
        i = int(np.searchsorted(ends, start, 'right'))
        out = []
        while i < len(starts) and starts[i] < end:
            out.append((int(starts[i]), int(ends[i])))
            i += 1
        return out

    def is_exonic(self, chrom, pos):
        return bool(self.exons_overlapping(chrom, pos, pos + 1))

    def novel_intervals(self, chrom, start, end):
        '''[start, end) minus the reference exons -- the parts of a block that are novel.'''
        out, cur = [], start
        for s, e in self.exons_overlapping(chrom, start, end):
            if s > cur:
                out.append((cur, min(s, end)))
            cur = max(cur, e)
        if cur < end:
            out.append((cur, end))
        return out

    def introns_near(self, chrom, donor, acceptor, tol):
        '''Annotated introns whose donor AND acceptor are within tol of (donor, acceptor) --
        the wild-type structures a read at this locus might really be showing.'''
        if chrom not in self.introns:
            return []
        d, a = self.introns[chrom]
        lo = int(np.searchsorted(d, donor - tol, 'left'))
        hi = int(np.searchsorted(d, donor + tol, 'right'))
        return [(int(d[i]), int(a[i])) for i in range(lo, hi)
                if abs(int(a[i]) - acceptor) <= tol]


#------------------------------------------------------------------------------
# Alignment helpers
#------------------------------------------------------------------------------

def merged_blocks(read, min_gap):
    '''Aligned segments of a read in genome space, merging gaps <= min_gap so small ONT
    indels do not split a read; gaps larger than min_gap (splice / deletion) remain as
    block boundaries. Returns [(start, end), ...].'''
    bl = read.get_blocks()
    if not bl:
        return []
    m = [list(bl[0])]
    for s, e in bl[1:]:
        if s - m[-1][1] <= min_gap:
            m[-1][1] = e
        else:
            m.append([s, e])
    return [(s, e) for s, e in m]


def candidate_gaps(read, gap_merge):
    '''Every deletion/skip a read could be showing at a locus: the RAW D/N ops (so a few-bp
    deletion registers) AND their COMPOSITES -- consecutive D/N ops separated by <= gap_merge
    aligned bp merged into one event. The composite is what matches a contig junction that
    carries an extra small deletion at its edge, and it rejoins an ONT-fragmented deletion.
    Returns {(start, end, deleted_bases)}.'''
    rpos, raw = read.reference_start, []
    for op, ln in (read.cigartuples or []):
        if op in (CIGAR_DEL, CIGAR_SKIP):
            raw.append((rpos, rpos + ln, ln))
        if op in CIGAR_REF_CONSUMING:
            rpos += ln
    if not raw:
        return ()
    out = set(raw)
    comp = [list(raw[0])]
    for s, e, ln in raw[1:]:
        if s - comp[-1][1] <= gap_merge:
            comp[-1][1], comp[-1][2] = e, comp[-1][2] + ln
        else:
            comp.append([s, e, ln])
    out.update((s, e, ln) for s, e, ln in comp)
    return out


def sa_alignments(read):
    '''(chrom, start, end) of each supplementary alignment in the SA tag, 0-based.'''
    if not read.has_tag('SA'):
        return []
    out = []
    for sa in str(read.get_tag('SA')).split(';'):
        f = sa.split(',')
        if len(f) < 4:
            continue
        try:
            start = int(f[1]) - 1
        except ValueError:
            continue
        span = sum(int(n) for n, op in _CIGAR_RE.findall(f[3]) if op in 'MDN=X')
        out.append((f[0], start, start + span))
    return out


def fetch_reads(case_bam, chrom, lo, hi):
    if chrom is None:
        return []
    try:
        return [r for r in case_bam.fetch(chrom, max(0, lo), hi) if not r.is_unmapped]
    except (ValueError, KeyError):
        return []


def junction_read_stats(case_bam, chrom, pos, window):
    '''One cheap pass: (depth_names, span_names). depth = distinct reads within +/-window;
    span = reads whose alignment envelope covers the breakpoint. Signature-independent, so
    spanning is always defined (matches IGV coverage).'''
    depth, span = set(), set()
    for r in fetch_reads(case_bam, chrom, pos - window, pos + window):
        depth.add(r.query_name)
        if r.reference_start <= pos <= r.reference_end:
            span.add(r.query_name)
    return depth, span


def gap_tolerances(length, tol):
    '''(position tolerance, size tolerance) for matching a read gap to a variant gap of
    `length` bp. The position tolerance shrinks with the event -- for a 17 bp deletion the
    default 35 bp is wider than the deletion itself, which is how ONT 1-2 bp indel noise used
    to be counted as support. The size tolerance is what stops a 1 bp indel standing in for a
    17 bp deletion, and a 1249 bp canonical intron standing in for a 1251 bp variant one.

    The constants are the best precision/recall trade-off on the badread simulation: tighter
    (max(5, 0.1L) / max(2, 0.05L)) costs 8% of the true DEL reads for +0.5 pp precision, and
    looser (position 35, or size max(5, 0.15L)) drops DEL precision from 98.0% to 96.9% / 95.1%.'''
    return (min(tol, max(10, int(0.2 * length))), min(25, max(3, int(0.1 * length))))


def expected_gap(variant_type, pos1, pos2):
    '''The gap a supporting read must show, in 0-based half-open genome coordinates.
    See the coordinate conventions in the module docstring.'''
    lo, hi = min(pos1, pos2), max(pos1, pos2)
    return (lo, hi + 1) if variant_type == 'DEL' else (lo, hi - 1)


#------------------------------------------------------------------------------
# Per-type supporting-read signatures. Each returns a set of supporting read names.
#------------------------------------------------------------------------------

def read_through_depth(seg_len, args):
    '''How far INTO the novel segment (retained intron / exon extension) a supporting read must
    stay aligned: the whole segment, capped at `through_depth`, less `bp_tol`.

    The subtraction is not cosmetic. The far end of the segment is itself a junction -- the donor
    into the next intron -- and the aligner wobbles there by a base or two. Demanding the segment's
    LAST base means one bp of wobble rejects every read: at PSCA the reads' blocks all ended at
    142681486 and the test wanted 142681487, so 41 of 43 variant reads failed and the variant
    reported 2.'''
    if not seg_len:
        return args.min_anchor
    return max(args.min_anchor, min(seg_len - args.bp_tol, args.through_depth))


def sig_through(case_bam, chrom, j, novel_right, ref_anchor, depth, min_gap, window):
    '''Read aligns straight through the boundary j between reference sequence and a novel
    segment: one merged block covers `ref_anchor` bp on the reference side and `depth` bp into
    the novel segment. Used for a retained intron (exon->intron boundary) and for an extended
    exon (exon->extension boundary). A wild-type read splices AT j, so no block covers it.'''
    left, right = (ref_anchor, depth) if novel_right else (depth, ref_anchor)
    return {r.query_name for r in fetch_reads(case_bam, chrom, j - window, j + window)
            if any(s <= j - left and e >= j + right for s, e in merged_blocks(r, min_gap))}


def _splice_in_side(case_bam, chrom, j, seg_right, window, anchor, tol, min_gap):
    '''Reads that SPLICE IN at j on the given side: a merged block whose near edge sits within
    `tol` of j and that extends >= anchor aligned bases into the segment. NE is a splice-IN
    variant -- the supporting read's aligned block BEGINS (acceptor) or ENDS (donor) AT the
    junction, it does not align through it -- so the signal is a block EDGE at j.

    Unchanged by the accuracy fix, deliberately: on the simulation NE is already 98.7% precise,
    and neither requiring a deeper read-through into the novel exon nor rejecting reads whose
    edge sits nearer an annotated exon edge improved it (98.8% at best), while both cost 4-5%
    of the true NE reads.'''
    lo, hi = (j, j + window) if seg_right else (j - window, j)
    out = set()
    for r in fetch_reads(case_bam, chrom, lo, hi):
        for b0, b1 in merged_blocks(r, min_gap):
            if seg_right and b0 <= j + tol and b1 >= j + anchor:
                out.add(r.query_name)
                break
            if not seg_right and b1 >= j - tol and b0 <= j - anchor:
                out.add(r.query_name)
                break
    return out


def sig_inside(case_bam, chrom, j, other, window, anchor, tol, min_gap):
    '''Novel exon with a known segment side: a read supports iff it splices in at j on the
    side toward `other`. Measured at the single counted junction only (measuring both
    boundaries would union two disjoint read sets for a large exon and inflate the count).'''
    return _splice_in_side(case_bam, chrom, j, other > j, window, anchor, tol, min_gap)


def sig_inside_dominant(case_bam, chrom, j, window, anchor, tol, min_gap):
    '''Single-point NE/EE (pos1 == pos2): the segment side is not encoded, so take whichever
    side carries the splice-in reads (the extension). A canonical read splices AT j with no
    aligned bases on the segment side.'''
    right = _splice_in_side(case_bam, chrom, j, True, window, anchor, tol, min_gap)
    left = _splice_in_side(case_bam, chrom, j, False, window, anchor, tol, min_gap)
    return right if len(right) >= len(left) else left


def sig_gap(case_bam, chrom, s_exp, e_exp, variant_type, tol, window, gap_merge, rivals):
    '''Read jumps s_exp -> e_exp the way the variant does, and is not better explained by one of
    `rivals` (the annotated junctions at this locus). How closely the read's gap must match the
    variant's depends on how big the event is, because the two error regimes are opposite:

    LARGE events (>= SMALL_EVENT_BP) -- splice junctions. The aligner snaps these to the same
    coordinates in the reads as in the contig, so both endpoints must land within `pos_tol` and
    the deleted length within `size_tol`. That strictness is what stops the canonical junction
    being counted when the variant is that junction shifted by a couple of bp.

    SMALL DELETIONS (< SMALL_EVENT_BP) -- here the ASSEMBLY and the reads routinely disagree about
    the size of the same event: the contig called 27 bp where the reads show 45, and 19 bp where
    the reads show 38 (ground truth confirms the reads). Demanding the contig's size discards the
    real evidence, so a read matches if its gap OVERLAPS the variant, is 0.5-2x its size, and is
    anchored at one end. ONT 1-2 bp indel noise still fails the ratio band (1 bp against a 17 bp
    deletion is 0.06x), which is the over-count this originally fixed.

    The loose branch is for DEL only. Applying it to short splice junctions too cost 13 extra
    over-counted PNJ on the simulation: a spliced read's coordinates are snapped by the aligner,
    so there is no mis-sizing to forgive, and the looseness only lets neighbouring junctions in.'''
    length = e_exp - s_exp
    pos_tol, size_tol = gap_tolerances(length, tol)
    small = variant_type == 'DEL' and length < SMALL_EVENT_BP
    lo, hi = min(s_exp, e_exp), max(s_exp, e_exp)
    out = set()
    for r in fetch_reads(case_bam, chrom, lo - window, hi + window):
        for gs, ge, glen in candidate_gaps(r, gap_merge):
            if small:
                if not 0.5 * length <= glen <= 2 * length:
                    continue
                if min(ge, e_exp) - max(gs, s_exp) <= 0:      # must overlap the variant itself
                    continue
                if abs(gs - s_exp) > pos_tol and abs(ge - e_exp) > pos_tol:
                    continue                                  # anchored at neither end
            else:
                if abs(gs - s_exp) > pos_tol or abs(ge - e_exp) > pos_tol:
                    continue
                if abs(glen - length) > size_tol:
                    continue
            d_var = abs(gs - s_exp) + abs(ge - e_exp)
            if any(abs(gs - a) + abs(ge - b) <= d_var for a, b in rivals):
                continue    # the read matches a canonical junction at least as well
            out.add(r.query_name)
            break
    return out


def sig_split(case_bam, chrom1, p1, chrom2, p2, window, split_tol):
    '''Read is chimeric: its primary alignment covers p1 and a supplementary alignment (SA tag)
    reaches p2 on chrom2 -- a split read bridging the two fusion partners. EITHER END of the
    supplementary alignment may sit at p2, depending on the orientation of the join, which the
    earlier start-only test missed (+2% true fusion reads on the simulation).

    Also requiring the primary to clip at p1 was tried and rejected: on the simulation it cost
    28% of the true fusion reads at identical precision (99.4%).'''
    out = set()
    for r in fetch_reads(case_bam, chrom1, p1 - window, p1 + window):
        for sa_chrom, sa_start, sa_end in sa_alignments(r):
            if sa_chrom == chrom2 and min(abs(sa_start - p2), abs(sa_end - p2)) <= split_tol:
                out.add(r.query_name)
                break
    return out


def sig_clip(case_bam, chrom, j, min_clip, tol, window):
    '''Read soft/hard-clips >= min_clip right at j (a breakpoint where it diverges).'''
    out = set()
    for r in fetch_reads(case_bam, chrom, j - window, j + window):
        ct = r.cigartuples
        if not ct:
            continue
        if ct[0][0] in CIGAR_CLIPS and ct[0][1] >= min_clip and abs(r.reference_start - j) <= tol:
            out.add(r.query_name)
        elif ct[-1][0] in CIGAR_CLIPS and ct[-1][1] >= min_clip and abs(r.reference_end - j) <= tol:
            out.add(r.query_name)
    return out


def sig_insertion(case_bam, chrom, j, size, tol, window):
    '''Read carries the insertion: the inserted bases within +/-tol of j sum to ~size (ONT
    aligners scatter one insertion across several I ops, so they are pooled rather than
    matched individually), OR the read carries a large soft-clip near j -- the aligner
    frequently clips an insertion it cannot place. Low confidence either way.

    The floor is MIN_INSERTION_BP, NOT min_clip. min_clip is a soft-CLIP threshold (25 bp);
    using it here meant a 20 bp insertion needed 25 inserted bases and could never be called at
    all -- at FBLN5 that reported 1 supporting read where 41 of the 42 variant reads carry a
    17-24 bp insertion. The two thresholds measure different things and are kept apart.

    Pooling alone is not enough to call a read, though: `ins_tol` spans 200 bp, and a read with a
    dozen 1-3 bp ONT hiccups across that window reaches the floor without carrying an insertion at
    all. Dropping the floor to 15 bp on its own cost 20 points of INS precision on the simulation.
    So a read must also show ONE substantial op -- the insertion may be fragmented, but not
    entirely into noise.'''
    need_lo = max(MIN_INSERTION_BP, int(0.5 * size)) if size else MIN_INSERTION_BP
    need_hi = max(2 * size, 2 * MIN_INSERTION_BP) if size else None
    need_single = max(10, int(0.25 * size)) if size else 10
    clip_need = max(MIN_INSERTION_BP, 0.4 * size) if size else MIN_INSERTION_BP
    out = set()
    for r in fetch_reads(case_bam, chrom, j - window, j + window):
        rpos, inserted, biggest = r.reference_start, 0, 0
        for op, ln in (r.cigartuples or []):
            if op == CIGAR_INS and abs(rpos - j) <= tol:
                inserted += ln
                biggest = max(biggest, ln)
            if op in CIGAR_REF_CONSUMING:
                rpos += ln
        hit = (inserted >= need_lo and biggest >= need_single
               and (need_hi is None or inserted <= need_hi))
        if not hit:
            ct = r.cigartuples or []
            if ct and ((ct[0][0] in CIGAR_CLIPS and ct[0][1] >= clip_need
                        and abs(r.reference_start - j) <= tol)
                       or (ct[-1][0] in CIGAR_CLIPS and ct[-1][1] >= clip_need
                           and abs(r.reference_end - j) <= tol)):
                hit = True
        if hit:
            out.add(r.query_name)
    return out


#------------------------------------------------------------------------------
# Boundary selection for the block types (RI / EE)
#------------------------------------------------------------------------------

def block_boundaries(ref, chrom, block_lo, block_hi, variant_type):
    '''Where a block-type variant is discriminating, as [(j, novel_right, novel_len), ...].

    A block spans [block_lo, block_hi) (0-based). The parts of it NOT covered by a reference
    exon are the novel ones -- for an EE the extension, for an RI the retained intron. A
    supporting read stays aligned from the reference side of the exon/novel boundary into the
    novel part. Both ENDS of the block sit at canonical edges shared with wild-type, and (for
    a one-sided EE) its midpoint usually sits INSIDE the reference exon, so neither of those
    is discriminating.

    Returns [] when the annotation cannot resolve a novel part (caller falls back).'''
    if ref is None or block_hi <= block_lo:
        return []
    novel = ref.novel_intervals(chrom, block_lo, block_hi)
    if not novel:
        return []
    n1, n2 = max(novel, key=lambda iv: iv[1] - iv[0])
    length = n2 - n1
    out = []
    # A boundary is usable when its reference side is exonic: either an exon inside the block
    # (a one-sided EE) or the exon the block abuts (an RI, or a block flush with an exon edge).
    if n1 > block_lo or ref.is_exonic(chrom, n1 - 1):
        out.append((n1, True, length))
    if n2 < block_hi or ref.is_exonic(chrom, n2):
        out.append((n2, False, length))
    if not out and variant_type in THROUGH_TYPES:
        out = [(n1, True, length), (n2, False, length)]
    return out


#------------------------------------------------------------------------------

def count_variant(row, case_bam, args, ref=None, support_writer=None, notes=None):
    vid = row['variant_id']
    vt = row['variant_type']
    rec = {'variant_id': vid, 'contig_id': row['contig_id'], 'variant_type': vt}

    chrom1, pos1, _ = parse_locus(row.get('pos1'))
    chrom2, pos2, _ = parse_locus(row.get('pos2'))
    if chrom1 is None or pos1 is None:
        rec.update({'supporting_read_count': pd.NA, 'spanning_reads': pd.NA,
                    'junction_VAF': pd.NA, 'candidate_depth': pd.NA, 'counted_locus': '',
                    'support_signature': 'none',
                    'supporting_read_count_reliability': 'low'})
        return rec
    try:
        varsize = abs(int(float(row.get('varsize')))) if pd.notna(row.get('varsize')) else 0
    except (TypeError, ValueError):
        varsize = 0

    low_conf = False
    two_point = (chrom2 == chrom1 and pos2 is not None and pos2 != pos1)
    seg_len, novel_right = 0, True

    # --- Block types (RI/EE): count at the boundary between reference exon and novel segment.
    boundaries = []
    if vt in THROUGH_TYPES | INTERIOR_TYPES and two_point:
        block_lo, block_hi = min(pos1, pos2) - 1, max(pos1, pos2)
        boundaries = block_boundaries(ref, chrom1, block_lo, block_hi, vt)
        if not boundaries and vt in INTERIOR_TYPES:
            # No resolvable novel part (no annotation, or the block is entirely exonic): fall
            # back to the old block-midpoint test, which cannot separate the variant from
            # wild-type when the midpoint sits inside a reference exon.
            boundaries = [((block_lo + block_hi) // 2, True,
                           max(1, (block_hi - block_lo) // 2))]
            low_conf = True
            if notes is not None:
                notes['ee_no_novel_region'] += 1

    if boundaries:
        # Two usable boundaries (a block flanked by exons) carry disjoint read sets for a long
        # segment, so count at the better-covered one -- never sum.
        best = None
        for j, side_right, length in boundaries:
            depth_names, span_names = junction_read_stats(case_bam, chrom1, j, args.window)
            if best is None or len(span_names) > len(best[4]):
                best = (j, side_right, length, depth_names, span_names)
        p_star, novel_right, seg_len, depth_names, span_names = best
        c_star = chrom1
        depth = len(depth_names)
        counted_locus = '%s:%d' % (c_star, p_star)
    elif vt in THROUGH_TYPES | INTERIOR_TYPES:
        # Single-point EE/RI, or a block whose partner locus is missing: keep pos1.
        c_star, p_star = chrom1, pos1
        depth_names, span_names = junction_read_stats(case_bam, c_star, p_star, args.window)
        depth = len(depth_names)
        counted_locus = '%s:%d' % (c_star, p_star)
        if vt in THROUGH_TYPES and two_point:
            seg_len = abs(pos2 - pos1)
            novel_right = pos2 > pos1
    else:
        junctions = [(chrom1, pos1)]
        if chrom2 is not None and (chrom2 != chrom1
                                   or abs((pos2 or 0) - (pos1 or 0)) > args.window):
            junctions.append((chrom2, pos2))
        per_j = [(c, p) + junction_read_stats(case_bam, c, p, args.window) for c, p in junctions]

        if vt in SPLIT_TYPES or len(per_j) == 1:
            # Fusions keep the union of both partner loci; single-junction variants land here.
            depth = len({n for _, _, dn, _ in per_j for n in dn})
            span_names = set().union(*[sn for _, _, _, sn in per_j]) if per_j else set()
            counted_locus = '|'.join('%s:%d' % (c, p) for c, p, _, _ in per_j)
            c_star, p_star = chrom1, pos1
        else:
            # Two junctions of one event -> count at the better-covered end (matches IGV).
            c_star, p_star, depth_names, span_names = max(per_j, key=lambda t: len(t[3]))
            depth = len(depth_names)
            counted_locus = '%s:%d' % (c_star, p_star)
    rec['candidate_depth'] = depth
    rec['counted_locus'] = counted_locus

    # --- Numerator: the read must reproduce the variant's alignment signature.
    if vt in THROUGH_TYPES:
        supp = sig_through(case_bam, c_star, p_star, novel_right, args.min_anchor,
                           read_through_depth(seg_len, args), args.min_gap, args.window)
        sig = 'THROUGH'
    elif vt in INTERIOR_TYPES:
        if boundaries:
            supp = sig_through(case_bam, c_star, p_star, novel_right, args.min_anchor,
                               read_through_depth(seg_len, args), args.min_gap, args.window)
            sig = 'INTERIOR'
            if seg_len < args.min_anchor:
                low_conf = True   # extension too short to anchor in -> not discriminating
        else:
            # Single-point EE (pos1 == pos2): there is no extension to read through, so there is
            # nothing to measure. The old fallback -- count whichever side carries more splice-in
            # reads -- reported 66 supporting reads at loci with 21 reads present, a number wrong
            # enough to mislead anyone who used it. Report no count rather than a number we know
            # is inflated; spanning_reads and counted_locus still say where to look in IGV.
            supp = None
            sig = 'INSIDE_UNRESOLVED'
            low_conf = True
    elif vt in INSIDE_TYPES:
        if two_point:
            other = pos2 if p_star == pos1 else pos1
            supp = sig_inside(case_bam, c_star, p_star, other, args.inside_window,
                              args.min_anchor, args.tol, args.min_gap)
        else:  # single-point NE: segment side not encoded -> take the dominant splice-in side
            supp = sig_inside_dominant(case_bam, c_star, p_star, args.inside_window,
                                       args.min_anchor, args.tol, args.min_gap)
        sig = 'INSIDE'
    elif vt in GAP_TYPES:
        if chrom2 == chrom1 and pos2 is not None:
            s_exp, e_exp = expected_gap(vt, pos1, pos2)
            rivals, nearest = [], None
            if ref is not None:
                rivals = [iv for iv in ref.introns_near(chrom1, s_exp, e_exp, args.tol)
                          if iv != (s_exp, e_exp)]
                for a, b in rivals:
                    d = abs(a - s_exp) + abs(b - e_exp)
                    nearest = d if nearest is None or d < nearest else nearest
            else:
                low_conf = True   # no annotation -> cannot rule out the canonical junction
            supp = sig_gap(case_bam, chrom1, s_exp, e_exp, vt, args.tol, args.window,
                           args.gap_merge, rivals)
            if nearest is not None and nearest <= args.near_canonical:
                low_conf = True   # a canonical junction sits a few bp away: likely an artefact
                logging.info('%s (%s) at %s:%d-%d is %d bp from an annotated junction -- '
                             'shifted canonical junction, reliability=low',
                             vid, vt, chrom1, s_exp, e_exp, nearest)
                if notes is not None:
                    notes['near_canonical'] += 1
            if vt == 'DEL' and varsize and varsize < 15:
                low_conf = True   # very small deletions are confounded with ONT indel noise
        else:  # a gap type with no usable second breakpoint: fall back to a clip
            supp = sig_clip(case_bam, c_star, p_star, args.min_clip, args.bp_tol, args.window)
            low_conf = True
        sig = 'GAP'
    elif vt in SPLIT_TYPES:
        if chrom2 is not None and pos2 is not None:
            supp = sig_split(case_bam, chrom1, pos1, chrom2, pos2, args.window, args.split_tol)
        else:
            supp = sig_clip(case_bam, chrom1, pos1, args.min_clip, args.bp_tol, args.window)
        sig = 'SPLIT'
        low_conf = True   # fusions/rearrangements flagged low by policy (warrant manual review)
    elif vt in INSERTION_TYPES:
        supp = sig_insertion(case_bam, chrom1, pos1, varsize, args.ins_tol,
                             max(args.window, args.ins_tol + 100))
        sig = 'INSERTION'
        low_conf = True   # ONT insertions align unreliably in both methods
    elif vt in CLIP_TYPES:
        supp = sig_clip(case_bam, c_star, p_star, args.min_clip, args.bp_tol, args.window)
        sig = 'CLIP'
    else:
        supp = sig_clip(case_bam, c_star, p_star, args.min_clip, args.bp_tol, args.window)
        sig = 'CLIP'
        low_conf = True

    if support_writer is not None and supp:
        for nm in sorted(supp):
            support_writer.write('%s\t%s\n' % (vid, nm))

    if supp is None:
        # Nothing measurable at this locus (single-point EE). spanning_reads and counted_locus
        # still stand, so the row remains a lead to look at; the count itself is left blank
        # rather than filled with a number that cannot be justified.
        rec['supporting_read_count'] = pd.NA
        rec['spanning_reads'] = len(span_names)
        rec['junction_VAF'] = pd.NA
    else:
        span = len(span_names | supp)
        sup = len(supp)
        rec['supporting_read_count'] = sup
        rec['spanning_reads'] = span
        rec['junction_VAF'] = round(sup / span, 4) if span else 0.0
    rec['support_signature'] = sig
    rec['supporting_read_count_reliability'] = 'low' if low_conf else 'high'
    return rec


def main():
    args = parse_args()
    init_logging(args.log)

    info = pd.read_csv(args.variant_info, sep='\t', low_memory=False)
    required = ['contig_id', 'variant_id', 'variant_type', 'pos1']
    missing = [c for c in required if c not in info.columns]
    if missing:
        logging.error('variant_info missing columns: %s', missing)
        sys.exit(2)
    info = info.drop_duplicates('variant_id')
    logging.info('Read %d variants across %d contigs.', len(info), info.contig_id.nunique())

    ref = None
    if args.tx_annotation:
        ref = RefAnnotation(args.tx_annotation)
    else:
        logging.warning('No --tx_annotation: the novel part of a block cannot be located and '
                        'reads cannot be arbitrated against annotated junctions. EE and '
                        'splice/deletion counts fall back to position-only matching and are '
                        'reported reliability=low.')

    case_bam = pysam.AlignmentFile(args.case_bam, 'rb')
    support_writer = open(args.support_reads_out, 'w') if args.support_reads_out else None
    if support_writer is not None:
        support_writer.write('variant_id\tread_name\n')

    notes = {'near_canonical': 0, 'ee_no_novel_region': 0}
    rows = []
    for n, (_, row) in enumerate(info.iterrows(), 1):
        rows.append(count_variant(row, case_bam, args, ref, support_writer, notes))
        if n % 500 == 0:
            logging.info('Processed %d/%d variants.', n, len(info))

    if support_writer is not None:
        support_writer.close()

    out = pd.DataFrame(rows, columns=OUT_COLUMNS) if rows else pd.DataFrame(columns=OUT_COLUMNS)
    out.to_csv(args.out, sep='\t', index=False)

    n_lc = int((out['supporting_read_count_reliability'] == 'low').sum()) if len(out) else 0
    logging.info('%d variants sit within %d bp of an annotated junction (flagged low); '
                 '%d EE rows had no resolvable novel region (fell back to the block midpoint).',
                 notes['near_canonical'], args.near_canonical, notes['ee_no_novel_region'])
    logging.info('Wrote %d variants to %s (%d reliability=low).', len(out), args.out, n_lc)
    print('Wrote %d variants to %s (%d reliability=low).' % (len(out), args.out, n_lc))


if __name__ == '__main__':
    main()
