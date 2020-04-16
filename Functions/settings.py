"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains various hard-coded Unicycler settings.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

# When aligning minimap reads to the graph (which is overlap-free), we still want to allow a tiny
# bit of overlap because minimap alignments are a bit course.
ALLOWED_MINIMAP_OVERLAP = 5

# When aligning minimap reads to the graph, we can exclude hits that are too much worse than the
# best hit.
MAX_TO_MIN_MINIMISER_RATIO = 10

# When testing various repeat counts using fully global alignment in Seqan, we use this band size
# to make the alignment faster.
SIMPLE_REPEAT_BRIDGING_BAND_SIZE = 50

# Illumina contigs are used as 'reads' in the miniasm and Racon steps. They are given this as a
# qscore at each base.
CONTIG_READ_QSCORE = 40

# This is the maximum number of times an assembly will be Racon polished
RACON_POLISH_LOOP_COUNT_HYBRID = 2
RACON_POLISH_LOOP_COUNT_LONG_ONLY = 4


# This is the number of times assembly graph contigs are included in the Racon polish reads. E.g.
# if 6, then each contig is included 6 times as a read (3 forward strand 3 reverse).
RACON_CONTIG_DUPLICATION_COUNT = 1

CONTIG_SEARCH_END_SIZES = [5000, 2500, 1000, 500]
CONTIG_SEARCH_MIN_IDENTITY = 95.0
FOUND_CONTIG_MIN_RATIO = 0.9
FOUND_CONTIG_MAX_RATIO = 1.11111
FOUND_CONTIG_MAX_OVERLAP_SIZE = 250

# Unicycler will only work with read alignments if they are long enough. This values specifies
# the minimum threshold.
MIN_LONG_READ_ALIGNMENT_LENGTH = 50

# If less than this fraction of a read was aligned in the first aligning pass, Unicycler will try
# again using much more sensitive alignment settings. This helps to align reads which come from
# particularly difficult repetitive regions.
MIN_READ_FRACTION_ALIGNED = 0.9

# This is how much overlap is allowed between two alignments in a single read, relative to the
# graph's overlap. For example, if the graph has an overlap of 95 and this value is 1.1,
# then alignments within a read can go up to 105 bp, but alignments with more overlap will be
# filtered out.
ALLOWED_ALIGNMENT_OVERLAP = 1.1

# Unicycler will not use the lowest quality alignments for making bridges. This setting specifies
# the threshold. E.g. if it is 5, then any alignment with a scaled score of less than the 5th
# percentile scaled score will be thrown out.
MIN_SCALED_SCORE_PERCENTILE = 5.0

# Unicycler-align can automatically determine a low score threshold. It does this by randomly
# aligning 100 bp sequences with the current scoring scheme and determining the mean and standard
# deviation of such random alignments. The threshold is then set to a certain number of standard
# deviations above the mean (this setting). This should ensure that any alignment which passes
# the threshold is at least a little bit better than a random sequence alignment.
AUTO_SCORE_STDEV_ABOVE_RANDOM_ALIGNMENT_MEAN = 7

# When Unicycler is searching for paths connecting two graph segments which matches a read
# consensus sequence, it will only consider paths which have a length similar to the expected
# sequence (based on the consensus sequence length). These settings define the acceptable range.
# E.g. if they are 0.7 and 1.3, Unicycler will consider graph paths that range from 70% to 130%
# of the expected sequence length.
MIN_RELATIVE_PATH_LENGTH = 0.9
MAX_RELATIVE_PATH_LENGTH = 1.1
RELATIVE_PATH_LENGTH_BUFFER_SIZE = 100

# These settings are used when Unicycler is exhaustively searching for paths connecting two graph
# segments. If the number of working paths or the number of final paths gets too high during the
# search (exceeds these thresholds), Unicycler will give up and instead try a progressive path
# search.
ALL_PATH_SEARCH_MAX_WORKING_PATHS = 10000
ALL_PATH_SEARCH_MAX_FINAL_PATHS = 500

# These settings are used when Unicycler is progressively searching for paths connecting two graph
# segments. When its number of working paths reaches PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS, it
# will cull them down by scoring the alignment of each. Paths which have a score within the
# PROGRESSIVE_PATH_SEARCH_SCORE_FRACTION of the best are kept.
PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS = 100
PROGRESSIVE_PATH_SEARCH_SCORE_FRACTION = 0.995

# These settings are used for Unicycler's copy number determination - the process by which it
# tries to figure out the depth of constituent components of each segment.
#   * INITIAL_SINGLE_COPY_TOLERANCE controls how much excess depth is acceptable for the first
#     single copy assignment pass.
#   * COPY_PROPAGATION_TOLERANCE controls how much discrepancy is allowed when propagating copy
#     depths from one segment to the next.
#   * MIN_SINGLE_COPY_LENGTH is how short of a segment can be called single copy when adding
#     additional single copy segments.
#   * MAX_COPY_DEPTH_DISTRIBUTION_ARRANGEMENTS caps the number of possible ways to redistribute a
#     segment's copy depths to its neighbours. If there are more possibilities than this,
#     Unicycler won't bother trying.
INITIAL_SINGLE_COPY_TOLERANCE = 0.1
COPY_PROPAGATION_TOLERANCE = 0.5
MIN_SINGLE_COPY_LENGTH = 1000
MAX_COPY_DEPTH_DISTRIBUTION_ARRANGEMENTS = 10000
COPY_DEPTH_PROPAGATION_TABLE_ROW_WIDTH = 35

# When Unicycler is cleaning up the graph after bridging, it can delete graph paths and graph
# components which are mostly (but not entirely) used up in bridges. This value controls the
# threshold. Making it smaller will make Unicycler more willing to delete stuff. Making it larger
# will make Unicycler less willing to delete stuff.
CLEANING_USEDUPNESS_THRESHOLD = 0.5

# When making a consensus sequence, Unicycler will use up to this many read sequences. It is
# capped because a consensus with too many reads will take too long and likely not be much
# better. E.g. a 100 read consensus will likely give a similar sequence as a 25 read consensus,
# but it will take much longer.
MAX_READS_FOR_CONSENSUS = 25

# The different bridging modes have different minimum bridge quality thresholds.
CONSERVATIVE_MIN_BRIDGE_QUAL = 25.0
NORMAL_MIN_BRIDGE_QUAL = 10.0
BOLD_MIN_BRIDGE_QUAL = 1.0

# These control how often progress lines are updated in the stdout. They define the percentage step
# used. E.g. if set to 5.0, the progress will go 5%, 10%, 15%, etc.
# Setting these to higher values helps to prevent excessive progress updates, which is a pain when
# piping Unicycler output to file.
LOADING_REFERENCES_PROGRESS_STEP = 1.0
LOADING_READS_PROGRESS_STEP = 1.0
LOADING_ALIGNMENTS_PROGRESS_STEP = 1.0

# These settings control how willing Unicycler is to make bridges that don't have a graph path.
# This depends on whether one or both of the segments being bridged ends in a dead end and
# whether we have any expected linear sequences (i.e. whether real dead ends are expected).
# If we don't expect any linear sequences, then bridging two dead ends with a pathless bridge is
# great (not at all penalised), bridging a dead end with a non-dead end is okay, and bridging two
# non-dead ends is very much discouraged (quite penalised). If we do expect linear sequences, then
# we are less willing to make bridges between dead ends because it's possible those dead ends are
# genuine and should remain dead ends.
PATHLESS_BRIDGE_QUAL_TWO_DEAD_ENDS = 1.0
PATHLESS_BRIDGE_QUAL_ONE_DEAD_END = 0.7
PATHLESS_BRIDGE_QUAL_NO_DEAD_ENDS = 0.2
PATHLESS_BRIDGE_QUAL_TWO_DEAD_ENDS_WITH_LINEAR_SEQS = 0.6
PATHLESS_BRIDGE_QUAL_ONE_DEAD_END_WITH_LINEAR_SEQS = 0.4
PATHLESS_BRIDGE_QUAL_NO_DEAD_ENDS_WITH_LINEAR_SEQS = 0.2

# If the user doesn't set the thread manually, it will use either the CPU count or this value,
# whichever is smaller. This is to prevent Unicycler from grabbing too many cores by default.
# E.g. if it was run on a large machine with 80 cores, it shouldn't use all 80 unless the user
# explicitly asks for it!
MAX_AUTO_THREAD_COUNT = 8

# The default sequence line wrapping length (e.g. for use in FASTA files).
BASES_PER_FASTA_LINE = 70

# Pilon is run multiple times to polish things up as nice as possible. It will stop when no more
# changes are made or this limit is hit.
MAX_PILON_POLISH_COUNT = 10

MINIASM_BRIDGE_QUAL_WITH_GRAPH_PATH = 1.0
MINIASM_BRIDGE_QUAL_WITH_DEAD_END = 1.0
MINIASM_BRIDGE_QUAL_WITHOUT_PATH_OR_DEAD_END = 0.7
MINIASM_BRIDGE_SCALED_SCORE_TO_USE_GRAPH_PATH = 95.0
MINIASM_BRIDGE_HALF_QUAL_LENGTH = 5000

LONG_READ_BRIDGE_HALF_QUAL_LENGTH = 2000


# If the miniasm assembly is too small, we won't even consider using for bridging in hybrid
# assembly. This size is relative to the estimated genome size from the short read assembly.
REQUIRED_MINIASM_ASSEMBLY_SIZE_FOR_BRIDGING = 0.5

# Unicycler trims short read contigs which end in dead ends based on miniasm trimming. This setting
# limits the amount of trimming it's willing to do. I.e. if miniasm trimmed more than this from a
# contig, Unicycler won't.
MAX_MINIASM_DEAD_END_TRIM_SIZE = 100
