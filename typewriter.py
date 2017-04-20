#!/usr/bin/env python

# ============================================================================ #
# Copyright (c) 2015 University of Oxford                                      #
# Tanya Golubchik                                                              #
# golubchi@well.ox.ac.uk                                                       #
# March 2015 (mk_typewriter: April 2012)                                       #
# ============================================================================ #

"""
TYPEWRITER

  Searches a given list of blast databases or fasta files for either:
    - a panel of one or more sequences (eg resistance genes) from a multifasta query file, or
    - the MLST profile for each sample
  and produces combined output for the dataset. Similarity is reported as
  relative coverage, calculated as (sequence identity * query coverage).

  With --save_results also outputs the hits for each query sequence in fasta format.

  With --do_mutations reports sequence differences (eg. drug resistance mutations)
  in the form XNY, where X is the reference (subject) residue, N is the position in
  the reference sequence, and Y is the query residue at this position.

  Example: MLST on assembled data from a single sample:
  {typewriter.py} --infastapath /path/to/my/sample/contigs.fasta
                  --outprefix my_mlst_results
                  --outdir    /path/to/outut/directory
                  --mlst True
                  --word_size 11
                  --allelefile  /path/to/my/allelefile.tsv
                  --stfile      /path/to/my/STfile.tsv
"""
# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

from __future__ import division

import sys, os, time
import subprocess as sp
import argparse
import cStringIO as StringIO
import cPickle as pickle

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Blast.Applications import NcbitblastnCommandline as tblastn
import Bio.Application

import itertools
from collections import namedtuple, deque, Counter
import uuid, shutil
import numpy as np
import re

# ============================================================================ #
# GLOBAL VARIABLES                                                             #
# ============================================================================ #
_args = None
_tmpdir = None
_blast_paths = None
_alleles = None
_STs = None
_outfile_prefix = None
_progname=os.path.basename(sys.argv[0])

# ============================================================================ #
# LOGGING                                                                      #
# ============================================================================ #
def loginfo(s):
    sys.stdout.write('  {0}\n'.format(s))
def logerr(s):
    sys.stderr.write('  Warning: {0}\n'.format(s))
def stoperr(s, errcode=1):
    errword = 'Finished' if not errcode else 'Error'
    sys.stderr.write('  {0}: {1}\n'.format(errword, s))
    sys.exit(errcode)

# ============================================================================ #
# PROGRAM USAGE                                                                #
# ============================================================================ #

def Initialise():
    """
    Parse command-line arguments.
    """
    global _args

    def str2bool(s):
        """Return the Python boolean value True or False dependingon what the string is."""
        return (s is not None and s.lower() in ['1', 't', 'true', 'y', 'yes'])

    parser = argparse.ArgumentParser(description=__doc__, epilog='Example: {progname} --blastdb my_contigs --outdir my_mlst_results/ --outprefix CdifMLST --mlst True --allelefile my-allele-file.tsv --stfile --my-st-file.tsv --orgid Cdif --save_results True'.format(progname=_progname))

    # paths for input and output
    parser.add_argument( '--outdir', default=None, help='Output directory for results', required=True)
    parser.add_argument( '--outprefix', default=None, help='Prefix for naming output files (Required)', required=True)
    parser.add_argument( '--querypath', default=None, help='Fasta file with one or more query sequences. Not required for MLST (default: None)')
    parser.add_argument( '--blastdb', help='Path to a single blast database (prefix before the .nhr etc extensions)' )
    parser.add_argument( '--blastdbfile', help='File with list of blast databases of velvet assemblies to search' )
    parser.add_argument( '--inbampath', help='Path to a single BAM file to search against. Temporary blast database will be created from this.' )
    parser.add_argument( '--inbampathfile', help='File with list of BAM files to search against; temporary blast databases will be made for these. Note files must have unique names.' )
    # optional program-specific
    parser.add_argument( '--infastqpath', help='Path to a single FastQ file to search against, eg sample_contigs.fa. Temporary blast database will be created from this.' )
    parser.add_argument( '--infastqpathfile', help='File with list of FastQ files to search against; temporary blast databases will be made for these. Note files must have unique names.' )
    parser.add_argument( '--infastapath', help='Path to a single Fasta file to search against, eg sample_contigs.fa. Temporary blast database will be created from this.' )
    parser.add_argument( '--infastapathfile', help='File with list of Fasta files to search against; temporary blast databases will be made for these. Note files must have unique names.' )
    # optional program-specific
    parser.add_argument( '--blastprog', choices=('blastn', 'megablast', 'tblastn'), default='blastn', help='Blast algorithm. Use blastn (more sensitive) or megablast (faster) for nucleotide queries, tblastn for protein queries (default: blastn)')
    parser.add_argument( '--min_cov', type=int, metavar='INT', default=90, help='For tblastn only: report mutations for matches with at least this percent coverage (default: 90)' )
    parser.add_argument( '--min_contig_kcov', type=int, metavar='INT', default=5, help='Ignore potentially contaminating low-coverage contigs (default: 5)' )
    parser.add_argument( '--min_contig_length', type=int, metavar='INT', default=300, help='Ignore short possibly spurious contigs (default: 300)' )
    parser.add_argument( '--word_size', type=int, metavar='INT', default=17, help='Size of anchor (word) for blastn searches (default: 17)' )
    parser.add_argument( '--ignore_poor_contigs', type=str2bool, metavar='BOOL', default=False, help='Ignore low-confidence contigs (kcov < min_contig_kcov or length < min_contig_length) (default: False)')
    parser.add_argument( '--consensus_overlaps', type=str2bool, metavar='BOOL', default=False, help='If there are hits that overlap with mismatches, take the consensus for the mismatched sites. Default behaviour is to trust the sequence of the best hit. Only use this if working with raw reads, where many high-quality overlapping hits are expected. (default: False)')
    parser.add_argument( '--do_mutations', type=str2bool, metavar='BOOL', default=False, help='Report specific mutations. (Default: True for tblastn (protein), False for blastn (nucleotide))' )
    parser.add_argument( '--hitonly', type=str2bool, metavar='BOOL', default=False, help='Do not report results for queries where there were no hits for any sample (default: False)' )
    parser.add_argument( '--save_results', type=str2bool, metavar='BOOL', default=False, help='Save all blast results in fasta format (default: False)' )
    parser.add_argument( '--save_split_cov', type=str2bool, metavar='BOOL', default=False, help='Save coverage and identity output separately, in addition to reporting these as "relative coverage" (default: False)' )
    parser.add_argument( '--save_markerstatus', type=str2bool, metavar='BOOL', default=False, help='Save tab file for loading to ogre MARKERSTATUS table (default: False unless running in MLST mode)' )
    parser.add_argument( '--save_qual', type=str2bool, metavar='BOOL', default=False, help='Save csv file of contig qualities (default: False). Assumes contig length and coverage is present in fasta headers of contigs file.' )
    parser.add_argument( '--genepospath', default=None, help='Path to file containing names of genes/queries of interest and the specific positions within these to examine for assembly quality. Format is a fasta-like file with ">gene_name" for each gene of interest, followed by a line of comma-separated positions, starting from 1 or mutations in the format XNY. Must match the type and length of the query sequences (ie amino acid or nucleotide). (default: None)')
    parser.add_argument( '--mlst', type=str2bool, metavar='BOOL', default=False, help='Treat search as in-silico MLST (default: False)' )
    parser.add_argument( '--allelefile', default=None, help='Path to tab-separated file containing alleles for MLST, each line in the form (genecd alleleid dnaseq).')
    parser.add_argument( '--stfile', default=None, help='Path to tab-separated file containing MLST profiles, each line in the form (ST a1 a2 a3 a4 a5 a6 a7).')
    parser.add_argument( '--orgid', default=None, help='Organism ID; must be specified for MLST (default: None).')

    parser.add_argument( '--refid', default=None, help='Reference ID; must be specified for MLST (default: None).')
    parser.add_argument( '--outmlst', default='outmlst', help='output file for MLST result with the sequence type')
    parser.add_argument( '--outmarkerstatus', default='outmarkerstatus', help='output file for MLST results with detailed information of presence of gene')
    parser.add_argument( '--outmarkerdenovo', default='outmarkerdenovo', help='output file for MLST results with information needed on de novo sequence type')
    # paths to third-party software
    parser.add_argument( '--makeblastdbpath', help='Path to makeblastdb program.', default='makeblastdb')
    parser.add_argument( '--samtoolspath', help='Path to samtools.', default='samtools')

    # general optional
    parser.add_argument( '--dbuserid', default=None, help='User ID for creating database-compatible output (default: unknown).')
    parser.add_argument( '--verbose', type=str2bool, metavar='BOOL', default=True, help='Produce detailed logging output (default: True)')
    parser.add_argument( '--overwrite', type=str2bool, metavar='BOOL', default=False, help='Overwrite existing output files if using the same prefix (default: False)')

    _args = parser.parse_args()

    loginfo('Started {0}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

    return


def _ungap(strm):
    """
    Remove gap characters (hyphens) from query file - required for makeblastdb.
    """
    for s in strm:
        s.seq = s.seq.ungap('-')
        yield s

def Bam2Fasta(inbampath):
    """
    Create query file with the first max_queries reads from the downsampled BAM file.
    """
    global _tmpdir, _args
    # extract only reads that have not failed QC and are not PCR duplicates
    cmd1 = '{samtools} view -F 1536 {inbampath}'.format(samtools=_args.samtoolspath, inbampath=inbampath)
    cmd2 = 'cut -f1,10'
    loginfo('Extracting reads from {0}.'.format(inbampath))
    loginfo(' | '.join((cmd1, cmd2)))
    p1 = sp.Popen(cmd1.split(), stdout=sp.PIPE, stderr=sp.PIPE)
    p2 = sp.Popen(cmd2.split(), stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    p1.stdout.close() # mimic pipe
    stdout, stderr = p2.communicate()
    if not stdout:
        stoperr('Failed to get reads from {b}: {e}'.format(b = inbampath, e = stderr if stderr else 'None'))

    try:
        outfastapath = os.path.join(_tmpdir, '{0}_queries.fasta'.format(os.path.splitext(os.path.basename(inbampath))[0]))
    except:
        stoperr('Failed to create query file : {e}'.format(e = stderr if stderr else 'None'))
    with open(outfastapath, 'w') as o:
        o.writelines('>{0}'.format(line.replace('\t','\n')) for line in StringIO.StringIO(stdout) if line.split('\t'))

    return outfastapath


def Fastq2Fasta(infastqpath):
    """
    Dump reads from fastq input file.
    """
    global _tmpdir, _args
    cmd = ['awk', 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print};if(P==4)P=0;P++}', infastqpath]
    loginfo('Extracting reads from {0}.'.format(infastqpath))
    loginfo(' '.join(cmd))
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if not stdout:
        stoperr('Failed to get reads from {b}: {e}'.format(b = infastqpath, e = stderr if stderr else 'None'))
    try:
        outfastapath = os.path.join(_tmpdir, '{0}_queries.fasta'.format(os.path.splitext(os.path.basename(infastqpath))[0]))
    except:
        stoperr('Failed to create query file : {e}'.format(e = stderr if stderr else 'None'))
    with open(outfastapath, 'w') as o:
        o.write(stdout)
    return outfastapath


def Get_Paths():
    ''' Get dictionary of blastdbfile, indexed by path, and MLST info as needed.'''
    blast_paths = {}
    global _alleles, _STs, _args

    # check all other blastdbs and add paths
    blastdblist = deque()
    if _args.blastdb:
        blastdblist.append(_args.blastdb)

    if _args.blastdbfile:
        with open(_args.blastdbfile) as inf:
            blastdblist.extend(line.split()[0] for line in inf if line.split())

    # create blast databases if we only got raw files: first save fasta files
    fastapathlist = deque()
    if _args.inbampath:
        fastapathlist.append(Bam2Fasta(_args.inbampath))
    if _args.infastqpath:
        fastapathlist.append(Fastq2Fasta(_args.infastqpath))
    if _args.infastapath:
        fastapathlist.append(_args.infastapath)

    if _args.inbampathfile:
        with open(_args.inbampathfile) as inf:
            fastapathlist.extend(Bam2Fasta(line.split()[0]) for line in inf if line.split())
    if _args.infastqpathfile:
        with open(_args.infastqpathfile) as inf:
            fastapathlist.extend(Fastq2Fasta(line.split()[0]) for line in inf if line.split())
    if _args.infastapathfile:
        with open(_args.infastapathfile) as inf:
            fastapathlist.extend(line.split()[0] for line in inf if line.split())

    # now create the blast databases
    if fastapathlist:
        for fastapath in fastapathlist:
            blastdb = os.path.join(_tmpdir, os.path.basename(fastapath))
            if os.path.exists(blastdb):
                os.remove(blastdb)
            try:
                SeqIO.write(_ungap(SeqIO.parse(fastapath, 'fasta')), blastdb, 'fasta')
            except:
                stoperr('Could not create blastdb {0}'.format(blastdb))
            print _tmpdir, blastdb
            cmd = '{makeblastdb} -dbtype nucl -in {fastapath} -out {blastdb}'.format(makeblastdb=_args.makeblastdbpath, fastapath=fastapath, blastdb=blastdb)
            loginfo('Making blast database: {cmd}'.format(cmd=cmd))
            p = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()

            if stderr != '' and os.path.isfile('{blastdb}.nhr'.format(blastdb=blastdb)):
                logerr(stderr)
                stoperr('Failed to create blast database from {0}.'.format(fastapath))
            else:
                blastdblist.append(blastdb)
                if _args.verbose:
                    loginfo('Added {blastdb}.'.format(blastdb=blastdb))

    # save paths to unique blastdbs
    for blastdb in blastdblist:
        if os.path.isfile(blastdb + '.nhr'):
            if not blast_paths.has_key(blastdb) and blastdb not in blast_paths.values():
                blast_paths[blastdb] = blastdb
            else:
                loginfo('Database {db} already in list.'.format(db=blastdb))
        else:
            stoperr('File does not exist : {db}.nhr'.format(db=blastdb))

    if not blast_paths:
        stoperr('No blast databases to add. Nothing to do, exiting.', errcode=0)

    return blast_paths


def Clean_Commandline():
    """
    Print errors, and exit if necessary, on bad input data.
    """
    global _tmpdir, _outfile_prefix, _blast_paths, _alleles, _STs

    if not os.path.exists(_args.outdir):
        try:
            os.mkdir(_args.outdir)
        except:
            stoperr('Unable to create output directory {0}. Exiting.'.format(_args.outdir))

    # set temporary directory as one level below outdir
    _tmpdir = os.path.join(_args.outdir, 'tmp_typewriter')
    if not os.path.isdir(_tmpdir):
        os.mkdir(_tmpdir)

    # get input paths and set up output prefix
    if (not _args.querypath and not _args.mlst) or (_args.querypath and not os.path.isfile(_args.querypath)):
        stoperr('Query FASTA file {0} does not exist.'.format(_args.querypath))
    if not os.path.isdir(_args.outdir):
        stoperr('Output directory {0} does not exist or is inaccessible'.format(_args.outdir))
    if _args.outprefix:
        from string import ascii_letters, digits
        if not frozenset(_args.outprefix).issubset(frozenset(ascii_letters + digits + '-_')):
            stoperr('--outprefix must contain only letters/digits/underscore/dash characters.')
    if _args.querypath:
        try:
            seq_ix = SeqIO.index(_args.querypath, 'fasta')
            assert seq_ix
        except ValueError:
            stoperr('Query file contains duplicate fasta headers.')
        loginfo('File contains {0} queries.'.format(len(seq_ix)))
        # Ungap query file for running blast if original has gaps
        if sum(1 for x in seq_ix if '-' in seq_ix[x].seq):
            loginfo('Query file contains gaps - ungapping.')
            ungapped_path = os.path.join(_tmpdir, os.path.basename(_args.querypath))
            nq = SeqIO.write(_ungap(SeqIO.parse(_args.querypath, 'fasta')), ungapped_path, 'fasta')
            if nq != len(seq_ix):
                stoperr('Number of ungapped queries ({0}) does not match number of queries ({1}) in original file.'.format(nq, len(seq_ix)))
            _args.querypath = ungapped_path
    if _args.genepospath and (not os.path.isfile(_args.genepospath) or not os.path.getsize(_args.genepospath)):
        stoperr('--genepospath file is not accessible.')

    # set generic outprefix
    if _args.mlst and not _args.consensus_overlaps:
        suffix = 'mlst'
    elif _args.mlst and _args.consensus_overlaps:
        suffix = 'rawmlst'
    else:
        suffix = os.path.splitext(os.path.basename(_args.querypath))[0]
    _outfile_prefix = os.path.join(_args.outdir, '{prefix}-{suffix}-'.format(prefix=_args.outprefix, suffix=suffix))
    if _args.verbose:
        loginfo('Output filepaths will start with "{0}."'.format(_outfile_prefix))

    # blast flavour
    if not _args.blastprog and not _args.mlst:
        stoperr('No blast flavour specified (blastn, megablast or tblastn).')
    elif _args.mlst:
        _args.blastprog = 'blastn'

    # decide whether to look at mutations depending on blast flavour
    if _args.do_mutations or (_args.blastprog == 'tblastn' and not _args.do_mutations):
        _args.do_mutations = True
    else:
        _args.do_mutations = False

    # check word size makes sense
    if _args.word_size != 17:
        if _args.blastprog == 'tblastn':
            stoperr('word_size has no effect for tblastn.')
        try:
            assert int(_args.word_size) > 0
        except:
            stoperr('word_size must be a positive integer.')

    # set minimum query coverage required to report mutations
    if _args.min_cov:
        try:
            _args.min_cov = float(_args.min_cov)
            assert 0. <_args.min_cov <= 100.
        except:
            stoperr('min_cov must be between 1 and 100.')

    if _args.min_contig_kcov != 5:
        try:
            _args.min_contig_kcov = int(_args.min_contig_kcov)
            assert _args.min_contig_kcov > 0
        except:
            stoperr('min_contig_kcov must be a positive integer (default: 5).')

    if _args.min_contig_length != 300:
        try:
            _args.min_contig_kcov = int(_args.min_contig_kcov)
            assert _args.min_contig_kcov > 0
        except:
            stoperr('min_contig_length must be a positive integer (default: 300).')

    # refuse incompatible argument combinations
    if _args.hitonly and _args.mlst:
        stoperr('Must report all query results in MLST mode.')
    if _args.mlst and (not (_args.allelefile) or (not _args.stfile) or not os.path.isfile(_args.allelefile) or not os.path.isfile(_args.stfile)):
        stoperr('Paths to existing allele file and ST file are required for MLST.')

    # confirm we have either a list of fasta files or databases to search
    if not _args.blastdb and not _args.blastdbfile and not _args.inbampath and not _args.inbampathfile and not _args.infastqpath and not _args.infastqpathfile and not _args.infastapath and not _args.infastapathfile:
        stoperr('No blast databases or BAM/FastQ/Fasta input files supplied.')
    if _args.blastdbfile and not os.path.isfile(_args.blastdbfile) or (_args.inbampathfile and not os.path.isfile(_args.inbampathfile)):
        stoperr('File containing input paths is not accessible.')

    # try opening blast databases if these are given, also fetch MLST alleles if needed
    _blast_paths = Get_Paths()
    if _args.verbose:
        loginfo('Got {n} blast databases.'.format(n=len(_blast_paths)))

    # set up MLST queries
    if _args.mlst:
        if _args.allelefile and os.path.isfile(_args.allelefile) and _args.stfile and os.path.isfile(_args.stfile):
            _alleles = np.loadtxt(_args.allelefile, dtype=np.str, delimiter='\t')
            _STs = dict((tuple(line[1:]), line[0]) for line in np.loadtxt(_args.stfile, dtype=np.str, delimiter='\t'))
             # check that we have MLST loci for the organism we want
            allele_orgid = set(allelecd.split('-')[0] for allelecd in _alleles[:,1])
            if _args.orgid not in allele_orgid:
                stoperr('Allele file does not contain any data for {orgid}'.format(_args.orgid))
            else:
                # remove any rows in the input file that have data for other organisms
                @np.vectorize
                def correct_orgid(a): return a.startswith(_args.orgid)
                _alleles = _alleles[correct_orgid(_alleles).sum(axis=1)>0]
        elif _alleles==None or _STs==None or not len(_alleles) or not len(_STs):
            stoperr('Paths to existing allele file and ST file are required for MLST.')

        if _args.verbose:
            loginfo('Retrieved {n} allele entries.'.format(n=len(_alleles)))

        # save sequence of allele 1 for each locus, to use as queries for the blast search
        # TODO: Change this to consensus of each locus
        loci = np.unique(_alleles[:,0])
        first_alleles = (_to_record(*_alleles[_alleles[:,0]==locus][0]) for locus in loci)
        _args.querypath = os.path.join(_tmpdir, _args.outprefix + '-mlst_queries.fasta')
        SeqIO.write(first_alleles, _args.querypath, 'fasta')

    final_outfiles = ['relcov.csv', 'numhits.csv', 'contig_names.csv']
    if _args.do_mutations:
        final_outfiles.extend(['mutations.csv', 'mutfreq.csv', 'mutlist.txt'])
    if _args.mlst:
        final_outfiles.extend(['MARKER_DENOVO.tab', 'MLST.csv'])
    if _args.save_markerstatus:
        final_outfiles.append('MARKERSTATUS.tab')
    if _args.save_qual:
        final_outfiles.append('contig_quality.csv')
    if _args.save_split_cov:
        final_outfiles.extend(['cov.csv', 'ident.csv'])
    # Check if final output already exists
    Check_Final_Output(final_outfiles)
    return

def Check_Final_Output(outfiles):
    """
    Check if final outpaths exist in outdir; if so, don't do anything else, just exit.
    """
    if _args.overwrite:
        return

    all_present = True
    for outfile in outfiles:
        if not os.path.isfile(_outfile_prefix + outfile):
            all_present = False

    if all_present:
        if _args.verbose:
            stoperr('Nothing to do; all final output files exist in {outdir} : {outpaths}'.format(outdir=_args.outdir, outpaths=', '.join(outfiles)), 0)
    return


# ============================================================================ #
# DATA PROCESSING FUNCTIONS                                                    #
# ============================================================================ #

def Get_Gene_Sites():
    """
    Parse file in _args.genepospath and return a dictionary of genes, each with a sorted list of unique positions of interest.
    """
    with open(_args.genepospath) as g_fp:
        sites_dic = dict(l.strip().split('\n') for l in  g_fp.read().split('>') if l.strip())
    mut_pat = re.compile(r'[A-Za-z]+(\d+)')
    def _get_mut_pos(mut_string):
        # format was mutations in standard XNY notation
        l = frozenset(int(m) for m in re.findall(mut_pat, mut_string))
        if not l:
            # format was numeric amino acid positions
            l = frozenset(int(m.strip()) for m in mut_string.split(',') if m.strip().isdigit())
        return sorted(l)
    for gene_name in sites_dic:
        sites_dic[gene_name] = _get_mut_pos(sites_dic[gene_name])
    if _args.verbose:
        loginfo('Sites of interest: {0}'.format(sites_dic))
    return sites_dic


def _to_record(name, desc, seq_string):
    """
    Format sequence for writing with Biopython.
    """
    return SeqRecord(id=name, name=name, description=desc, seq=Seq(seq_string))


def _drop_ins_in_range(x, y, inslist, pat=re.compile(r'ins(\d+)')):
    """
    Select insertions outside given range [x,y).
    """
    for ins in inslist:
        match = re.match(pat, ins)
        if  match and not (x <= int(match.group(1)) < y):
            yield ins
        else:
            loginfo('dropping insertion {0}.'.format(ins))

def Run_Blast(min_contig_kcov=5, min_contig_length=300, min_cov=90):
    """
    Blast nucleotide queries against given databases.
    """
    save_split_cov = True if _args.save_split_cov else False
    task = _args.blastprog # prefer use blastn instead of megablast; better sensitivity
    blastprog=eval(task) if task in ('blastn', 'tblastn') else eval('blastn')
    min_hit_length = 20 if blastprog == tblastn else 50

    # create dictionary of query numbers and sequences indexed by query name
    markergene = namedtuple('markergene', ('pos', 'seq'))
    if _args.verbose:
         loginfo('Queries : {0}.'.format(_args.querypath))
    panel = dict((s.name, markergene(i, s)) for i, s in enumerate(_ungap(SeqIO.parse(_args.querypath, 'fasta'))))
    rel_cov_array = np.zeros((len(_blast_paths), len(panel)), dtype=np.float)
    num_hits_array = np.zeros((len(_blast_paths), len(panel)), dtype=np.int)
    if save_split_cov:
        cov_array = np.zeros((len(_blast_paths), len(panel)), dtype=np.float)
        ident_array = np.zeros((len(_blast_paths), len(panel)), dtype=np.float)
    contig_id_dic = dict.fromkeys(panel)
    for k in panel:
        contig_id_dic[k] = dict.fromkeys(_blast_paths)
    if _args.save_qual:
        contig_qual_array = np.zeros((len(_blast_paths), len(panel), 2), dtype=np.float)
    if _args.mlst:
        mlst_array = np.zeros((len(_blast_paths), len(panel)), dtype='S20')
    if _args.do_mutations:
        mutation_dic = dict.fromkeys(panel)
        for k in panel:
            mutation_dic[k] = dict.fromkeys(_blast_paths)
    if _args.genepospath:
        # get dictionary of unique positions to examine for quality, by gene
        sites_dic = Get_Gene_Sites()
        # dictionary in which to save the results, per blastdb
        gene_pos_dic = dict.fromkeys(panel)
        for k in panel:
            # check sites are valid
            if sites_dic.has_key(k) and len(sites_dic[k]) and not (1<=max(sites_dic[k])<=len(panel[k].seq)):
                stoperr('Site {0} is invalid for query {1}.'.format(max(sites_dic[k]), k))
            gene_pos_dic[k] = dict.fromkeys(_blast_paths)

    # list of matching sequences for queries, with the query as first sequence in each
    # result records is a full list of all sequence matches with headers showing the source
    # conflated records is the best guess at a single best match for each genome
    result_records = dict.fromkeys(panel)
    conflated_records = dict.fromkeys(panel)
    for g in panel.iterkeys():
        result_records[g] = [panel[g].seq]
        conflated_records[g] = [panel[g].seq]

    # iterate over databases, one set of results for each db
    db_names = sorted(_blast_paths.keys())

    # set blast parameters to match NBCI recommendations
    gapopen = 5 if blastprog == blastn else 11
    gapextend = 2 if blastprog == blastn else 1
    word_size = _args.word_size if blastprog == blastn else 3

    # scaffolding/masking marker
    scaff = set('X') if blastprog == tblastn else set('N')

    blastn_cline = blastprog(query=_args.querypath, word_size=word_size, gapopen=gapopen, gapextend=gapextend, culling_limit=1, evalue=0.001, outfmt=5)
    if task != 'tblastn':
        blastn_cline.task = task
    loginfo(blastn_cline)
    for i, db_name in enumerate(db_names):
        loginfo('Processing {0}.'.format(db_name))
        blastn_cline.db=_blast_paths[db_name]
        try:
            blast_xml, blast_err = blastn_cline()
        except Bio.Application.ApplicationError, err:
            logerr(err)
            loginfo('Trying -culling_limit 2.')
            blastn_cline.culling_limit=2
            try:
                blast_xml, blast_err = blastn_cline()
            except:
                stoperr('Could not run blast on {0}.'.format(db_name))
                continue

        if blast_err:
            logerr('Blast Error:' + blast_err)
            continue

        # each blast_rec is result from one query sequence; iterate over these
        for query_blast_rec in NCBIXML.parse(StringIO.StringIO(blast_xml)):
            noncontig_warn_flag = False
            lowconf_warn_flag = False
            q_name = query_blast_rec.query.split()[0]
            if not panel.has_key(q_name):
                logerr('query {q} not found in query file.'.format(q=q_name))
                continue
            q_len = query_blast_rec.query_length

            # relative coverage across query sequence
            q_cov = np.zeros(q_len, dtype=np.float)
            if save_split_cov:
                q_rawcov = np.zeros(q_len, dtype=np.float)
                q_rawident = np.zeros(q_len, dtype=np.float)

            # best guess at a matching sequence for this query, with hit parts amalgalmated
            fillchar = 'X' if blastprog == tblastn else '-'
            q_seq = np.repeat(fillchar, q_len)
            unknownchar = 'X' if blastprog == tblastn else 'N'
            # per-position hit quality, each row is [contig_coverage, contig_length]
            q_qual = np.zeros((q_len, 2))
            # contig IDs of good hits for this query
            q_contig_ids = []
            # dictionary of sites where mismatched overlaps have been found
            if _args.consensus_overlaps:
                mismatch_dic = {}

            # total number of matches, counting all hsp parts individually
            q_numhits = 0
            new_seq_from, new_seq_to = q_len, 0
            for aln_num, alignment in enumerate(query_blast_rec.alignments):
                contig_desc = alignment.hit_def
                contig_length = alignment.length
                try:
                    contig_cov = float(contig_desc.split('cov_')[-1])
                except:
                    contig_cov = None
                # ignore low-confidence contigs
                if (contig_cov and contig_cov < min_contig_kcov) or contig_length < min_contig_length:
                    if not lowconf_warn_flag:
                        logerr('Match to {q} in low-confidence contig: {desc}.'.format(q=q_name, desc=contig_desc))
                        lowconf_warn_flag = True
                    if _args.ignore_poor_contigs:
                        loginfo('ignoring this match.'.format(q=q_name))
                        continue

                aln_seq = np.repeat('-', q_len)

                # calculate relative query coverage for this match, across all parts of the match
                for hsp_num, hsp in enumerate(alignment.hsps):
                    frompos, topos = hsp.query_start, hsp.query_end
                    perc_id = 100.*(float(hsp.identities) / len(hsp.query))
                    hit_length = topos - frompos + 1
                    if hit_length < min_hit_length:
                        continue

                    # if we get this far, we have found a good hit
                    rel_cov = (perc_id * hsp.align_length)/hsp.align_length
                    curr_cov = np.mean(q_cov[frompos-1:topos])
                    # proceed only if hit is better than what we have, or there are gaps in the current match
                    if rel_cov > curr_cov or (q_seq[frompos-1:topos]==fillchar).any():
                        q_numhits += 1
                        q_cov[frompos-1:topos] = rel_cov
                        if save_split_cov:
                            q_rawcov[frompos-1:topos] = 100.*(float(hsp.align_length)/hit_length)
                            q_rawident[frompos-1:topos] = perc_id
                        # forget insertions previously found in this query range
                        if _args.do_mutations and mutation_dic[q_name][db_name]:
                            mutation_dic[q_name][db_name] = list(_drop_ins_in_range(frompos,topos,mutation_dic[q_name][db_name]))
                        # select positions where query doesn't have a gap (excludes insertions in subject sequence)
                        if hsp.align_length == hit_length:
                            # alignment without gaps, save as is
                            aln_seq[frompos-1:topos] = np.asarray(list(hsp.sbjct), dtype='S1')
                        elif hsp.align_length - hsp.sbjct.count('-') == hit_length:
                            # aligned with gaps, but length matches query length, so can just ungap this
                            # this almost always happens because of incorrect gap placement by blastn
                            aln_seq[frompos-1:topos] = np.asarray(list(hsp.sbjct.replace('-','')))
                        else:
                            # gappy alignment, need to cut down to size
                            sel = np.where(np.asarray(list(hsp.query), dtype='S1') != '-')[0]
                            aln_seq[frompos-1:topos] = np.asarray(list(hsp.sbjct), dtype='S1')[sel]
                            if '-' in hsp.query:
                                # insertions in subject - may cause frameshifts etc
                                pat = re.compile(r'-+')
                                for insertion in re.finditer(pat, hsp.query):
                                    ins_value = hsp.sbjct[insertion.start():insertion.end()]
                                    ins_pos = frompos-1 + len(hsp.query[:insertion.start()].replace('-',''))
                                    if _args.verbose:
                                        loginfo('Insertion of {0} residues {1} at position {2} in {3} match to {4}.'.format((insertion.end()-insertion.start()), ins_value, ins_pos, db_name, q_name))
                                    # look only at mutations above coverage/identity threshold and not consisting entirely of scaffolding character
                                    if _args.do_mutations and rel_cov > min_cov and set(ins_value) != scaff:
                                        if not mutation_dic[q_name][db_name]: mutation_dic[q_name][db_name] = []
                                        mutation_dic[q_name][db_name].append('ins{pos}{ins}'.format(pos=ins_pos, ins=ins_value))

                        # check if we're overwriting something, and if we are, report this in the mutations file
                        if not(q_seq[frompos-1:topos] == fillchar).all():
                            if not noncontig_warn_flag:
                                loginfo('Match to {0} is composed of multiple non-contiguous parts.'.format(q_name))
                                noncontig_warn_flag = True
                            if _args.do_mutations and (not mutation_dic[q_name][db_name] or 'non-contiguous' not in mutation_dic[q_name][db_name]):
                                if not mutation_dic[q_name][db_name]: mutation_dic[q_name][db_name] = []
                                mutation_dic[q_name][db_name].append('non-contiguous')
                            # identify mismatches in overlaps
                            mismatch_positions = np.where((q_seq[frompos-1:topos] != aln_seq[frompos-1:topos]) & (q_seq[frompos-1:topos] != fillchar) & (aln_seq[frompos-1:topos] != fillchar))[0]+frompos-1
                            if mismatch_positions.any():
                                logerr('Mismatches found in overlapping hits, treat mutations with caution.')
                                # save the bases seen at mismatching positions, for resolving consensus later
                                if _args.consensus_overlaps:
                                    for mpos in mismatch_positions:
                                        if not mismatch_dic.has_key(mpos):
                                            mismatch_dic[mpos] = deque(q_seq[mpos] * (q_numhits-1))
                                        mismatch_dic[mpos].append(aln_seq[mpos])
                            # record values at any previously observed mismatched positions
                            if _args.consensus_overlaps:
                                for prev_mpos in mismatch_dic:
                                    if prev_mpos not in mismatch_positions:
                                        mismatch_dic[prev_mpos].append(q_seq[mpos])

                        # add to conflated sequences
                        q_seq[frompos-1:topos] = aln_seq[frompos-1:topos]
                        # add matching quality info
                        q_qual[frompos-1:topos] = (contig_cov, contig_length)

                        # save names of contigs that make up this hit
                        q_contig_ids.append('{0} ({1}..{2})'.format(contig_desc, hsp.sbjct_start, hsp.sbjct_end))

                        # all matching subject sequences aligned at starting position
                        if _args.save_results:
                            new_seq_tag = '{0}_{1}-{2}'.format(db_name, frompos, topos)
                            new_seq_desc = '{0} length={1} cov={2}'.format(contig_desc, contig_length, contig_cov)
                            new_seq = '-'*(frompos-1) + hsp.sbjct
                            result_records[q_name].append(_to_record(new_seq_tag, new_seq_desc, new_seq))

            # generate conflated query sequence string from array of chars
            q_seq = q_seq.tostring()
            # add in mismatches if wanting a consensus
            if _args.consensus_overlaps:
                for mpos, mbases in mismatch_dic.iteritems():
                    basecounts = Counter(mbases)
                    loginfo('Mismatches at position {0}: {1}'.format(mpos, basecounts))
                    toptwo = basecounts.most_common()[:2]
                    if toptwo[0][1] == toptwo[1][1]:
                        loginfo('Mixture is 50/50, replacing with {0}'.format(unknownchar))
                        most_common = unknownchar
                    else:
                        most_common = toptwo[0][0]
                    logerr('Replacing mismatch with consensus {0}'.format(most_common))
                    q_seq = q_seq[:mpos] + most_common + q_seq[mpos+1:]


            # calculate relative coverage over query range
            relcov = sum(q_cov)/q_len if q_cov!=None else 0.
            if save_split_cov:
                rawcov = sum(q_rawcov)/q_len if q_rawcov!=None else 0.
                # identity should be undiluted by remaining uncovered sequence, so
                # sum over covered sequence only
                covered_len = (q_rawident != 0).sum()
                rawident = sum(q_rawident)/covered_len if (covered_len and q_rawident!=None) else 0.

            # store relative coverage, number of hits, contig quality and contig name(s) for this query
            rel_cov_array[i, panel[q_name].pos] = relcov
            num_hits_array[i, panel[q_name].pos] = q_numhits
            contig_id_dic[q_name][db_name] = q_contig_ids
            if save_split_cov:
                cov_array[i, panel[q_name].pos] = rawcov
                ident_array[i, panel[q_name].pos] = rawident
            if _args.save_qual:
                contig_qual_array[i, panel[q_name].pos] = q_qual.mean(axis=0) # mean contig_cov, mean contig_length

            # if required, get quality info at specific positions of interest in this sequence
            if _args.genepospath:
                if sites_dic.has_key(q_name):
                    for spos in sites_dic[q_name]:
                        if not gene_pos_dic[q_name][db_name]: gene_pos_dic[q_name][db_name] = []
                        gene_pos_dic[q_name][db_name].append('{ccov:.2f};{clen:.2f}'.format(ccov=q_qual[spos-1,0], clen=q_qual[spos-1,1]))

            # save concatenated sequences for this query
            if _args.save_results:
                conflated_name = '{0}_{1}-{2}'.format(db_name, new_seq_from, new_seq_to)
                conflated_desc = 'hits={0} relcov={1}'.format(q_numhits, relcov)
                conflated_records[q_name].append(_to_record(conflated_name, conflated_desc, q_seq))

            # create list of substitutions in this (conflated) sequence
            if _args.do_mutations:
                if not mutation_dic[q_name][db_name]: mutation_dic[q_name][db_name] = []
                sbj_seq = panel[q_name].seq.seq.tostring()
                if q_seq == sbj_seq and not mutation_dic[q_name][db_name]:
                    mutation_dic[q_name][db_name] = ['wt']
                else:
                    startpos = q_seq.rfind(q_seq.lstrip(fillchar))
                    endpos = len(q_seq.rstrip(fillchar))
                    if startpos == len(q_seq):
                        mutation_dic[q_name][db_name] = ['absent']
                    elif relcov < min_cov:
                        mutation_dic[q_name][db_name] = ['identity <{0:.0f}'.format(min_cov)]
                    else:
                        if startpos != 0:
                            mutation_dic[q_name][db_name].append("truncated 5' {0}aa".format(startpos))
                        for aapos in xrange(startpos, endpos):
                            if q_seq[aapos] != sbj_seq[aapos] and sbj_seq[aapos] != fillchar:# and q_seq[aapos] != 'X':
                                # mutation format is "XNY (cov;len)" where X = wild-type residue, N = position from 1, Y = substituted residue, cov=contig coverage, len=contig length
                                mutation_dic[q_name][db_name].append('{subj}{pos}{query}({ccov:.2f};{clen:.2f})'.format(subj=sbj_seq[aapos], pos=aapos+1, query=q_seq[aapos], ccov=q_qual[aapos,0], clen=q_qual[aapos,1]))
                        if endpos != len(q_seq):
                            mutation_dic[q_name][db_name].append("truncated 3' {0}aa".format(len(q_seq)-endpos))
                        if not mutation_dic[q_name][db_name]:
                            # sequence matches and no mutations, so must differ only at masked characters
                            mutation_dic[q_name][db_name] = ['wt']

            # if mlst, save allele number
            if _args.mlst:
                pat = re.compile(r'({seq})'.format(seq=q_seq.upper().replace('N', '.').replace('-', '.')))
                locus_alleles = _alleles[_alleles[:,0]==q_name]
                possible_alleles = [alpos for alpos, a in enumerate(re.match(pat, s) for s in locus_alleles[:,2]) if a]
                if len(possible_alleles) > 1:
                    loginfo('Sample {db} has {n} potential {locus} matches: {matches}.'.format(db=db_name, n=len(possible_alleles), locus=q_name, matches=locus_alleles[possible_alleles, 1]))
                    # if not found, create allele name like 'Saur-arcC-Notfound'
                    allele_id = locus_alleles[0,1].rpartition('-')[0] + '-Notfound'
                elif possible_alleles:
                    # save allele match
                    allele_id = locus_alleles[possible_alleles[0],1]
                else:
                    logerr('Could not determine MLST allele {q} for sample {db}.\n\tSequence: {seq}'.format(q=q_name, db=db_name, seq=q_seq))
                    # if not found, create allele name like 'Saur-arcC-Notfound'
                    allele_id = locus_alleles[0,1].rpartition('-')[0] + '-Notfound'
                mlst_array[i, panel[q_name].pos] = allele_id

    # sort columns alphabetically if doing MLST, otherwise in order given in the query file
    queries = sorted(panel.keys(), key=lambda x: panel[x].pos) if not _args.mlst else sorted(panel.keys())
    query_order = [panel[q].pos for q in queries]
    loginfo('Queries were {0}'.format(', '.join(queries)))
    loginfo('Query order {0}'.format(query_order))

    # drop columns where there were no hits
    if _args.hitonly:
        queries = np.asarray(queries)
        query_order = np.asarray(query_order)
        # require all columns to have at least one sample with relcov>1%
        hitcols = (rel_cov_array.max(axis=0)>1.).astype(np.bool)
        loginfo('No hits for {n} queries.'.format(n=(len(queries)-hitcols.sum())))
        if _args.verbose:
            loginfo('Queries without hits in any sample : {0}.'.format(', '.join(queries[:,~hitcols])))
        queries = queries[:,hitcols].tolist()
        query_order = query_order[:,hitcols].tolist()

    if _args.verbose:
        loginfo('Final output has {n} queries.'.format(n=len(queries)))

    # dump arrays to file with row and column labels
    header_line = ['BLASTDB'] + queries
    row_labels = np.asarray(db_names).reshape((len(db_names),1))
    np.savetxt(_outfile_prefix+'relcov.csv', np.vstack((header_line, np.hstack((row_labels, rel_cov_array[:,query_order].astype('S20'))))), fmt='%s', delimiter=',')
    np.savetxt(_outfile_prefix+'numhits.csv', np.vstack((header_line, np.hstack((row_labels, num_hits_array[:,query_order].astype('S20'))))), fmt='%s', delimiter=',')
    # save split raw coverage and raw percent identity if required
    if save_split_cov:
        np.savetxt(_outfile_prefix+'cov.csv', np.vstack((header_line, np.hstack((row_labels, cov_array[:,query_order].astype('S20'))))), fmt='%s', delimiter=',')
        np.savetxt(_outfile_prefix+'ident.csv', np.vstack((header_line, np.hstack((row_labels, ident_array[:,query_order].astype('S20'))))), fmt='%s', delimiter=',')
    # save contig quality
    if _args.save_qual:
        with open(_outfile_prefix+'contig_quality.csv', 'w') as contig_qual_fp:
            contig_qual_fp.write(','.join(header_line) + '\n')
            for row_label, row in zip(db_names, contig_qual_array[:,query_order]):
                contig_qual_fp.write('{0},{1}\n'.format(row_label, ','.join('{0:.2f};{1:.2f}'.format(ccov, clen) for (ccov, clen) in row)))
    # save contig identifiers for good hits
    with open(_outfile_prefix+'contig_names.csv', 'w') as o:
        o.write(','.join(header_line) + '\n')
        for db_name in db_names:
            o.write(db_name)
            for q_name in queries:
                contig_id_list = contig_id_dic[q_name][db_name]
                if contig_id_list:
                    o.write(',' + ';'.join(contig_id_list))
                else:
                    o.write(',')
            o.write('\n')
    if _args.mlst:
        with open(_outfile_prefix+'MLST.csv', 'w') as o:
            o.writelines('{row_label},{stid},{alleles}\n'.format(
                row_label=row_labels[i][0],
                stid=(_STs[tuple(row)] if _STs.has_key(tuple(row)) else 'NF'),
                alleles=','.join(tuple(row))) for i, row in enumerate(mlst_array[:,query_order]))
        os.rename(_outfile_prefix+'MLST.csv', _args.outmlst)

    # save quality at specific sites
    if _args.genepospath:
        with open(_outfile_prefix+'posquality.csv', 'w') as o:
            o.write('BLASTDB,')
            for q_name in queries:
                if sites_dic.has_key(q_name):
                    for site in sites_dic[q_name]:
                        o.write('{0} {1},'.format(q_name, site))
            o.write('\n')
            for db_name in db_names:
                o.write(db_name)
                for q_name in queries:
                    slist = gene_pos_dic[q_name][db_name]
                    if slist:
                        o.write(',' + ','.join(slist))
                o.write('\n')

    # save mutations if relevant
    if _args.do_mutations:
        with open(_outfile_prefix+'mutations.csv', 'w') as o:
            o.write(','.join(header_line) + '\n')
            for db_name in db_names:
                o.write(db_name)
                for q_name in queries:
                    mutlist = mutation_dic[q_name][db_name]
                    if mutlist:
                        o.write(',' + ' '.join(mutlist))
                    else:
                        o.write(',absent')
                o.write('\n')
        # mutation frequencies
        with open(_outfile_prefix+'mutfreq.csv', 'w') as o:
            with open(_outfile_prefix+'mutlist.txt', 'w') as mutlist_fp:
                o.write('#QUERY,SAMPLES_AFFECTED,MUTATION,COUNT\n')
                for q_name in mutation_dic:
                    all_mutations = sorted(m.split('(')[0] for m in itertools.chain(*(mutation_dic[q_name].itervalues())))
                    mutlist_fp.write('>{0}\n{1}\n'.format(q_name, ','.join(set(all_mutations))))
                    # frequencies of all unique mutations, in descending order, excluding truncations and the markers 'low cov' & 'absent'
                    mut_counts = sorted(((m, sum(1 for _ in g)) for m, g in itertools.groupby(all_mutations) if m != 'low cov' and m != 'absent'), key=lambda x: -x[-1])
                    num_samples = sum(1 for mutlist in mutation_dic[q_name].itervalues() if len(mutlist) >= 1)
                    o.writelines('{0},{1},{2},{3}\n'.format(q_name, num_samples, m, count) for m,count in mut_counts)



    if _args.save_results:
        for q in queries:
            if len(result_records[q]) > 1:
                SeqIO.write(result_records[q], '{prefix}{gene}.fasta'.format(prefix=_outfile_prefix, gene=q), 'fasta')
            if len(conflated_records[q]) > 1:
                SeqIO.write(conflated_records[q], '{prefix}{gene}-conflated.fasta'.format(prefix=_outfile_prefix, gene=q), 'fasta')

    mlst_return = mlst_array[:,query_order] if _args.mlst else None

    return queries, rel_cov_array[:,query_order], num_hits_array[:,query_order], mlst_return


# ============================================================================ #
# OUTPUT FOR DB TABLES                                                         #
# ============================================================================ #
def Write_Tabs(queries, rel_cov_array, num_hits_array, mlst_array):
    """
    Save DB-compatible output files for loading into database.
    """
    dbbirthdt = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

    if _args.mlst:
        mlst_outfile = open(_outfile_prefix + 'MARKER_DENOVO.tab', 'w')
        mlst_outfile.write(('\t'.join(('blastdb', 'genecd', 'alleleid', 'status', 'comment', 'dbbirthdt')) + '\n'))

    if _args.save_markerstatus or _args.mlst:
        # output format suitable for loading into ogre table MARKERSTATUS
        with open(_outfile_prefix + 'MARKERSTATUS.tab',  'w') as o:
            o.write(('\t'.join(('blastdb', 'genename', 'status', 'blastnrelcov', 'blastnhits', 'comment', 'dbbirthdt')) + '\n'))
            for i, db_name in enumerate(sorted(_blast_paths.keys())):
                for j, genename in enumerate(queries):
                    if num_hits_array[i,j] <= 0 or rel_cov_array[i,j] <= 50:
                        status = 'Absent'
                    elif 50 < rel_cov_array[i,j] <= 98:
                        status = 'Partial'
                    elif rel_cov_array[i,j] > 98:
                        status = 'Present'
                    o.write('{guid}\t{genename}\t{status}\t{relcov}\t{hits}\t\t{dbbirthdt}\n'.format(
                        guid=db_name,
                        genename=genename,
                        status=status,
                        relcov=rel_cov_array[i,j],
                        hits=num_hits_array[i,j],
                        dbbirthdt=dbbirthdt))

                    # output format suitable for loading into ogre table MARKER_DENOVO
                    if _args.mlst:
                        if mlst_array[i,j].endswith('Notfound'):
                            allele_status = 'NotFound'
                        elif num_hits_array[i, j] == 1:
                            allele_status = 'Exact'
                        else:
                            allele_status = 'Partial'
                        mlst_outfile.write('{guid}\t{genecd}\t{alleleid}\t{status}\t\t{dbbirthdt}\n'.format(
                            guid=db_name,
                            genecd=genename,
                            alleleid=mlst_array[i,j],
                            status=allele_status,
                            dbbirthdt=dbbirthdt))
    return


# ============================================================================ #
# MAIN                                                                         #
# ============================================================================ #

if __name__ == '__main__':
    Initialise()
    Clean_Commandline()
    queries, rel_cov_array, num_hits_array, mlst_array = Run_Blast(min_contig_kcov = _args.min_contig_kcov,
                                                               min_contig_length = _args.min_contig_length,
                                                               min_cov=_args.min_cov)
    Write_Tabs(queries, rel_cov_array, num_hits_array, mlst_array)
    loginfo('Success; {tmpdir} will now be deleted if it is empty.'.format(tmpdir=_tmpdir))
    try:
        os.rmdir(_tmpdir)
    except:
        pass

# ============================================================================ #
