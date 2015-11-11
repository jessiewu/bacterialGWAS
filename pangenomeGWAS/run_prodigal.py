#Jane Charlesworth
#jane.charlesworth@ndm.ox.ac.uk
#15/05/15
#
#code to construct a pan-genome set from contigs or genomes.
#standalone code to run Prodigal (gene finder) on set of contig files

import argparse
import os
from subprocess import call
from Bio import SeqIO
from Bio.Blast import NCBIXML

#get and check options
opts = argparse.ArgumentParser(description="Construct pan-genome from contigs, using CD-hit")
opts.add_argument('input_seqs', action="store",nargs='?', help="List of contig or genome files to make a pan-genome from, in fasta or gzipped fasta format.")

options = opts.parse_args()

#check sequences exist
assert os.path.exists(options.input_seqs), "list of sequences not found"
seq_list = [line.rstrip() for line in open(options.input_seqs)]
for _file in seq_list:
        assert os.path.isfile(_file), 'sequence file does not exist: %r' % _file

#define functions
def prodigal(isolate):
	"""annotate open reading frames on a set of contigs, using Prodigal"""
	isolate_name= isolate.split(".")[0]
        out = "%s_prodigal.faa" % isolate
	call(["/dipro/mmm/gorm/v2/ana/Saur/jchar/pangenome/scripts/Prodigal-GoogleImport/prodigal",\
		 "-i", isolate, "-a", out, "-c", "-m", "-o", "prodigal.out", "-q"])


#call genes with prodigal
for _file in seq_list:
	if _file.endswith(".gz"):
		if _file.rstrip(_file.split("/")[-1]) != os.getcwd():
			call(["cp", _file, "./"])
			sep ="/"
			_file = sep.join((os.getcwd(),_file.split("/")[-1]))
			call(["gzip","-f", "-d", _file])
		else:
			call(["gzip","-f", "-d", _file])
		_file = _file.rstrip(".gz")
	assert len(list(SeqIO.parse(_file, 'fasta'))) > 0, 'Sequence file is not in fasta format: %r' % _file
	prodigal(_file)
	call(["rm", _file])

#clean up
call(["rm", "prodigal.out"])

