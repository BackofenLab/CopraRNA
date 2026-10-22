"""Build a small offline multi-genome job for full legacy/Python comparison."""

from pathlib import Path
import gzip
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

ROOT = Path(__file__).resolve().parents[1]


def build(directory):
    directory.mkdir(parents=True, exist_ok=True)
    rng = random.Random(5017)
    length = 18000
    dna = "".join(rng.choices("ACGT", k=length))
    query = "".join(rng.choices("ACGT", k=35))
    rrna = "".join(rng.choices("ACGT", k=130))
    proteins = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=80)) for _ in range(30)]
    accessions = ["NC_000913", "NC_000964", "NC_002516"]
    input = []
    for j, accession in enumerate(accessions):
        genome = list(dna)
        for pos in rng.sample(range(length), j * 130):
            genome[pos] = rng.choice("ACGT")
        record = SeqRecord(
            Seq("".join(genome)),
            id=accession + ".1",
            name=accession,
            description="Synthetic offline pipeline regression",
        )
        record.annotations = dict(
            molecule_type="DNA", topology="circular", organism=f"Test species {j}"
        )
        for i in range(30):
            tag = f"org{j}_gene{i:03d}"
            location = SimpleLocation(
                500 + i * 500, 740 + i * 500, strand=1 if i % 2 == 0 else -1
            )
            qualifiers = dict(
                locus_tag=[tag],
                gene=[f"gene{i}"],
                db_xref=[f"GeneID:{10000 + j * 100 + i}"],
                product=[f"test product {i}"],
                translation=[proteins[i]],
            )
            record.features.append(
                SeqFeature(
                    location,
                    type="gene",
                    qualifiers={
                        k: v for k, v in qualifiers.items() if k != "translation"
                    },
                )
            )
            record.features.append(
                SeqFeature(location, type="CDS", qualifiers=qualifiers)
            )
        sequence = list(record.seq)
        rrna_seq = list(rrna)
        for pos in rng.sample(range(130), j * 4):
            rrna_seq[pos] = rng.choice("ACGT")
        sequence[16500:16630] = rrna_seq
        record.seq = Seq("".join(sequence))
        record.features.append(
            SeqFeature(
                SimpleLocation(16500, 16630, strand=1),
                type="rRNA",
                qualifiers={"product": ["16S ribosomal RNA"]},
            )
        )
        with gzip.open(directory / (accession + ".gb.gz"), "wt") as stream:
            SeqIO.write(record, stream, "genbank")
        q = list(query)
        if j:
            q[j * 3] = rng.choice("ACGT")
        input.append(f">{accession}\n" + "".join(q) + "\n")
    (directory / "input_sRNA.fa").write_text("".join(input))
    mapping = (ROOT / "coprarna_aux/kegg2refseqnew.csv").read_text().splitlines()
    codes = [
        next(line.split("\t")[0] for line in mapping if accession in line)
        for accession in accessions
    ]
    table = ["#cluster\t" + "\t".join(codes) + "\n"]
    for i in range(30):
        table.append(
            str(i + 1)
            + "\t"
            + "\t".join(f"{code}:org{j}_gene{i:03d}(1)" for j, code in enumerate(codes))
            + "\n"
        )
    (directory / "cluster.tab").write_text("".join(table))


if __name__ == "__main__":
    import sys

    build(
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else ROOT / "tests/fixtures/pipeline/input"
    )
