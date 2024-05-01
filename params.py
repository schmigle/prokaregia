import argparse
import sys
import torch
import pkg_resources

from tiara.constants import (
    DESCRIPTION,
    MIN_SEQ_LEN,
    DEFAULT_PROB_CUTOFF,
    FIRST_NNET_KMER2PARAMS,
    SECOND_NNET_KMER2PARAMS,
)

from tiara.src.classification import Classification


def get_classifier(args: argparse.Namespace) -> Classification:
    torch.set_num_threads(args.threads)

    first_stage_model_params = FIRST_NNET_KMER2PARAMS[args.first_stage_kmer]
    second_stage_model_params = SECOND_NNET_KMER2PARAMS[args.second_stage_kmer]

    if len(args.prob_cutoff) == 1:
        prob_cutoff_1, prob_cutoff_2 = args.prob_cutoff[0], args.prob_cutoff[0]
    else:
        prob_cutoff_1, prob_cutoff_2 = args.prob_cutoff

    first_stage_model_params["prob_cutoff"] = prob_cutoff_1
    second_stage_model_params["prob_cutoff"] = prob_cutoff_2
    first_stage_model_params["fragment_len"] = 5000
    second_stage_model_params["fragment_len"] = 5000
    first_stage_model_params["dim_out"] = 5
    second_stage_model_params["dim_out"] = 3

    first_stage_tfidf_fname = f"k{first_stage_model_params['k']}-first-stage"
    second_stage_tfidf_fname = f"k{second_stage_model_params['k']}-second-stage"

    params = [first_stage_model_params, second_stage_model_params]
    nnet_weights = [
        pkg_resources.resource_filename(__name__, "models/nnet-models/" + path)
        for path in [
            first_stage_model_params["fname"],
            second_stage_model_params["fname"],
        ]
    ]
    tfidfs = [
        pkg_resources.resource_filename(__name__, "models/tfidf-models/" + path)
        for path in [first_stage_tfidf_fname, second_stage_tfidf_fname]
    ]

    return Classification(
        min_len=args.min_len,
        nnet_weights=nnet_weights,
        params=params,
        tfidf=tfidfs,
        threads=args.threads,
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-i",
        "--input",
        metavar="input",
        help="A path to a fastq/a file.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output",
        help="A path to output file. If not provided, the result is printed to stdout.",
        default=None,
    )
    parser.add_argument(
        "-m",
        "--min_len",
        help=f"""Minimum length of a sequence. Sequences shorter than min_len are discarded. 
        Default: {MIN_SEQ_LEN}.""",
        type=int,
        default=MIN_SEQ_LEN,
    )
    parser.add_argument(
        "--first_stage_kmer",
        "--k1",
        help=f"k-mer length used in the first stage of classification. Default: 6.",
        type=int,
        default=6,
    )
    parser.add_argument(
        "--second_stage_kmer",
        "--k2",
        help=f"k-mer length used in the second stage of classification. Default: 7.",
        type=int,
        default=7,
    )
    parser.add_argument(
        "-p",
        "--prob_cutoff",
        metavar="cutoff",
        help=f"""Probability threshold needed for classification to a class. 
        If two floats are provided, the first is used in a first stage, the second in the second stage
        Default: {DEFAULT_PROB_CUTOFF}.""",
        nargs="+",
        type=float,
        default=DEFAULT_PROB_CUTOFF,
    )

    parser.add_argument(
        "--to_fastq/a",
        "--tf",
        metavar="class",
        help="""Write sequences to fastq/a files specified in the arguments to this option.
        The arguments are: mit - mitochondria, pla - plastid, bac - bacteria, 
        arc - archaea, euk - eukarya, unk - unknown, pro - prokarya, 
        all - all classes present in input fasta (to separate fasta files).""",
        nargs="+",
        type=str,
        default=[],
    )

    parser.add_argument(
        "-t", "--threads", help="Number of threads used.", type=int, default=1
    )

    parser.add_argument(
        "--probabilities",
        "--pr",
        action="store_true",
        help="""Whether to write probabilities of individual classes for each sequence to the output.""",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="""Whether to display some additional messages and progress bar during classification.""",
    )

    parser.add_argument(
        "--gzip",
        "--gz",
        action="store_true",
        help="""Whether to gzip results or not.""",
    )
    parser.add_argument(
        "--fq", 
        "--fastq", 
        action="store_true", 
        help = "input and output are in fastq format"
    )
    parser.add_argument(
        "--fa", 
        "--fasta", 
        action="store_true", 
        help = "input and output are in fasta format"
    )

    return parser.parse_args(args=None if sys.argv[1:] else ["--help"])