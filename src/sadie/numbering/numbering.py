from functools import lru_cache
import logging
import warnings
import sys

import pandas as pd
from numpy import nan

from sadie.antibody.exception import LongHCDR3Error, AnarciDecreasing
from .germlines import all_germlines
from .schemes import (
    number_kabat_heavy,
    number_kabat_light,
    number_chothia_heavy,
    number_chothia_light,
    number_martin_heavy,
    number_martin_light,
    number_imgt,
    number_aho,
    number_wolfguy_heavy,
    number_wolfguy_light,
)

logger = logging.getLogger("NUMBERING")


class Numbering:
    def __init__(self):
        pass

    @lru_cache(maxsize=None)
    def get_identity(self, state_sequence, germline_sequence):
        """
        Get the partially matched sequence identity between two aligned sequences.
        Partial in the sense that gaps can be in the state_sequence.
        """
        # Ensure that the sequences are the expected length
        assert len(state_sequence) == len(germline_sequence) == 128
        n, m = 0, 0
        for i in range(128):
            if germline_sequence[i] == "-":
                continue
            if state_sequence[i].upper() == germline_sequence[i]:
                m += 1
            n += 1

        if not n:
            return 0
        return float(m) / n

    @staticmethod
    def validate_numbering(xxx_todo_changeme, name_seq=[]):
        """
        Wrapper to do some basic validation of the numbering.

        Further validation could be done but at the moment we just check that the numbering indices are incremental (they should be)
        """
        (numbering, start, end) = xxx_todo_changeme
        name, seq = name_seq
        last = -1
        nseq = ""

        for (index, _), a in numbering:
            if index < last:
                raise AnarciDecreasing(name, "decreasing sequence count in the numbering")

                # , "Numbering was found to decrease along the sequence %s. Please report." % name
            last = index
            nseq += a.replace("-", "")

        assert nseq in seq.replace("-", ""), (
            "The algorithm did not number a contiguous segment for sequence %s. Please report" % name
        )

        return numbering, start, end

    @staticmethod
    def number_sequence_from_alignment(state_vector, sequence, scheme="imgt", chain_type=None):
        """
        Given you have an alignment. Give back the numbering

        @param state_vector: List of states from the hmm. Effectively these are imgt columns but CDR3 has not been redone.
        @param sequence: The original sequence string or list.
        @param scheme: The numbering scheme to apply
        @param chain_type: The type of chain to apply numbering for. Some schemes do not require this (IMGT). Others (e.g. Chothia/Wolfguy) do.

        @return: A list of numbering identifier / amino acids tuples over the domain that has been numbered. The indices of the start (inclusive) and end point (exclusive) in the sequence for the numbering
        """
        scheme = scheme.lower()
        if scheme == "imgt":
            return number_imgt(state_vector, sequence)
        elif scheme == "chothia":
            if chain_type == "H":
                return number_chothia_heavy(state_vector, sequence)
            elif chain_type in "KL":
                return number_chothia_light(state_vector, sequence)
            else:
                raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        elif scheme == "kabat":
            if chain_type == "H":
                return number_kabat_heavy(state_vector, sequence)
            elif chain_type in "KL":
                return number_kabat_light(state_vector, sequence)
            else:
                raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        elif scheme == "martin":
            if chain_type == "H":
                return number_martin_heavy(state_vector, sequence)
            elif chain_type in "KL":
                return number_martin_light(state_vector, sequence)
            else:
                raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        elif scheme == "aho":
            return number_aho(
                state_vector, sequence, chain_type
            )  # requires the chain type to heuristically put the CDR1 gap in position.
        elif scheme == "wolfguy":
            if chain_type == "H":
                return number_wolfguy_heavy(state_vector, sequence)
            elif chain_type in "KL":
                return number_wolfguy_light(state_vector, sequence)
            else:
                raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))

    def run_germline_assignment(self, state_vector, sequence, chain_type, allowed_species=None):
        """
        Find the closest sequence identity match.
        """
        genes = {"v_gene": [None, None], "j_gene": [None, None]}

        # Extract the positions that correspond to match (germline) states.
        state_dict = dict(((i, "m"), None) for i in range(1, 129))
        state_dict.update(dict(state_vector))
        state_sequence = "".join(
            [sequence[state_dict[(i, "m")]] if state_dict[(i, "m")] is not None else "-" for i in range(1, 129)]
        )

        # Iterate over the v-germline sequences of the chain type of interest.
        # The maximum sequence identity is used to assign the germline
        if chain_type in all_germlines["V"]:
            if allowed_species is not None:
                _allowed = []
                for sp in allowed_species:
                    if sp not in all_germlines["V"][chain_type]:
                        logger.debug(f"removeing {sp} from all types since it does not exists for {chain_type}")
                        continue
                    else:
                        _allowed.append(sp)
            else:
                allowed_species = _allowed
            seq_ids = {}
            for species in allowed_species:
                if species not in all_germlines["V"][chain_type]:
                    continue  # Previously bug.
                for gene, germline_sequence in all_germlines["V"][chain_type][species].items():
                    seq_ids[(species, gene)] = self.get_identity(state_sequence, germline_sequence)
            genes["v_gene"][0] = max(seq_ids, key=lambda x: seq_ids[x])
            genes["v_gene"][1] = seq_ids[genes["v_gene"][0]]

            # Use the assigned species for the v-gene for the j-gene.
            # This assumption may affect exotically engineered abs but in general is fair.
            species = genes["v_gene"][0][0]
            if chain_type in all_germlines["J"]:
                if species in all_germlines["J"][chain_type]:
                    seq_ids = {}
                    for gene, germline_sequence in all_germlines["J"][chain_type][species].items():
                        seq_ids[(species, gene)] = self.get_identity(state_sequence, germline_sequence)
                    genes["j_gene"][0] = max(seq_ids, key=lambda x: seq_ids[x])
                    genes["j_gene"][1] = seq_ids[genes["j_gene"][0]]

        return genes

    def number_sequences_from_alignment(
        self,
        sequences,
        alignments,
        scheme="imgt",
        allow=set(["H", "K", "L", "A", "B", "G", "D"]),
        assign_germline=False,
        allowed_species=None,
    ):
        """
        Given a list of sequences and a corresponding list of alignments from run_hmmer apply a numbering scheme.
        """

        # Iteration over the sequence alignments performing the desired numbering
        numbered = []
        alignment_details = []
        hit_tables = []

        for i in range(len(sequences)):

            # Unpack
            hit_table, state_vectors, detailss = alignments[
                i
            ]  # We may have multiple domains per sequence (e.g. single chain fvs).

            # Iterate over all the domains in the sequence that have been recognised (typcially only 1 with the current hmms available)
            hit_numbered, hit_details = [], []
            for di in range(len(state_vectors)):
                state_vector = state_vectors[di]
                details = detailss[di]
                details["scheme"] = scheme
                details["query_name"] = sequences[i][0]

                # Only number things that are allowed. We still keep the alignment details and hit_table
                if state_vector and details["chain_type"] in allow:
                    try:
                        # Do the numbering and validate (for development purposes)
                        try:
                            hit_numbered.append(
                                self.validate_numbering(
                                    self.number_sequence_from_alignment(
                                        state_vector,
                                        sequences[i][1],
                                        scheme=scheme,
                                        chain_type=details["chain_type"],
                                    ),
                                    sequences[i],
                                )
                            )
                        except LongHCDR3Error as e:
                            e.sequence_name = details["query_name"]
                            raise e
                        except AnarciDecreasing:
                            warnings.warn(f"Skipping {details['query_name']}", UserWarning)
                            continue
                            # raise e

                        if assign_germline:
                            details["germlines"] = self.run_germline_assignment(
                                state_vector,
                                sequences[i][1],
                                details["chain_type"],
                                allowed_species=allowed_species,
                            )
                        hit_details.append(details)
                    except AssertionError as e:  # Handle errors. Those I have implemented should be assertion.
                        print(str(e), file=sys.stderr)
                        raise e  # Validation went wrong. Error message will go to stderr. Want this to be fatal during development.
                    except Exception as e:
                        print(
                            "Error: Something really went wrong that has not been handled",
                            file=sys.stderr,
                        )
                        print(str(e), file=sys.stderr)
                        raise e

            if hit_numbered:
                numbered.append(hit_numbered)
                alignment_details.append(hit_details)
            else:
                numbered.append(None)
                alignment_details.append(None)
            hit_tables.append(hit_table)

        return numbered, alignment_details, hit_tables

    def parsed_output(self, sequences, numbered, details):
        """
        Write numbered sequences to dataframe
        Kappa and Lambda chains are written to the same file
        The sequences will written aligned to the numbering scheme. Gaps in the sequences with respect to the alignment are written
        as a '-'
        @param sequences: List of name, sequence tuples
        @param numbered: Numbered sequences in the same order as the sequences list.
        @param details: List of alignment details in the same order as the sequences list.
        """

        chain_types = {}
        pos_ranks = {}
        all_pos = {}
        _lc = {"K": "KL", "L": "KL"}

        # Divide the set into chain types and find how to order the numbering for each type.
        for i in range(len(sequences)):  # Iterate over entries
            if numbered[i] is None:
                continue

            for j in range(len(numbered[i])):  # Iterate over domains.
                # Record the chain type index
                c = details[i][j]["chain_type"]
                c = _lc.get(c, c)  # Consider lambda and kappa together.
                chain_types.setdefault(c, []).append((i, j))
                if c not in pos_ranks:
                    pos_ranks[c] = {}
                    all_pos[c] = set()

                # Update the insertion order for the scheme. i.e. is it A B C or C B A (e.g. imgt 111 and 112 repectively)
                tmp_p = -1
                r = 0
                for p, _ in numbered[i][j][0]:
                    if p[0] != tmp_p:
                        tmp_p = p[0]
                        r = 0
                    else:
                        r += 1
                    pos_ranks[c][p] = max(r, pos_ranks[c].get(p, r))
                    all_pos[c].add(p)

        summary_dataframes = []
        # Write a new file for each chain type. Kappa and lambda are written together as light chains.
        for cts in ["H", "KL", "A", "B", "G", "D"]:
            if cts in chain_types:

                # Sort the positions by index and insertion order
                positions = sorted(all_pos[cts], key=lambda p: (p[0], pos_ranks[cts][p]))

                # Header line
                fields = [
                    "Id",
                    "sequence",
                    "domain_no",
                    "hmm_species",
                    "chain_type",
                    "e-value",
                    "score",
                    "seqstart_index",
                    "seqend_index",
                    "identity_species",
                    "v_gene",
                    "v_identity",
                    "j_gene",
                    "j_identity",
                ]
                numbering_ = [("%d%s" % (p)).strip() for p in positions]
                # print(",".join(fields), file=out)

                # Iterate over the domains identified
                for i, j in chain_types[cts]:
                    if details[i][j].get("germlines", {}).get("j_gene", [0])[0]:
                        j_gene = details[i][j].get("germlines", {}).get("j_gene", [["", ""], 0])[0][1]
                        j_gene_score = "%.2f" % details[i][j].get("germlines", {}).get("j_gene", [["", ""], 0])[1]
                    else:
                        j_gene = nan
                        j_gene_score = nan
                    line = [
                        sequences[i][0].replace(",", " "),
                        sequences[i][1],
                        str(j),
                        details[i][j].get("species", ""),
                        details[i][j].get("chain_type", ""),
                        str(details[i][j].get("evalue", "")),
                        str(details[i][j].get("bitscore", "")),
                        str(numbered[i][j][1]),
                        str(numbered[i][j][2]),
                        details[i][j].get("germlines", {}).get("v_gene", [["", ""], 0])[0][0],
                        details[i][j].get("germlines", {}).get("v_gene", [["", ""], 0])[0][1],
                        "%.2f" % details[i][j].get("germlines", {}).get("v_gene", [["", ""], 0])[1],
                        j_gene,
                        j_gene_score,
                    ]

                    # Hash the numbering. Insertion order has been preserved in the positions sort.
                    d = dict(numbered[i][j][0])
                    seq_ = [d.get(p, "-") for p in positions]
                    numbering_ = [[int(p[0]), p[1]] for p in positions]
                    assert len(line) == len(fields)
                    _df = pd.Series(line, index=fields)
                    _df["Chain"] = cts
                    _df["Numbering"] = list(map(lambda x: x[0], numbering_))
                    _df["Insertion"] = list(map(lambda x: x[1].strip(), numbering_))
                    _df["Numbered_Sequence"] = seq_
                    summary_dataframes.append(_df)

        return summary_dataframes
