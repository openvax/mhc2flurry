"""
Select allele sequences for pan-class II models by analyzing distances between
each MHC residue and the peptide across a set of structures from PDB.
"""
from __future__ import print_function

import sys
import argparse
import collections
import os
import operator

import numpy
import pandas

import tqdm

import atomium


from mhc2flurry.fasta import read_fasta_to_dataframe


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "alpha_aligned_fasta",
    metavar="FASTA",
    help="Aligned sequences")

parser.add_argument(
    "beta_aligned_fasta",
    metavar="FASTA",
    help="Aligned sequences")

parser.add_argument(
    "pdb_dir",
    metavar="DIR",
    help="Directory containing PDB structures")

parser.add_argument(
    "--reference-allele",
    nargs=2,
    help="Alpha and beta alleles to use for position numbering.")

parser.add_argument(
    "--reference-structure",
    action="append",
    default=[],
    help="Structures to write out with b-factors rewritten according to "
    "inclusion in pseudosequences(for visualization).")

parser.add_argument(
    "--out-csv",
    help="Result file for sequences")

parser.add_argument(
    "--out-aux-dir",
    help="Result DIR for extra information")

parser.add_argument(
    "--cutoffs",
    default=[2.0, 4.0, 6.0, 8.0, 10.0],
    nargs="+",
    type=float,
    metavar="X",
    help="Cutoff distances to evaluate. Default: %(default)s.")

parser.add_argument(
    "--criteria",
    nargs=3,
    type=float,
    action="append",
    default=[],
    required=True,
    metavar="X",
    help="Criteria for selecting a position. Triple of: min minor allele "
    "fraction, cutoff distance, fraction of structures with a contact at "
    "the given cutoff. May be specified any number of times.")

parser.add_argument(
    "--peptide-chain-min-length",
    default=5,
    metavar="N",
    type=int,
    help="Default: %(default)s.")

parser.add_argument(
    "--peptide-chain-max-length",
    default=50,
    metavar="N",
    type=int,
    help="Default: %(default)s.")

parser.add_argument(
    "--subsample-pdb",
    metavar="N",
    type=int,
    help="Subsample to at most N PDB structures. For debugging.")


def make_position_to_aligned_position_dict(aligned_sequence):
    result = {}
    position = 0
    for (i, char) in enumerate(aligned_sequence):
        if char != "-":
            result[position] = i
            position += 1
    return result


def make_aligned_position_to_position_dict(aligned_sequence):
    result = {}
    position = 0
    for (i, char) in enumerate(aligned_sequence):
        if char != "-":
            result[i] = position
            position += 1
    return result


def run():
    args = parser.parse_args(sys.argv[1:])
    print(args)

    alpha_aligned_df = read_fasta_to_dataframe(
        args.alpha_aligned_fasta,
        full_descriptions=True)
    alpha_aligned_df["kind"] = "alpha"

    beta_aligned_df = read_fasta_to_dataframe(
        args.beta_aligned_fasta,
        full_descriptions=True)
    beta_aligned_df["kind"] = "beta"

    aligned_df = pandas.concat(
        [alpha_aligned_df, beta_aligned_df], ignore_index=True)
    aligned_df["unaligned"] = aligned_df.sequence.str.replace("-", "")
    aligned_df = aligned_df.rename(columns={
        "sequence": "aligned_sequence",
    }).set_index("sequence_id")

    non_pdb_aligned_df = aligned_df.loc[
        ~aligned_df.index.str.startswith("pdb")
    ].copy()
    minor_allele_fraction_df = []
    for kind, sub_df in non_pdb_aligned_df.groupby("kind"):
        print("Calculating minor allelic fractions: ", kind)
        (length,) = sub_df.aligned_sequence.str.len().unique()
        for pos in tqdm.tqdm(range(length)):
            s = sub_df.aligned_sequence.str.get(pos)
            mode = s.mode()[0]
            maf = (s != mode).mean()
            minor_allele_fraction_df.append((kind, pos, mode, maf))
    minor_allele_fraction_df = pandas.DataFrame(
        minor_allele_fraction_df,
        columns=[
            "mhc_chain_kind",
            "mhc_residue_aligned",
            "major_allele",
            "minor_allele_fraction",
        ])
    minor_allele_fraction_df = minor_allele_fraction_df.set_index(
        ["mhc_chain_kind", "mhc_residue_aligned"])
    print(minor_allele_fraction_df)

    pdb_aligned_df = aligned_df.loc[
        aligned_df.index.str.startswith("pdb")
    ].copy()
    pdb_aligned_df["accession"] = pdb_aligned_df.index.str.split(".").str.get(
        1).str.split("_").str.get(0)
    pdb_aligned_df["chain"] = pdb_aligned_df.index.str.split("_").str.get(-1)

    if args.subsample_pdb:
        keep_accessions = list(
            pandas.Series(
                pdb_aligned_df.accession.unique()).sample(
                n=args.subsample_pdb)) + args.reference_structure
        pdb_aligned_df = pdb_aligned_df.loc[
            pdb_aligned_df.accession.isin(keep_accessions)
        ].copy()

    info_by_accession = {}
    contacts_df = []
    for accession, sub_df in tqdm.tqdm(
            pdb_aligned_df.groupby("accession"),
            total=pdb_aligned_df.accession.nunique()):

        sub_df = sub_df.set_index("chain")
        alpha_chains = sub_df.loc[sub_df.kind == "alpha"].index.values
        beta_chains = sub_df.loc[sub_df.kind == "beta"].index.values

        mhc_chain_to_kind = {}
        for chain in alpha_chains:
            mhc_chain_to_kind[chain] = "alpha"
        for chain in beta_chains:
            mhc_chain_to_kind[chain] = "beta"

        if len(alpha_chains) != len(beta_chains):
            print(
                "Skipping", accession,
                "because num chains for alpha != beta",
                len(alpha_chains), len(beta_chains))
            continue

        structure = atomium.open(
            os.path.join(
                args.pdb_dir, "%s.cif.gz" % accession)).model

        peptides = [
            c for c in structure.chains()
            if len(c) >= args.peptide_chain_min_length and
               len(c) <= args.peptide_chain_max_length
        ]

        if len(peptides) == 0:
            print("Skipping", accession, "because no peptides")
            continue

        structure.optimise_distances()

        if accession in args.reference_structure:
            # Save for later
            info_by_accession[accession] = {
                "structure": structure,
                "peptides": peptides,
                "mhc_chain_to_kind": mhc_chain_to_kind,
                "aligned_df": sub_df.copy(),
            }

        mhc_chain_to_position_map = {}
        for chain in mhc_chain_to_kind:
            mhc_chain_to_position_map[chain] = make_position_to_aligned_position_dict(
                sub_df.loc[chain, "aligned_sequence"])

        for peptide in peptides:
            seen = set()
            for cutoff in sorted(args.cutoffs):
                nearby = [
                    r for r in peptide.nearby_hets(
                        cutoff=cutoff, residues=True, ligands=False)
                    if r not in seen
                ]
                seen.update(nearby)
                for residue in nearby:
                    kind = mhc_chain_to_kind.get(residue.chain.id)
                    if kind is not None:
                        index = residue.chain.residues().index(residue)
                        row = sub_df.loc[residue.chain.id]
                        numpy.testing.assert_equal(
                            residue.code,
                            row.unaligned[index])
                        aligned_position = (
                            mhc_chain_to_position_map[residue.chain.id][index])
                        numpy.testing.assert_equal(
                            residue.code,
                            row.aligned_sequence[aligned_position])
                        contacts_df.append((
                            accession,
                            cutoff,
                            peptide.id,
                            residue.chain.id,
                            kind,
                            index,
                            aligned_position,
                            residue.code))

    contacts_df = pandas.DataFrame(
        contacts_df, columns=[
            "accession",
            "cutoff",
            "peptide_chain",
            "mhc_chain",
            "mhc_chain_kind",
            "mhc_residue_unaligned",
            "mhc_residue_aligned",
            "mhc_residue",
        ])

    num_accessions = contacts_df.accession.nunique()
    positional_contact_rates_df = contacts_df.groupby(
        ["mhc_chain_kind", "mhc_residue_aligned", "cutoff"]
    ).accession.nunique().unstack().reindex(
        sorted(args.cutoffs), axis=1).fillna(0.0).cumsum(1) / num_accessions

    positional_df = minor_allele_fraction_df.merge(
        positional_contact_rates_df,
        how="left",
        left_index=True,
        right_index=True).fillna(0)

    # Criteria name -> alpha or beta -> list of positions
    criteria_to_positions = collections.OrderedDict()
    for (maf, cutoff, fraction) in args.criteria:
        name = "maf_%s_and_%s_within_%s_angstrom" % (maf, fraction, cutoff)
        positional_df[name] = (
            (positional_df.minor_allele_fraction >= maf) &
            (positional_df[cutoff] >= fraction)
        )
        positions = positional_df.loc[
            positional_df[name]
        ].index.to_frame().reset_index(drop=True).groupby(
            "mhc_chain_kind"
        ).mhc_residue_aligned.unique().map(sorted).to_dict()

        criteria_to_positions[name] = positions
        print("Criteria", name, "selected:")
        for (k, v) in criteria_to_positions[name].items():
            print(k, len(v))

    pseudosequences_df = non_pdb_aligned_df.copy()
    for (criteria, d) in criteria_to_positions.items():
        for kind in ["alpha", "beta"]:
            positions = d.get(kind, [])
            sub = pseudosequences_df.loc[
                pseudosequences_df.kind == kind,
            ]
            pseudosequences_df.loc[
                sub.index,
                criteria
            ] = sub.aligned_sequence.map(
                operator.itemgetter(*positions)
            ).map("".join).str.replace("-", "X")
    pseudosequences_df.index = pseudosequences_df.index.str.split().str.get(1)
    assert pseudosequences_df.index.value_counts().max() == 1

    main_result_df = pseudosequences_df[
        list(criteria_to_positions) + ["kind"]
    ].copy()

    main_result_df.to_csv(args.out_csv, index=True)
    print("Wrote %s: " % str(main_result_df.shape), args.out_csv)

    if args.out_aux_dir:
        if not os.path.exists(args.out_aux_dir):
            os.mkdir(args.out_aux_dir)

        filename = os.path.join(args.out_aux_dir, "aligned_sequences.csv")
        pseudosequences_df.to_csv(filename, index=True)
        print("Wrote: ", filename)

        filename = os.path.join(args.out_aux_dir, "contacts.csv")
        contacts_df.to_csv(filename, index=True)
        print("Wrote: ", filename)

        # Positional. We add reference allele position numbering and amino acids.
        if args.reference_allele:
            write_df = positional_df.copy()
            (alpha_reference, beta_reference) = args.reference_allele
            reference_name = "%s/%s" % (alpha_reference, beta_reference)
            reference_alleles = {
                "alpha": alpha_reference,
                "beta": beta_reference,
            }
            for kind in ["alpha", "beta"]:
                reference_allele = reference_alleles[kind]
                reference_sequence = pseudosequences_df.loc[
                    reference_allele, "aligned_sequence"
                ]
                position_map = make_aligned_position_to_position_dict(
                    reference_sequence)
                write_df.loc[
                    kind,
                    reference_name + " position"
                ] = write_df.loc[
                    kind
                ].index.map(position_map)
                write_df.loc[
                    kind,
                    reference_name + " aa"
                ] = write_df.loc[
                    kind
                ].index.map(lambda pos: reference_sequence[pos])
        filename = os.path.join(args.out_aux_dir, "positional.csv")
        write_df.to_csv(filename, index=True)
        print("Wrote: ", filename)

        # Reference structures
        # Write out reference structures with the "bvalue" atom property used
        # to indicate minor allele fractions / fraction of residues within a
        # given distance of the peptide / inclusion in pseudosequences.
        # This can be used to generate colored renderings showing these
        # properties, e.g. in pymol.
        # This "b-factor" hack is commonly used to store arbitrary user data
        # in a PDB file. There may be a better way for CIF files but I don't
        # know of one.
        for accession in args.reference_structure:
            positional_with_residues_df = positional_df.copy()
            positional_with_residues_df[
                "residues"
            ] = positional_with_residues_df.index.map(lambda i: [])

            info = info_by_accession.get(accession)
            if not info:
                print("No info for reference structure", accession)
                continue

            structure = info['structure']

            for chain, row in info['aligned_df'].iterrows():
                position_map = make_position_to_aligned_position_dict(
                    row.aligned_sequence)
                residues_df = pandas.DataFrame({
                    "residue": structure.chain(chain).residues(),
                })
                residues_df["aligned_position"] = residues_df.index.map(
                    position_map)

                for _, residue_row in residues_df.iterrows():
                    positional_with_residues_df.loc[
                        (row.kind, residue_row.aligned_position),
                        "residues"
                    ].append(residue_row.residue)

            positional_with_residues_df = positional_with_residues_df.loc[
                positional_with_residues_df.residues.str.len() > 0
            ]

            quantitative_columns = positional_with_residues_df.dtypes.loc[
                (positional_with_residues_df.dtypes == float) |
                (positional_with_residues_df.dtypes == bool)
            ].index

            for atom in structure.atoms():
                atom.bvalue = 0

            for col in quantitative_columns:
                # Assign bfactors based on the particular column.
                for _, row in positional_with_residues_df.iterrows():
                    for residue in row.residues:
                        for atom in residue.atoms():
                            atom.bvalue = float(row[col]) * 100.0

                # Write out the file with modified bvalues.
                filename = os.path.join(
                    args.out_aux_dir,
                    "%s.%s.cif" % (accession, col))
                structure.save(filename)
                print("Wrote:", filename)



if __name__ == '__main__':
    run()
