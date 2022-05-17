import re

from sadie.hmmer.clients import G3


def test_g3_hmm(fixture_setup):
    g3 = G3()

    species_list = [
        "human",
        "mouse",
        "rat",
        "rabbit",
    ]
    chains = [
        "H",
        "L",
        "K",
    ]

    numbering_stockholm_folder = fixture_setup.numbering_data

    for species in species_list:
        for chain in chains:

            try:
                stockholm_pairs = g3.get_stockholm_pairs(species=species, chain=chain)
            except ValueError:
                print(f"G3 cannot find this species/chain, {species} {chain}")
                continue

            print(f"{species} {chain}")

            g3_name2align = {}
            numbering_name2align = {}

            for name, align in stockholm_pairs:
                name = "-".join([name_seg for name_seg in re.split(r"_|-", name) if "*" in name_seg])
                g3_name2align[name] = align

            if (numbering_stockholm_folder / f"{species}_{chain}.stockholm").is_file() is False:
                print(f"{species} {chain} numbering sto not found")
                continue

            with open(numbering_stockholm_folder / f"{species}_{chain}.stockholm") as f:
                for line in f.read().split("\n"):
                    if not line:
                        continue
                    if line[0] in ["#", "/"]:
                        continue
                    name, align = line.split()
                    name = "-".join([name_seg for name_seg in re.split(r"_|-", name) if "*" in name_seg])
                    numbering_name2align[name] = align

            for name, numbering_align in numbering_name2align.items():
                g3_align = g3_name2align.get(name)
                if not g3_align:
                    print("missing", name)  # TODO: Only looking at Ig for now
                    continue
                # seqs can change but the alignment positions from end of V to J looks highly conserved so we check there only.
                # TODO: fuzzy match?
                assert numbering_align[-5:] == g3_align[-5:]
                assert numbering_align[-35:].count("-") == g3_align[-35:].count("-")
