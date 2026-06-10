RECEPTORS = ["IG", "TR"]
IMGT_DEF_nt = {
    "FW1": {"start": 1, "end": 78},
    "CDR1": {"start": 79, "end": 114},
    "FW2": {"start": 115, "end": 165},
    "CDR2": {"start": 166, "end": 195},
    "FW3": {"start": 196, "end": 312},
    "CDR3": {"start": 312, "end": ""},
    "V-REGION": {"start": 0, "end": ""},
}

IMGT_GB_LOOKUP = {
    "Canis_lupus_familiaris_boxer": "dog",
    "Felis_catus_Abyssinian": "cat",
    "Canis_lupus_familiaris_Canis_lupus_familiaris_boxer": "dog",
    "Rattus_norvegicus_BN;_Sprague-Dawley": "rat",
    "Rattus_norvegicus_BN/SsNHsdMCW": "rat",
    "Rattus_norvegicus": "rat",
    "Mus_musculus_C57BL/6": "mouse",
    "Mus_musculus_BALB/c": "mouse",
    "Mus_musculus_C57BL/6J": "mouse",
    "Mus_musculus_MRL/lpr": "mouse",
    "Mus_musculus": "mouse",
    "Mus_musculus_A/J": "mouse",
    "Mus_musculus_C57BL/10": "mouse",
    "Mus_musculus_129/Sv": "mouse",
    "Mus_musculus_NZB": "mouse",
    "Mus_musculus_I/St": "mouse",
    "Mus_musculus_NFS": "mouse",
    "Mus_musculus_BALB.K": "mouse",
    "Mus_musculus_C3H": "mouse",
    "Mus_musculus_NZB/BINJ": "mouse",
    "Mus_musculus_CE/J": "mouse",
    "Mus_musculus_PERU": "mouse",
    "Mus_musculus_AKR": "mouse",
    "Mus_musculus_domesticus": "mouse",
    "Mus_musculus_O20/A": "mouse",
    "Mus_musculus_castaneus": "mouse",
    "Mus_musculus_molossinus_MOLF/Ei": "mouse",
    "Mus_musculus_musculus": "mouse",
    "Mus_musculus_castaneus_CAST/Ei": "mouse",
    "Mus_musculus_C58": "mouse",
    "Mus_musculus_SK": "mouse",
    "Mus_musculus_PERA": "mouse",
    "Mus_musculus_MRL": "mouse",
    "Mus_musculus_MBK": "mouse",
    "Mus_musculus_PWK": "mouse",
    "Mus_musculus_MAI": "mouse",
    "Homo_sapiens": "human",
    "Vicugna_pacos": "alpaca",
}


IMGT_LOOKUP = {
    "human": "Homo_sapiens",
    "cow": "Bos_taurus",
    "camel": "Camelus_dromedarius",
    "dog": "Canis_lupus_familiaris",
    "cat": "Felis_catus",
    "junglefowl": "Gallus_gallus",
    "night_monkey": "Aotus_nancymaae",
    "goat": "Capra_hircus",
    "sharks": "Chondrichthyes",
    "zebrafish": "Danio_rerio",
    "horse": "Equus_caballus",
    "cod": "Gadus_morhua",
    "catfish": "Ictalurus_punctatus",
    "crabmacaque": "Macaca_fascicularis",
    "macaque": "Macaca_mulatta",
    "mouse": "Mus_musculus",
    "ferret": "Mustela_putorius_furo",
    "nhp": "Nonhuman_primates",
    "trout": "Oncorhynchus_mykiss",
    "platypus": "Ornithorhynchus_anatinus",
    "rabbit": "Oryctolagus_cuniculus",
    "sheep": "Ovis_aries",
    "rat": "Rattus_norvegicus",
    "salmon": "Salmo_salar",
    "boar": "Sus_scrofa",
    "teleosts": "Teleostei",
    "dolphin": "Tursiops_truncatus",
    "alpaca": "Vicugna_pacos",
}

REVERSE_IMGT_LOOKUP = {v: k for k, v in IMGT_LOOKUP.items()}

BLAST_CONVENTION = {"IG": "Ig", "TR": "TCR"}


SEGMENTS_INTERNAL_DATA = {
    "TR": ["TRAV", "TRBV", "TRDV", "TRGV"],
    "IG": ["IGHV", "IGKV", "IGLV"],
}
SEGMENTS = {
    "TR": ["TRAV", "TRBD", "TRBJ", "TRBV", "TRDD", "TRDJ", "TRDV", "TRGV", "TRGJ"],
    "IG": ["IGHD", "IGHJ", "IGHV", "IGKJ", "IGKV", "IGLJ", "IGLV"],
}


# Per-species regular expressions for the conserved J-gene framework-4 (FWR4) anchor.
#
# Each pattern marks the start of FWR4 in a *translated* germline J segment: the
# IMGT invariant J-TRP 118 for the heavy locus (canonical ``WG.G`` == W-G-x-G) and
# J-PHE 118 for the light loci (canonical ``FG.G`` == F-G-x-G). The non-canonical
# variants -- more specific (``WGQG`` platypus, ``WG.GT`` rabbit), less specific
# (``W[GD].G`` horse), or sideways (``[WL]G[TQK][VG]`` alpaca) -- are NOT taken from
# any published source. They were curated empirically by inspecting the J alleles
# available at IMGT/OGRDB for each organism; the per-species ``ignore`` lists (gene
# ids such as ``IGKJ3`` and literal pseudogene peptides) are hand exceptions noted
# during that curation. They are therefore reliable on the in-database alleles they
# were fitted to, but not a dependable general caller for de-novo J alleles
# (over-specific patterns miss divergent alleles; loose ones can mis-anchor).
#
# Historically (from the 2020-11-17 pibody init commit 4d846853 through commit
# d81db6ae on 2021-08-21) these motifs WERE the active FWR4 caller: the now-removed
# ``genesegment.py`` matched each pattern against the translated J allele
# (``re.findall(motif, aa)[0]`` -- the FIRST match) to compute the CDR3-end /
# FWR4-start nucleotide index, and ``aux_file.py`` wrote that index into the IgBLAST
# ``*_gl.aux`` files. That logic moved to G3, leaving this table as a curation
# reference (it was an active caller, not orphaned at birth).
#
# These motifs are a curation reference and are NOT used by SADIE's current FWR4
# pipeline: FWR4 boundaries come from the IgBLAST auxiliary files
# (``data/germlines/aux_db/imgt/*_gl.aux``) and the numbering schemes / G3. The
# empirical accuracy of each pattern against the shipped germline J data is pinned
# by tests/unit/reference/test_motif_lookup.py. See GitHub issue #267 for context.
MOTIF_LOOKUP = {
    "mouse": {
        "IGHJ": r"WG.G",
        "IGKJ": r"F[SG].G",
        "IGLJ": r"FG.G",
        "TRAJ": r"[FLWC][GSVA].[GEWE]",
        "TRBJ": r"[FH][GA].G",
        "TRGJ": r"FA[EK]G",
        "ignore": ["IGLJ3P", "TRAJ41", "TRAJ60", "TRAJ61", "TRBJ2-6"],
    },
    "rat": {"IGHJ": r"WG.G", "IGKJ": r"FG.G", "IGLJ": r"[FL]G.G", "ignore": ["IGKJ3"]},
    "human": {
        "IGHJ": r"WG.G",
        # Anchored to the full F-G-x-G FR4 consensus (FGQG/FGPG/FGGG across all
        # shipped human IGKJ alleles). A bare ``FG`` matches the first Phe-Gly
        # dipeptide anywhere and can mis-anchor FR4 on de-novo alleles (#267).
        "IGKJ": r"FG.G",
        "IGLJ": r"FG.G",
        "TRAJ": r"[FWC][GA].[GEN]",
        "TRBJ": r"[FVG][GR].[G]",
        "TRGJ": r"F[GA].G",
        "TRDJ": r"FG.G",
        "ignore": [],
    },
    "macaque": {
        "IGHJ": r"WG.G",
        "IGKJ": r"FG.G",
        "IGLJ": r"F[GC].GT",
        "TRAJ": r"[FWS][GV].[GSRE][TMVS]",
        "TRBJ": r"[F][G].[G]",
        "TRDJ": r"FG.G",
        "TRGJ": r"F[GA].[G*]",
        "ignore": [],
    },
    "nhp": {
        "TRAJ": r"[F][G].[G][T]",
        "TRBJ": r"[F][G].[G][TS]",
        "TRDJ": r"F[GDA].[GT]T",
        "TRGJ": r"F[GDA].[GT]T",
        "ignore": ["TRAJ34", "TRBJ1-1", "TRBJ1-4", "TRBJ1-5", "TRBJ2-5"],
    },
    "platypus": {"IGHJ": r"WGQG", "ignore": []},
    "rabbit": {
        "IGHJ": r"WG.GT",
        "IGKJ": r"[FLR]G.[GE]T",
        "IGLJ": r"F[GS].[RG]T",
        "TRAJ": r"[FLW][GE].[CGRK][TMS]",
        "TRBJ": r"[F]G.G[TS]",
        "TRDJ": r"[F]G.G[TS]",
        "TRGJ": r"[F]G.GT",
        "ignore": [],
    },
    "night_monkey": {
        "TRAJ": r"[F][G].[G][T]",
        "TRBJ": r"[F][G].[G][TS]",
        "TRDJ": r"F[GA].[G]T",
        "TRGJ": r"F[GA].[G]T",
        "ignore": ["TRAJ34", "TRBJ1-1", "TRBJ1-4", "TRBJ1-5", "TRBJ2-5", "TRDJ2"],
    },
    "boar": {
        "ignore": [],
        "IGHJ": r"WG.G",
        "IGKJ": r"FG.GT",
        "IGLJ": r"FG.GT",
        "TRBJ": r"FG.G",
    },
    "cow": {
        "ignore": [],
        "IGHJ": r"[WC][SG][QPSR].",
        "IGKJ": r"[FL]G.[GR]T..E",
        "IGLJ": r"[FL][GI][SG][GR]T",
        "TRAJ": r"[FWL][GAS].[GK][TS]",
        "TRDJ": r"FG.[GE]",
        "TRGJ": r"[FLY][GN][VEK][GA]",
    },
    "crabmacaque": {"ignore": [], "IGHJ": r"WG.G"},
    "dolphin": {
        "ignore": [],
        "TRAJ": r"[FCLWSY][GS].[GRLK]",
        "TRGJ": r"[FCLWSY]G.[GRL]",
        "TRDJ": r"[FCLWSY][RG].[GRL]",
    },
    "ferret": {"ignore": [], "TRBJ": r"[F][GA].G", "TRAJ": r"[F][GA].G"},
    "camel": {
        "ignore": ["TRBJ3-4"],
        "IGHJ": r"WG.G",
        "IGKJ": r"[FL]G.GT",
        "IGLJ": r"FG.GT",
        "TRBJ": r"FG.G",
        "TRGJ": r"FG.G",
    },
    "goat": {"ignore": [], "IGKJ": r"[FL]G.GT", "IGLJ": r"[FL]G.GT"},
    "horse": {"ignore": [], "IGKJ": r"[F]G.GT", "IGHJ": r"[W][GD].G"},
    "dog": {
        "IGHJ": r"WG.G",
        "IGKJ": r"F[GS].G",
        "IGLJ": r"FG.G",
        "TRAJ": r"[FSW][GW].[GRLER]",
        "TRBJ": r"[F][GA].[G]",
        "TRGJ": r"[LFM][GTA].[GDV]",
        "TRDJ": r"FG.[GL]",
        "ignore": [],
    },
    "cat": {
        "IGHJ": r"WG.G",
        "IGKJ": r"F[G].G",
        "IGLJ": r"F[GNS].G",
        "TRAJ": r"[FWL][ERGW].[GRCKEQ]",
        "TRBJ": r"[F][TG].[G]",
        "TRGJ": r"[SF][TGAD].[G]",
        "TRDJ": r"FG.[G]",
        "ignore": [],
    },
    "alpaca": {"IGHJ": r"[WL]G[TQK][VG]", "ignore": []},
    "salmon": {"IGHJ": r"[W*][EG][KQ]GT", "ignore": []},
    "sharks": {
        "ignore": [
            "PEKGVGTVLTVR",
            "SYEYGGGTVVTVNP",
            "RHGLLGTRDHGDGDC",
            "ACGDGTFVTVNP",
            "YGADTVVTVNP",
            "YAACGAGTAVTVNP",
            "LSRLLGTRDHGDGDC",
            "YGSGTVLTVNP",
            "YGGGTVVTVNP",
            "HHGLLGTRDHGDGDF",
            "GLLGTRDHGDGDC",
            "SFDEYGGGTVVT",
            "SPNYWGGGSMVTVTC",
            "YAAVGDGTAVTVNP",
            "YAACGDATAVTVNP",
            "DYKGGDTLLTVK",
            "SYEYGGGTVVT",
            "CGDNTAVTVNP",
            "ERPGTALTVK",
            "QLCCMRRRHCRD",
            "HHGLLGTRDHGDGDC",
            "MLHAEMALRDCES",
            "YEKGAGTVLTVK",
            "NEKGAGTVLTVK",
            "DEEGAGTVLTVK",
            "GGAGTVLTVK",
            "YGGGTGVTVNP",
            "LPRLLGTRDHGDGDC",
        ],
        "IGHJ": r"[WCH][G].[RGS][TK]",
    },
    "sheep": {
        "IGLJ": r"[LF]G[GS]GT",
        "IGHJ": r"W[GD].[GR]",
        "IGKJ": r"FG[PGQ]GT",
        "TRAJ": r"[LFWE][GAKW].[GQ][TDRA]",
        "TRBJ": r"[F][G].[G][TS]",
        "TRDJ": r"FG.[EG]T",
        "ignore": ["FGDFYFLRGEGRRLAVV", "RPGAALTYGAGSGLAAG"],
    },
    "trout": {
        "ignore": ["SGAYAAYFGEXTKLTVL", "SYSEAYXXAGXKLTVL"],
        "IGHJ": r"WG.G",
        "TRBJ": r"FG.G[ATS]",
    },
    "zebrafish": {
        "ignore": ["IGIJ1", "IGIJ2", "IGIJ3", "IGIJ5", "IGIJ6S1", "IGIJ7S1", "IGIJ8S1"],
        "IGHJ": r"WG.GT",
        "TRAJ": r"[FM][GAST].G[TVSM]",
        "TRDJ": r"FG.P",
    },
}


# Machine-readable provenance / curation metadata for ``MOTIF_LOOKUP`` -- one entry
# per species key, recording where each set of motifs came from and how far it can
# be trusted. ``imgt_validated`` is True only where the motif is the canonical
# Lefranc IMGT-ONTOLOGY FR4 anchor end-to-end (W-G-x-G heavy / F-G-x-G light); every
# other species is an empirical fit to the IMGT/OGRDB J alleles shipped with SADIE,
# with no published per-species citation in repo history. See the comment above
# ``MOTIF_LOOKUP`` and GitHub issue #267.
_CANONICAL_FR4_SOURCE = (
    "Canonical IMGT FR4 anchor (Lefranc IMGT-ONTOLOGY): J-TRP 118 (W-G-x-G, heavy) " "/ J-PHE 118 (F-G-x-G, light)."
)
_EMPIRICAL_FR4_SOURCE = (
    "Empirical fit to the IMGT/OGRDB J alleles shipped with SADIE; introduced in the "
    "original pibody init commit 4d846853 (2020-11-17, J. Willis). No published "
    "per-species citation exists in repo history."
)
MOTIF_PROVENANCE = {
    "mouse": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "IGHJ/IGLJ are canonical W-G-x-G / F-G-x-G; IGKJ F[SG].G is a non-canonical empirical variant.",
    },
    "rat": {
        "source": _CANONICAL_FR4_SOURCE,
        "imgt_validated": True,
        "last_reviewed": "2026-06-10",
        "notes": "IGHJ W-G-x-G and IGKJ/IGLJ F-G-x-G match the canonical IMGT FR4 anchor.",
    },
    "human": {
        "source": _CANONICAL_FR4_SOURCE,
        "imgt_validated": True,
        "last_reviewed": "2026-06-10",
        "notes": "IGKJ corrected FG->FG.G (#267) to avoid anchoring on an upstream Phe-Gly.",
    },
    "macaque": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "IGL F[GC].GT misses the divergent IGLJ7*02 allele (see KNOWN_MOTIF_MISSES).",
    },
    "nhp": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "TR-only motifs fitted to the shipped non-human-primate J alleles.",
    },
    "platypus": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Over-specific WGQG; will miss divergent de-novo alleles.",
    },
    "rabbit": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Motifs carry the rabbit-specific ...T tail; fitted to the shipped J alleles.",
    },
    "night_monkey": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "TR-only motifs fitted to the shipped Aotus (night monkey) J alleles.",
    },
    "boar": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped pig J alleles; de-novo reliability untested.",
    },
    "cow": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Highly bovine-specific heavy motif [WC][SG][QPSR].; needs IMGT validation.",
    },
    "crabmacaque": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Heavy-only WG.G fitted to the shipped crab-eating macaque J alleles.",
    },
    "dolphin": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Permissive TR motifs fitted to the shipped dolphin J alleles.",
    },
    "ferret": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "TR-only motifs fitted to the shipped ferret J alleles.",
    },
    "camel": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped camel J alleles; de-novo reliability low.",
    },
    "goat": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Light-only motifs fitted to the shipped goat J alleles.",
    },
    "horse": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Non-canonical W[GD].G heavy motif; fitted to the shipped horse J alleles.",
    },
    "dog": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped dog J alleles; de-novo reliability untested.",
    },
    "cat": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped cat J alleles; de-novo reliability untested.",
    },
    "alpaca": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Over-specific [WL]G[TQK][VG]; de-novo reliability low.",
    },
    "salmon": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Over-specific [W*][EG][KQ]GT; de-novo reliability low.",
    },
    "sharks": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Highly permissive [WCH][G].[RGS][TK]; needs IMGT validation.",
    },
    "sheep": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped sheep J alleles; de-novo reliability untested.",
    },
    "trout": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped trout J alleles; de-novo reliability low.",
    },
    "zebrafish": {
        "source": _EMPIRICAL_FR4_SOURCE,
        "imgt_validated": False,
        "last_reviewed": "2026-06-10",
        "notes": "Fitted to the shipped zebrafish J alleles; de-novo reliability low.",
    },
}

J_SEGMENTS = {"IG": ["IGHJ", "IGKJ", "IGKJ"], "TR": ["TRAJ", "TRBJ", "TRGJ", "TRDJ"]}
RENAME_DICT = {
    "CDR1-IMGT": "cdr1_nt",
    "FR1-IMGT": "fwr1_nt",
    "FR2-IMGT": "fwr2_nt",
    "FR3-IMGT": "fwr3_nt",
    "CDR2-IMGT": "cdr2_nt",
    "CDR3-IMGT": "cdr3_nt",
    "FW1": "fwr1_nt",
    "FW2": "fwr2_nt",
    "FW3": "fwr3_nt",
    "CDR1": "cdr1_nt",
    "CDR2": "cdr2_nt",
    "CDR3": "cdr3_nt",
    "V-REGION": "v_gene_nt",
    "D-REGION": "d_gene_nt",
    "J-REGION": "j_gene_nt",
}
RENAME_DICT_TRANSLATE = {
    "cdr1_nt": "cdr1_aa",
    "fwr1_nt": "fwr1_aa",
    "fwr2_nt": "fwr2_aa",
    "fwr3_nt": "fwr3_aa",
    "cdr2_nt": "cdr2_aa",
    "cdr3_nt": "cdr3_aa",
}
