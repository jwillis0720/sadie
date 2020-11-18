from pybody import antibody
import tempfile


kappa_chain_nt = antibody.KappaChainNT(
    fwr1_nt="GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC",
    cdr1_nt="CACTCAATCTACCGGAAT",
    fwr2_nt="CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC",
    cdr2_nt="AACGCCAAC",
    fwr3_nt="TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC",
    cdr3_nt="CAACAGTACTATCGAGGATGGACG",
    fwr4_nt="TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG",
    v_gene="IGKV22S4",
    j_gene="IGKJ1",
    species="rat",
)

heavy_chain_nt = antibody.HeavyChainNT(
    fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
    cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
    fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
    cdr2_nt="ATTTGGTCCAACGGTGGCACC",
    fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
    cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
    fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
    v_gene="IGHV2-47",
    j_gene="IGHJ3",
    species="rat",
)

lambda_chain_aa = antibody.LambdaChainAA(
    name="fezakinumab",
    fwr1_aa="QAVLTQPPSVSGAPGQRVTISCTGS",
    cdr1_aa="SSNIGAGYG",
    fwr2_aa="VHWYQQLPGTAPKLLIY",
    cdr2_aa="GDS",
    fwr3_aa="NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC",
    cdr3_aa="QSYDNSLSGYV",
    fwr4_aa="FGGGTQLTVL",
    v_gene="IGLV1-40*01",
    j_gene="IGLJ7*01",
    species="human",
)
heavy_chain_aa = antibody.HeavyChainAA(
    name="Antibody2",
    fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
    cdr1_aa="GLSLTSNS",
    fwr2_aa="VSWIRQPPGKGLEWMGV",
    cdr2_aa="IWSNGGT",
    fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
    cdr3_aa="ASIYYYDADYLHWYFDF",
    fwr4_aa="WGPGTMVTVSS",
    v_gene="IGHV2-47",
    j_gene="IGHJ3",
    species="rat",
)


def test_construction():
    antibody_object_aa = antibody.AntibodyAA(heavy_chain_aa, lambda_chain_aa)
    assert (
        antibody_object_aa.get_segmented_vdj_aa()
        == "QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS\n\nQAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSGYV FGGGTQLTVL"
    )
    antibody_object_nt = antibody.AntibodyNT(heavy_chain_nt, kappa_chain_nt)
    assert (
        antibody_object_nt.get_segmented_alignment_nt()
        == "IGHV2-47|IGHJ3  CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGT\nheavy_chain     ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC\n\nIGHV2-47|IGHJ3  TATCATTAACCAGCAATAGT GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGA\nheavy_chain     ....GC.C........CTCC ..C........A..T.....G.....C..A.................T..G ..T...T\n\nIGHV2-47|IGHJ3  GTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTT\nheavy_chain     CC..C..T.....C ..C..C..C..C......G.GAG...CT..TCT....A.C.C........T..ATC...G.....\n\nIGHV2-47|IGHJ3  CTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA----------------------\nheavy_chain     ...G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTTATTACTATGACGCTGACTAC\n\nIGHV2-47|IGHJ3  ----ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG\nheavy_chain     CTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C\n\n\n\nIGKV22S4|IGKJ1  GACATCCAGATGACCCAGTCTCCTTCATTCCTGTCTGCATCTGTGGGAGACAGAGTCACTATCAACTGCAAAGCAAGT C\nkappa_chain     ........A.....A..T..G........G...AG...G........T...C.C....G.C.G........G..CTCC .\n\nIGKV22S4|IGKJ1  AGAATATTAACAGGTAC TTAAACTGGTACCAGCAAAAGCTTGGAGAAGCTCCCAAACTCCTGATATAT AATGCAAAC \nkappa_chain     .CTCA..CT..C..A.T C.GGC.........A..G..A..C..T..G.....A.....A..C..C..C ..C..C... \n\nIGKV22S4|IGKJ1  AGTTTGCAAACGGGCATCCCATCAAGGTTCAGTGGCAGTGGATCTGGTACTGATTTCACACTCACCATCAGCAGCCTGCA\nkappa_chain     TC.C....G..A..A.....G..T..A..T..C..ATCC..C..C.....C..C.....C..G.....T...TC......\n\nIGKV22S4|IGKJ1  GCCTGAAGATGTTGCCACATATTTCTGC TTGCAGCATAATAGTTGGCCGGTGGACG T--TCGGTGGAGGCACCAAGCT\nkappa_chain     ...C..G.....G..G..C......... CAA...T.CT.----.C.AG.A...... .GGA..T.C.GT.G.GGT.CAA\n\nIGKV22S4|IGKJ1  GGAATTGAAAC\nkappa_chain     A.CTGGAGCTG\n\n"
    )


def test_antibody_io():
    ##Antibody AA
    antibody_object_aa = antibody.AntibodyAA(heavy_chain_aa, lambda_chain_aa)
    reconstructed = antibody.AntibodyAA.from_json(antibody_object_aa.get_json())
    assert antibody_object_aa == reconstructed
    with tempfile.NamedTemporaryFile(suffix="json") as t_file:
        antibody_object_aa.to_json(t_file.name)
        second_class = antibody.AntibodyAA.read_json(t_file.name)
        assert antibody_object_aa == second_class

    antibody_object_nt = antibody.AntibodyNT(heavy_chain_nt, kappa_chain_nt)
    reconstructed = antibody.AntibodyNT.from_json(antibody_object_nt.get_json())
    assert antibody_object_nt == reconstructed
    with tempfile.NamedTemporaryFile(suffix="json") as t_file:
        antibody_object_nt.to_json(t_file.name)
        second_class = antibody.AntibodyNT.read_json(t_file.name)
        assert antibody_object_nt == second_class
