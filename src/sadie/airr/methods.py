#     if self.infer:
#         self._table["vdj_igl"] = self._table.apply(self._get_igl, axis=1)
#         self._table.loc[self._table[self._table["vdj_igl"].isna()].index, "note"] = "liable"
#         self._table["igl_mut_aa"] = self._table[["vdj_aa", "vdj_igl"]].apply(self._get_diff, axis=1)

# def _get_igl(self, row: pd.Series) -> str:
#     """Get infered germline sequenxe

#     Parameters
#     ----------
#     row : pd.Series
#         A row from the airr table
#     Returns
#     -------
#     str
#         the igl sequecne
#     """

#     # get germline components
#     v_germline = row.v_germline_alignment_aa
#     full_germline = row.germline_alignment_aa
#     if isinstance(v_germline, float):
#         return
#     cdr3_j_germline = full_germline[len(v_germline) :]

#     # get mature components
#     v_mature = row.v_sequence_alignment_aa
#     full_mature = row.sequence_alignment_aa
#     cdr3_j_mature = full_mature[len(v_mature) :]

#     # if the mature and cdr3 are not the same size
#     # this will happen on non-productive
#     if len(cdr3_j_mature) != len(cdr3_j_germline):
#         logger.debug(f"{row.name} - strange iGL")
#         return

#         # # quick aligment
#         # _alignments = align.globalxs(cdr3_j_mature, cdr3_j_germline, -10, -1)
#         # cdr3_j_mature, cdr3_j_germline = _alignments[0][0], _alignments[0][1]

#     iGL_cdr3 = ""
#     for mature, germline in zip(cdr3_j_mature, cdr3_j_germline):
#         if germline == "X" or germline == "-":
#             iGL_cdr3 += mature
#             continue
#         iGL_cdr3 += germline

#     full_igl = v_germline.replace("-", "") + iGL_cdr3.replace("-", "")
#     return full_igl
